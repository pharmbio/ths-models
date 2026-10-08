'''The ADMET-AI runtime: loading the model and predicting a list of SMILES.

Wraps admet_ai.ADMETModel, the package's own predictor; every value is the
one it computes. What this adds is what an API needs and ADMETModel does not
do:

- One row per input SMILES, in input order. ADMETModel drops a SMILES RDKit
  cannot parse and indexes the rest by SMILES string, so its rows no longer
  line up with the input; and given no parsable SMILES at all it raises.
- Torch threads sized to the container's CPU limit, not the host's cores.
- No Lightning progress bars or banners in the server log on every call.

server.py serves it; validate.py checks it against the predictions the
ADMET-AI authors ship for DrugBank approved drugs.
'''

import logging
import math
import os
import threading
import types
from typing import List

# ADMETModel wraps each step in a tqdm bar. tqdm reads TQDM_* defaults when it
# is first imported, which happens below, so this has to come first.
os.environ.setdefault("TQDM_DISABLE", "1")

import numpy as np
import torch
from lightning import pytorch as pl
from rdkit import Chem, RDLogger

import admet_ai.admet_model
from admet_ai import ADMETModel, __version__ as ADMET_AI_VERSION
from admet_ai.admet_info import get_admet_info

# Unparsable SMILES are reported in the response; RDKit's own lines would only
# fill the log. Lightning announces the (absent) GPU on every Trainer.
RDLogger.DisableLog("rdApp.*")
logging.getLogger("lightning.pytorch").setLevel(logging.WARNING)
logging.getLogger("lightning.fabric").setLevel(logging.WARNING)
# Chemprop warns it is "dropping" the last molecule when a batch is 1 more
# than a multiple of 64. Lightning sets drop_last=False when predicting, so
# nothing is dropped (and Predictor.predict checks the row count).
logging.getLogger("chemprop.data.dataloader").setLevel(logging.ERROR)


def _cpu_limit() -> int:
    '''CPUs this process may use: the cgroup quota if there is one, else the CPU affinity.

    torch otherwise starts a thread per host core; in a container limited to
    2 vCPU on a 64-core node, those threads are throttled and a prediction
    runs many times slower.
    '''

    try:
        with open("/sys/fs/cgroup/cpu.max") as fh:              # cgroup v2
            quota, period = fh.read().split()[:2]
        if quota != "max":
            return max(1, math.ceil(int(quota) / int(period)))
    except (OSError, ValueError):
        pass
    try:
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as fh:  # cgroup v1
            quota = int(fh.read())
        with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as fh:
            period = int(fh.read())
        if quota > 0:
            return max(1, math.ceil(quota / period))
    except (OSError, ValueError):
        pass
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:                                       # macOS
        return os.cpu_count() or 1


class _QuietTrainer(pl.Trainer):
    '''The Trainer ADMETModel builds for every ensemble on every call, without its progress bar.

    ADMETModel passes enable_progress_bar=True, which an environment variable
    cannot override; without this, each request prints ten progress bars.
    Only the bar and the model summary change, not the prediction.
    '''

    def __init__(self, *args, **kwargs):
        kwargs.update(enable_progress_bar=False, enable_model_summary=False)
        super().__init__(*args, **kwargs)


# ADMETModel looks up `pl.Trainer` when it predicts; give it the quiet one.
# Scoped to admet_ai's module, so Lightning itself is left as it is.
admet_ai.admet_model.pl = types.SimpleNamespace(Trainer=_QuietTrainer)


# Property metadata

def property_catalogue() -> List[dict]:
    '''What each predicted property is, from admet.csv in the admet_ai package.

    Returns:
    ----------
    properties (list): one dict per property, in the package's order, with
        its id (the column name), name, category, task type, units, species,
        training-set size and the test-set metrics the authors report.
    '''

    def text(x):
        # admet.csv writes "-", or leaves the cell empty, where nothing applies.
        return None if (isinstance(x, float) and not math.isfinite(x)) or x == "-" else x

    table = get_admet_info()
    return [{
        "id": row["id"],
        "name": row["name"],
        "category": row["category"],
        "task_type": row["task_type"],
        "units": text(row["units"]),
        "species": text(row["species"]),
        "training_size": int(row["size"]) if math.isfinite(row["size"]) else None,
        "test_metrics": {name: float(row[name]) for name in ("AUROC", "AUPRC", "R^2", "MAE")
                         if math.isfinite(row[name])} or None,
        "url": text(row["url"]),
    } for _, row in table.iterrows()]


# The model

class Predictor:
    '''ADMET-AI behind one lock, returning one row per input SMILES.'''

    def __init__(self, num_threads: int = 0):
        '''Load the ADMET-AI model and fix the output columns.

        Parameters:
        ---------
        num_threads (int): torch threads; 0 sizes them to the CPU limit.
        '''

        self.num_threads = num_threads or _cpu_limit()
        torch.set_num_threads(self.num_threads)
        self.version = ADMET_AI_VERSION
        # Without a DrugBank reference set, ADMETModel returns the predicted
        # values only, not the <property>_drugbank_approved_percentile columns.
        self.model = ADMETModel(drugbank_path=None)
        # Concurrent calls would share the ten Chemprop models, which each
        # Trainer attaches itself to; one prediction at a time.
        self._lock = threading.Lock()
        # Predicting one molecule fixes the column order, so a batch where
        # nothing parses still gets every column, and proves the model runs
        # before the server opens its port.
        self.columns = list(self.model.predict(smiles=["CCO"]).columns)

    def predict(self, smiles_list: List[str]):
        '''Predict a list of SMILES, keeping one row per input.

        A SMILES RDKit cannot parse, or that parses to no atoms, is not sent
        to the model: its row holds None in every column, and it is listed in
        the second return value.

        Parameters:
        ---------
        smiles_list (list): SMILES strings, in the order the rows should come back.

        Returns:
        ----------
        result (tuple): the rows, as lists of values in `self.columns` order, and the SMILES that did not parse.
        '''

        valid, unparsed = [], []
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None and mol.GetNumAtoms() > 0:
                valid.append(i)
            else:
                unparsed.append(smiles)

        rows = [[None] * len(self.columns) for _ in smiles_list]
        if valid:
            with self._lock:
                frame = self.model.predict(smiles=[smiles_list[i] for i in valid])
            # ADMETModel keeps every SMILES that parses, in order, so its rows
            # map back by position (its index is the SMILES, duplicates and all).
            if len(frame) != len(valid):
                raise RuntimeError("ADMET-AI returned {} rows for {} parsable SMILES.".format(
                    len(frame), len(valid)))
            frame = frame[self.columns]
            values = [_column_values(frame[name]) for name in self.columns]
            for position, i in enumerate(valid):
                rows[i] = [column[position] for column in values]
        return rows, unparsed


def _column_values(column) -> list:
    '''One output column as JSON-ready Python numbers.

    The Chemprop outputs are float32. Each is written with the fewest digits
    that read back as the same float32 (0.08982882, not 0.08982881903648376),
    as the DrugBank reference and a CSV of the same frame show them. Counts
    stay integers. NaN and infinity, which JSON cannot carry, become None.
    '''

    if column.dtype == np.float32:
        return [float(str(x)) if np.isfinite(x) else None for x in column.to_numpy()]
    if np.issubdtype(column.dtype, np.integer):
        return [int(x) for x in column.to_numpy()]
    return [float(x) if math.isfinite(x) else None for x in column.to_numpy(dtype=float)]
