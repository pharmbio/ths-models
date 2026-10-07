'''The RiskMix conformal QSAR runtime: descriptors, reading a model pickle, and
the conformal prediction itself.

Ported from chemsafe-agent's skills/qsar_modelling/scripts/utils.py, where it
was validated against the original 2021 stack (validate.py repeats that check
for this model). The arithmetic is unchanged; only the agent's input handling
and output files were left behind, and server.py replaces them.
'''

import pickle
from typing import List

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, GraphDescriptors, rdMolDescriptors


class PredictionError(RuntimeError):
    '''A prediction could not be produced, with a reason worth showing the caller.'''


# Descriptors
# The 119 RDKit descriptors the models were trained on, in the order the trees
# index them. Ported verbatim from the upstream descriptor.py -- the order and
# the exact set are part of the models, so neither may be changed.
_DESCRIPTOR_NAMES = [
    'MolLogP', 'MolMR', 'LabuteASA', 'TPSA', 'MolWt', 'ExactMolWt', 'CalcNumLipinskiHBA',
    'CalcNumLipinskiHBD', 'NumRotatableBonds',
    'CalcNumHBD', 'CalcNumHBA', 'CalcNumAmideBonds', 'NumHeteroatoms', 'HeavyAtomCount', 'NumAtoms',
    'CalcNumAtomStereoCenters',
    'CalcNumUnspecifiedAtomStereoCenters', 'RingCount', 'NumAromaticRings', 'NumSaturatedRings',
    'NumAliphaticRings',
    'NumAromaticHeterocycles', 'NumSaturatedHeterocycles', 'NumAliphaticHeterocycles',
    'NumAromaticCarbocycles',
    'NumSaturatedCarbocycles', 'NumAliphaticCarbocycles', 'FractionCSP3',
    'Chi0v', 'Chi1v', 'Chi2v', 'Chi3v', 'Chi4v', 'Chi1n', 'Chi2n', 'Chi3n', 'Chi4n',
    'HallKierAlpha', 'Kappa1', 'Kappa2', 'Kappa3',
    'SlogP_VSA1', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7',
    'SlogP_VSA8', 'SlogP_VSA9', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12',
    'SMR_VSA1', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8',
    'SMR_VSA9', 'SMR_VSA10',
    'PEOE_VSA1', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8',
    'PEOE_VSA9', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14',
    'MQNs_',
]


def compound_descriptors(smiles: str) -> np.ndarray:
    '''The 119 model descriptors for one SMILES, or raise if RDKit cannot parse it.

    Parameters:
    ---------
    smiles (str): one SMILES string.

    Returns:
    ----------
    descriptors (numpy array): the 119 descriptor values, in model order.
    '''

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise PredictionError("bad molecule {}!".format(smiles))
    Chem.Kekulize(mol)

    values = []
    for feature in _DESCRIPTOR_NAMES:
        if feature == 'NumAtoms':
            found = len(Chem.rdchem.Mol.GetAtoms(Chem.rdmolops.AddHs(mol)))
        else:
            found = next((getattr(source, feature)(mol)
                          for source in (Descriptors, rdMolDescriptors, GraphDescriptors)
                          if hasattr(source, feature)), None)
        if isinstance(found, list):
            values.extend(found)      # MQNs_ expands to 42 values
        else:
            values.append(found)
    return np.asarray(values)


def descriptor_matrix(smiles_list: List[str]):
    '''Descriptors for a list of SMILES, reporting which ones RDKit rejected.

    Parameters:
    ---------
    smiles_list (list): SMILES strings to featurize.

    Returns:
    ----------
    result (tuple): the (n_ok, 119) descriptor matrix, the indices of the molecules that failed, and those molecules.
    '''

    rows, bad_index, bad_smiles = [], [], []
    for i, smiles in enumerate(smiles_list):
        try:
            rows.append(compound_descriptors(smiles))
        except Exception:
            bad_index.append(i)
            bad_smiles.append(smiles)
    matrix = np.asarray(rows) if rows else np.empty((0, 119))
    return matrix, bad_index, bad_smiles


# Reading the model pickles

_REAL_ROOTS = {"numpy", "builtins", "collections", "functools", "copyreg", "_codecs"}


class _Stub:
    '''Placeholder for any class we refuse to import; keeps the pickled state.'''

    def __init__(self, *args, **kwargs):
        self._ctor_args = args

    def __setstate__(self, state):
        if isinstance(state, dict):
            self.__dict__.update(state)
        else:
            self.__dict__["_state"] = state


class _AnyCallable:
    '''Absorbs cloudpickle's function and code-object machinery.

    The pickles embed `condition = lambda x: x[1]` and
    `agg_func = lambda x: np.median(x, axis=2)` from Consensus.py as raw CPython
    3.5 bytecode. Both are reimplemented below, so the stored copies -- which
    could not be built or run anyway -- are discarded.
    '''

    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        return _AnyCallable()

    def __setstate__(self, state):
        pass


class _StubUnpickler(pickle.Unpickler):
    _classes = {}

    def find_class(self, module, name):
        root = module.split(".")[0]
        if root in _REAL_ROOTS:
            return super().find_class(module, name)
        if root == "cloudpickle":
            return _AnyCallable
        key = (module, name)
        if key not in self._classes:
            self._classes[key] = type(name, (_Stub,), {"__module__": module})
        return self._classes[key]


# The conformal predictor

class _Forest:
    '''A random forest flattened into padded arrays, one row per tree.'''

    __slots__ = ("left", "right", "feature", "threshold", "value")

    def __init__(self, trees):
        n_trees = len(trees)
        width = max(tree.node_count for tree in trees)
        n_classes = trees[0].values.shape[2]

        self.left = np.zeros((n_trees, width), dtype=np.int32)
        self.right = np.zeros((n_trees, width), dtype=np.int32)
        self.feature = np.zeros((n_trees, width), dtype=np.int32)
        self.threshold = np.zeros((n_trees, width), dtype=np.float64)
        self.value = np.zeros((n_trees, width, n_classes), dtype=np.float64)

        for i, tree in enumerate(trees):
            n = tree.node_count
            nodes = tree.nodes[:n]
            self.left[i, :n] = nodes["left_child"]
            self.right[i, :n] = nodes["right_child"]
            # Leaves keep feature == -2; clamp so it stays a valid column index.
            self.feature[i, :n] = np.maximum(nodes["feature"], 0)
            self.threshold[i, :n] = nodes["threshold"]

            # DecisionTreeClassifier.predict_proba normalises each leaf.
            values = tree.values[:n, 0, :].astype(np.float64)
            totals = values.sum(axis=1)
            totals[totals == 0.0] = 1.0
            self.value[i, :n] = values / totals[:, None]

    def predict_proba(self, x):
        '''Mean of the per-tree normalised leaf distributions.'''
        n_trees, n_samples = self.left.shape[0], x.shape[0]
        tree_idx = np.arange(n_trees)[:, None]
        sample_idx = np.arange(n_samples)[None, :]
        node = np.zeros((n_trees, n_samples), dtype=np.int32)

        while True:
            left = self.left[tree_idx, node]
            internal = left != -1                 # -1 marks a leaf (TREE_LEAF)
            if not internal.any():
                break
            feature = self.feature[tree_idx, node]
            go_left = x[sample_idx, feature] <= self.threshold[tree_idx, node]
            nxt = np.where(go_left, left, self.right[tree_idx, node])
            node = np.where(internal, nxt, node).astype(np.int32)

        return self.value[tree_idx, node].mean(axis=0)


class _Icp:
    '''One inductive conformal predictor: a forest plus its calibration scores.'''

    __slots__ = ("forest", "cal_scores", "classes")

    def __init__(self, forest, cal_scores, classes):
        self.forest = forest
        # Upstream stores them descending and reverses on every lookup.
        self.cal_scores = {k: np.asarray(v)[::-1].copy() for k, v in cal_scores.items()}
        self.classes = classes

    def p_values(self, x, smoothing, rng):
        proba = self.forest.predict_proba(x)
        p = np.zeros((x.shape[0], self.classes.size))

        for i, c in enumerate(self.classes):
            # InverseProbabilityErrFunc, in float32 as upstream computes it.
            prob = np.zeros(x.shape[0], dtype=np.float32)
            if int(c) < proba.shape[1]:
                prob[:] = proba[:, int(c)]
            nc = 1 - prob

            cal = self.cal_scores[c]              # label-conditional: condition = x[1]
            n_cal = cal.size
            idx_left = np.searchsorted(cal, nc, "left")
            idx_right = np.searchsorted(cal, nc, "right")
            n_gt = n_cal - idx_right
            n_eq = idx_right - idx_left + 1

            p[:, i] = n_gt / (n_cal + 1)
            if smoothing:
                # Upstream draws one uniform per sample, in sample order; one
                # vectorised call consumes the same MT19937 stream in the same
                # order, so seeded runs stay bit-identical to the original.
                p[:, i] += (n_eq * rng.uniform(0, 1, x.shape[0])) / (n_cal + 1)
            else:
                p[:, i] += n_eq / (n_cal + 1)

        return p


class ConformalModel:
    '''The 50 aggregated ICPs behind one endpoint.'''

    __slots__ = ("name", "icps", "classes")

    def __init__(self, name, icps, classes):
        self.name = name
        self.icps = icps
        self.classes = classes

    def predict_pvalues(self, x, smoothing=False, rng=None):
        '''p-values per compound and class, median-aggregated over the ICPs.

        Parameters:
        ---------
        x (numpy array): the (n_compounds, 119) descriptor matrix.
        smoothing (bool): add the upstream smoothing term to each p-value.
        rng (numpy RandomState): source of the smoothing draws.

        Returns:
        ----------
        p_values (numpy array): shape (n_compounds, 2), columns [inactive, active].
        '''

        if smoothing and rng is None:
            rng = np.random
        # scikit-learn's tree code casts the input to float32 before walking the
        # nodes, so thresholds are effectively compared against float32-rounded
        # descriptors. Match that or borderline splits take the other branch.
        x = np.ascontiguousarray(x, dtype=np.float32)
        stacked = np.dstack([icp.p_values(x, smoothing, rng) for icp in self.icps])
        return np.median(stacked, axis=2)         # Consensus' agg_func


def load_model(path: str) -> ConformalModel:
    '''Read one endpoint's pickle without importing scikit-learn or pandas.

    Parameters:
    ---------
    path (str): the model pickle.

    Returns:
    ----------
    model (ConformalModel): ready to predict on a descriptor matrix.
    '''

    with open(path, "rb") as fh:
        raw = _StubUnpickler(fh).load()

    icps = []
    for icp in raw.predictor.predictors:
        forest = icp.nc_function.model.model
        if type(forest).__name__ != "RandomForestClassifier":
            raise PredictionError("unexpected underlying model: {}".format(type(forest).__name__))
        if getattr(icp.nc_function, "normalizer", None) is not None:
            raise PredictionError("model uses a normalizer, which this runtime does not implement")
        if not icp.conditional:
            raise PredictionError("model is not label-conditional; check its condition lambda")
        icps.append(_Icp(_Forest([e.tree_ for e in forest.estimators_]),
                         icp.cal_scores, icp.classes))

    model = ConformalModel(raw.model_name, icps, raw.predictor.classes)
    del raw                                        # drop the embedded training frames
    return model


def region(p_inactive, p_active, significance):
    '''One compound's conformal prediction region at the given significance.'''
    if np.isnan(p_inactive) or np.isnan(p_active):
        return ""
    inactive, active = p_inactive > significance, p_active > significance
    if inactive and active:
        return "both"
    if inactive:
        return "inactive"
    if active:
        return "active"
    return "empty"
