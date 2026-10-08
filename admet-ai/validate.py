'''Check this model against the predictions the ADMET-AI authors ship.

The admet_ai package carries drugbank_approved.csv: its own predictions for
the 2,845 DrugBank approved drugs. This predicts those drugs through
admet.py, the same path server.py
serves, and requires every property to match: the RDKit properties exactly,
the model outputs within float32 noise. An unparsable SMILES is slipped in
between the drugs, so the rows must also come back in input order with a gap
where it was. The Dockerfile runs it, so an image that does not match is
never built.

    python validate.py
'''

import sys

import numpy as np
import pandas as pd

from admet import Predictor
from admet_ai.constants import DEFAULT_DRUGBANK_PATH
from admet_ai.physchem import PHYSCHEM_PROPERTY_TO_FUNCTION

# The model outputs are float32 networks, and how the arithmetic is ordered
# depends on the CPU: on arm64 they match the reference bit for bit, on amd64
# to within 3e-6 of each property's spread across the drugs. Allowed: 1e-4 of
# that spread. A changed model or featurizer moves them by 1e-2 or more.
MODEL_TOLERANCE = 1e-4
# The RDKit properties are float64 and differ only in the last digits (1e-13).
RDKIT_RTOL = 1e-9
# Today's RDKit computes QED differently from the one that made the reference
# for these two deuterated drugs (0.496 -> 0.521 and 0.805 -> 0.853). Every
# other property of theirs still has to match.
KNOWN_DIFFERENCES = {("DB16650", "QED"), ("DB12161", "QED")}
UNPARSABLE = "not_a_smiles("


def main():
    reference = pd.read_csv(DEFAULT_DRUGBANK_PATH)
    predictor = Predictor()
    properties = [name for name in predictor.columns if name in reference.columns]

    smiles = []
    for i, drug in enumerate(reference["smiles"]):
        smiles.append(drug)
        if i % 500 == 0:
            smiles.append(UNPARSABLE)
    rows, unparsed = predictor.predict(smiles)

    gaps = [i for i, s in enumerate(smiles) if s == UNPARSABLE]
    if unparsed != [UNPARSABLE] * len(gaps) or any(any(v is not None for v in rows[i]) for i in gaps):
        sys.exit("MISMATCH: the unparsable SMILES did not come back as empty rows in place")
    produced = pd.DataFrame([rows[i] for i, s in enumerate(smiles) if s != UNPARSABLE],
                            columns=predictor.columns)

    failures, worst = [], 0.0
    for name in properties:
        got = produced[name].to_numpy(dtype=float)
        expected = reference[name].to_numpy(dtype=float)
        known = reference["id"].isin([drug for drug, prop in KNOWN_DIFFERENCES if prop == name]).to_numpy()
        if name in PHYSCHEM_PROPERTY_TO_FUNCTION:
            bad = ~np.isclose(got, expected, rtol=RDKIT_RTOL, atol=0)
        else:
            spread = np.std(expected)
            bad = ~(np.abs(got - expected) <= MODEL_TOLERANCE * spread)
            worst = max(worst, np.max(np.abs(got - expected)[~known]) / spread)
        for i in np.flatnonzero(bad & ~known):
            failures.append("{} {} ({}): expected {}, got {}".format(
                reference.at[i, "id"], name, reference.at[i, "name"], expected[i], got[i]))

    if failures:
        sys.exit("MISMATCH with the reference, {} value(s):\n  {}".format(
            len(failures), "\n  ".join(failures[:20])))
    print("{} drugs x {} properties match the reference (model outputs within {:.1e} of their spread); "
          "{} unparsable SMILES came back as empty rows in place".format(
              len(reference), len(properties), worst, len(gaps)))


if __name__ == "__main__":
    main()
