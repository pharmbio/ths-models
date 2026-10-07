'''Check this model against p-values from the original 2021 stack.

reference.npy holds the p-values this model's pickle produced on the 100
compounds in descriptors.npy, under its real dependencies (scikit-learn
0.21.3, pandas 1.1.5, nonconformist, cloudpickle 1.6, executing the real
pickled lambdas), with smoothing disabled so the numbers are deterministic.
This recomputes them with conformal.py and requires an exact match. The
Dockerfile runs it, so an image that does not match is never built.

    python validate.py
'''

import os
import sys

import numpy as np

from conformal import load_model

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.environ.get("MODEL_DIR", HERE)


def main():
    x = np.load(os.path.join(HERE, "descriptors.npy"))
    expected = np.load(os.path.join(MODEL_DIR, "reference.npy"))
    produced = load_model(os.path.join(MODEL_DIR, "model.pkl")).predict_pvalues(x, smoothing=False)

    if not np.array_equal(produced, expected):
        sys.exit("MISMATCH with the original stack: max|diff| {:.3e}".format(
            np.abs(produced - expected).max()))
    print("{} compounds: p-values match the original stack exactly".format(x.shape[0]))


if __name__ == "__main__":
    main()
