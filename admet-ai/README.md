# ADMET-AI

[ADMET-AI](https://github.com/swansonk14/admet_ai) (Swanson et al., 2024)
served over HTTP, with the same routes, inputs and answer shape as the QSAR
models in this repository. From chemical structure it predicts 41 absorption,
distribution, metabolism, excretion and toxicity (ADMET) endpoints and
computes 11 physicochemical properties. This folder holds everything needed
to build, run and deploy it.

| | |
|---|---|
| **Predicts** | 41 ADMET endpoints and 11 physicochemical properties |
| **Model** | `admet_ai` 2.0.1: Chemprop-RDKit graph neural networks, two ensembles of five |
| **URL** | https://admet-ai.serve.scilifelab.se |
| **Image** | `ths-admet-ai` |
| **Port** | 8080 |

## Files

| File | What it is |
|---|---|
| `admet.py` | The runtime: loads `admet_ai.ADMETModel` and predicts a list of SMILES, one row per input. |
| `server.py` | The HTTP API. |
| `validate.py` | Checks the model against the predictions the ADMET-AI authors ship. The Dockerfile runs it. |
| `model.json` | The model's name, what it predicts and its Serve subdomain. |
| `requirements.txt` | The Python packages, pinned. |
| `start-script.sh` | Starts the API on port 8080. The container's command. |
| `Dockerfile` | Builds the image. |

The model weights, the property descriptions and the DrugBank predictions
`validate.py` checks against come inside the `admet-ai` package, so there is
no model file here.

## Use the API

| Route | What it does |
|---|---|
| `GET /` | The model, its version, how to call it, the citation, and under `properties` what every column is: units, task type, training-set size and the authors' test-set metrics. |
| `GET /health` | `200` once the model is loaded. |
| `GET /predict?smiles=<SMILES>` | Predicts one compound. Repeat `smiles=`, or comma-separate, for more. |
| `POST /predict` | Predicts from JSON, a CSV/TSV upload or a CSV/TSV body, as below. |

`/docs` has interactive documentation, where a request can be tried in the
browser.

```bash
URL=https://admet-ai.serve.scilifelab.se

# One SMILES, or several
curl "$URL/predict?smiles=CCO"
curl "$URL/predict" --json '{"smiles": "CCO"}'
curl "$URL/predict" --json '{"smiles": ["CCO", "c1ccccc1O"]}'

# A CSV or TSV file with a column whose name contains "smiles"
curl "$URL/predict" -F file=@compounds.csv
curl "$URL/predict" -H 'Content-Type: text/csv' --data-binary @compounds.csv

# A CSV back instead of JSON
curl "$URL/predict?format=csv" -F file=@compounds.csv -o ADMET-AI_results.csv
```

The answer has one row per SMILES, in input order. For ethanol, phenol and a
SMILES that RDKit cannot parse, with 9 of each row's 52 property columns
shown:

```json
{
  "model": "ADMET-AI",
  "predicts": "Absorption, distribution, metabolism, excretion and toxicity (ADMET) endpoints and physicochemical properties",
  "confidence": null,
  "n_compounds": 3,
  "unparsed": ["not_a_smiles("],
  "predictions": [
    {"smiles": "CCO", "endpoint": "ADMET-AI", "molecular_weight": 46.069, "logP": -0.0014000000000000123, "QED": 0.4068079656553945, "AMES": 0.11980573, "hERG": 0.0101208985, "DILI": 0.1729049, "BBB_Martins": 0.979434, "Solubility_AqSolDB": 0.92181903, "LD50_Zhu": 1.6471183},
    {"smiles": "c1ccccc1O", "endpoint": "ADMET-AI", "molecular_weight": 94.11299999999999, "logP": 1.3922, "QED": 0.514729544768675, "AMES": 0.05648624, "hERG": 0.047842465, "DILI": 0.17205094, "BBB_Martins": 0.8335414, "Solubility_AqSolDB": -0.49433836, "LD50_Zhu": 2.1948981},
    {"smiles": "not_a_smiles(", "endpoint": "ADMET-AI", "molecular_weight": null, "logP": null, "QED": null, "AMES": null, "hERG": null, "DILI": null, "BBB_Martins": null, "Solubility_AqSolDB": null, "LD50_Zhu": null}
  ]
}
```

- The columns are `admet_ai`'s own names, in its order, as
  `ADMETModel(drugbank_path=None).predict()` returns them in Python: the 11
  physicochemical properties, then the 41 ADMET endpoints. `GET /` describes
  every one under `properties`. The `<property>_drugbank_approved_percentile`
  columns that `ADMETModel()` adds by default are left out.
- A **classification** endpoint (`AMES`, `hERG`, `DILI`, the CYP and Tox21
  endpoints, ...) is the predicted probability of the positive label, from 0
  to 1. It is not a label, and no threshold is applied.
- A **regression** endpoint is a value in the units listed under
  `properties`, e.g. `Solubility_AqSolDB` in log(mol/L) and `LD50_Zhu` in
  log(1/(mol/kg)). Values are not clipped to a physical range: `Half_Life_Obach`
  can come back negative. Check each endpoint's test-set metrics in
  `properties` before relying on it. For `Half_Life_Obach` and `VDss_Lombardo`
  the authors report a negative R².
- There is **no applicability domain**. Every structure RDKit parses gets
  values, including salts, mixtures, metals and inorganics, so standardize
  structures first. `confidence` is `null`, because this is not a conformal
  model.
- A SMILES that RDKit cannot parse, or that parses to no atoms, gets a row of
  `null`s and is listed in `unparsed`. The rest of the request still succeeds.
- A request holds at most 2000 compounds and 50 MB; split larger inputs. With
  2 vCPU, 1000 compounds take about 12–27 s. A malformed request gets `400`
  (no SMILES, no `smiles` column, not a CSV or TSV), `413` (too large) or
  `415` (another content type), with the reason in `detail`.
- Values agree with `admet_ai` run in Python. The model computes in float32,
  so the last digit can vary with the CPU and with which other molecules
  share the request.

### Differences from the QSAR models

The routes, the inputs, the envelope (`model`, `predicts`, `confidence`,
`n_compounds`, `unparsed`, `predictions`), the order of the rows, the error
codes and the limits are the same as the QSAR models'. The client in
chemsafe-agent's `qsar_modelling/scripts/utils.py` sends to and checks this
API unchanged: `_request_predictions("ADMET-AI", smiles)` finds it at
`admet-ai` and gets back one row per SMILES. What differs is inside a row:

| | QSAR models | ADMET-AI |
|---|---|---|
| Columns after `smiles`, `endpoint` | `confidence`, `p_inactive`, `p_active`, `prediction` | the 52 property columns above |
| A row for an unparsable SMILES | `p_inactive`, `p_active` and `prediction` are `null` | every property column is `null` |
| `confidence` | the model's conformal confidence | `null` |
| `GET /` | | also `version` and `properties` |

So code that reads `row["prediction"]`, such as `predict_endpoint`, needs
its own reader for ADMET-AI rows.

## Build and run

```bash
docker build --platform linux/amd64 -t ths-admet-ai:1.0.0 .
docker run --rm -p 8080:8080 ths-admet-ai:1.0.0
curl "http://localhost:8080/predict?smiles=CCO"
```

The build runs `validate.py`, as the app's user. It predicts all 2,845
DrugBank approved drugs and compares them with the values the ADMET-AI
authors ship in the package. The RDKit properties must match to 1e-9. The
model outputs must match to within 1e-4 of each property's spread across the
drugs: on arm64 they match bit for bit, and on amd64 to within 3e-6. Two known
exceptions are allowed: today's RDKit computes QED differently for two
deuterated drugs. The check also requires unparsable SMILES placed among the
drugs to come back as empty rows in place. If anything differs, the build
fails.

torch comes from PyTorch's CPU-only index, at the version
`requirements.txt` pins. RDKit is pinned to 2025.9.6, the last 2025.09
release: `admet-ai` 2.0.1 needs at least 2025.9.5, so the QSAR models'
2025.9.3 cannot be used. On the DrugBank set it gives the same
physicochemical values as 2026.3.

The server sizes torch's threads to the container's CPU limit. Set
`ADMET_NUM_THREADS` to override that, and `MAX_COMPOUNDS` or `MAX_BODY_BYTES`
to change the request limits.

Without Docker, with Python 3.12: `pip install -r requirements.txt`, then
`./start-script.sh`. On Linux, first install torch from the CPU-only index as
the Dockerfile does, or pip fetches the CUDA build.

## Citation

> Swanson, K.; Walther, P.; Leitz, J.; Mukherjee, S.; Wu, J. C.;
> Shivnaraine, R. V.; Zou, J. ADMET-AI: a machine learning ADMET platform for
> evaluation of large-scale chemical libraries. *Bioinformatics* **2024**,
> *40* (7), btae416.
> DOI: [10.1093/bioinformatics/btae416](https://doi.org/10.1093/bioinformatics/btae416)
