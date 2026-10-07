# CAR_antagonist

A conformal QSAR model of constitutive androstane receptor (CAR) antagonism,
served over HTTP. It predicts from chemical structure, and is one of the
RiskMix thyroid hormone system models of Dracheva et al. (2022). This folder
holds everything needed to build, run and deploy it.

| | |
|---|---|
| **Predicts** | Constitutive androstane receptor (CAR) antagonism |
| **Confidence** | 0.8 |
| **URL** | https://car-antagonist.serve.scilifelab.se |
| **Image** | `ths-car-antagonist` |
| **Port** | 8080 |

## Files

| File | What it is |
|---|---|
| `model.pkl` | The model (256 MB): 50 inductive conformal predictors, each a 300-tree random forest on 119 RDKit descriptors. The original 2021 pickle, unchanged. |
| `model.json` | The model's name, what it predicts, its confidence and its Serve subdomain. |
| `server.py` | The HTTP API. |
| `conformal.py` | The runtime: descriptors, reading the pickle, the conformal prediction. |
| `validate.py` | Checks the model against the original 2021 stack. The Dockerfile runs it. |
| `descriptors.npy`, `reference.npy` | The 100 compounds `validate.py` predicts, and the p-values the original stack gave them. |
| `requirements.txt` | The Python packages, pinned. |
| `start-script.sh` | Starts the API on port 8080. The container's command. |
| `Dockerfile` | Builds the image. |

## Use the API

| Route | What it does |
|---|---|
| `GET /` | The model, what it predicts, its confidence, how to call it and the citation. |
| `GET /health` | `200` once the model is loaded. |
| `GET /predict?smiles=<SMILES>` | Predicts one compound. Repeat `smiles=`, or comma-separate, for more. |
| `POST /predict` | Predicts from JSON, a CSV/TSV upload or a CSV/TSV body, as below. |

`/docs` has interactive documentation, where a request can be tried in the
browser.

```bash
URL=https://car-antagonist.serve.scilifelab.se

# One SMILES, or several
curl "$URL/predict?smiles=CCO"
curl "$URL/predict" --json '{"smiles": "CCO"}'
curl "$URL/predict" --json '{"smiles": ["CCO", "c1ccccc1O"]}'

# A CSV or TSV file with a column whose name contains "smiles"
curl "$URL/predict" -F file=@compounds.csv
curl "$URL/predict" -H 'Content-Type: text/csv' --data-binary @compounds.csv

# A CSV back instead of JSON
curl "$URL/predict?format=csv" -F file=@compounds.csv -o CAR_antagonist_results.csv
```

The answer has one row per SMILES, in input order. For ethanol, phenol and a
SMILES that RDKit cannot parse:

```json
{
  "model": "CAR_antagonist",
  "predicts": "Constitutive androstane receptor (CAR) antagonism",
  "confidence": 0.8,
  "n_compounds": 3,
  "unparsed": ["not_a_smiles("],
  "predictions": [
    {"smiles": "CCO", "endpoint": "CAR_antagonist", "confidence": 0.8, "p_inactive": 0.32680320569902044, "p_active": 0.25, "prediction": "both"},
    {"smiles": "c1ccccc1O", "endpoint": "CAR_antagonist", "confidence": 0.8, "p_inactive": 0.201246660730187, "p_active": 0.38888888888888884, "prediction": "both"},
    {"smiles": "not_a_smiles(", "endpoint": "CAR_antagonist", "confidence": 0.8, "p_inactive": null, "p_active": null, "prediction": null}
  ]
}
```

- `prediction` is a conformal prediction region at confidence 0.8, not a
  probability:
  - `active` or `inactive`: the one label the model cannot rule out.
  - `both`: undecided at this confidence. A data gap, not a borderline result.
  - `empty`: outside the applicability domain. The model says nothing.
- `p_inactive` and `p_active` are conformal p-values and do not sum to 1. A
  label stays in the region when its p-value is above 0.2
  (1 − confidence).
- A SMILES that RDKit cannot parse gets a row of `null`s and is listed in
  `unparsed`. The rest of the request still succeeds.
- A request holds at most 2000 compounds and 50 MB; split larger inputs. A
  malformed request gets `400` (no SMILES, no `smiles` column, not a CSV or
  TSV), `413` (too large) or `415` (another content type), with the reason in
  `detail`.

## Build and run

```bash
docker build --platform linux/amd64 -t ths-car-antagonist:1.0.0 .
docker run --rm -p 8080:8080 ths-car-antagonist:1.0.0
curl "http://localhost:8080/predict?smiles=CCO"
```

The build runs `validate.py`. It predicts the 100 reference compounds and
requires p-values identical to those the original 2021 stack produced
(scikit-learn 0.21.3, nonconformist, the real pickled lambdas). If they differ,
the build fails.

Without Docker, with Python 3.12: `pip install -r requirements.txt`, then
`./start-script.sh`.

## Deploy on SciLifeLab Serve

1. Build for `linux/amd64`, which Serve runs, and push the image to a public
   registry:

   ```bash
   docker build --platform linux/amd64 -t ghcr.io/<user>/ths-car-antagonist:1.0.0 .
   docker push ghcr.io/<user>/ths-car-antagonist:1.0.0
   ```

2. In a Serve project, create a custom app with:
   - **Image**: `ghcr.io/<user>/ths-car-antagonist:1.0.0`
   - **Port**: `8080`
   - **Subdomain**: `car-antagonist`

   It needs no persistent storage and no settings. Serve's default 2 vCPU and
   4 GB are enough: the app peaks at about 0.5 GB while the model loads.
3. For an update, build and push a new tag (`1.0.1`, ...) and change the app's
   image. Serve does not redeploy a tag it has already run.

The image follows Serve's rules: user `qsar` with uid 1000, `start-script.sh`
in `WORKDIR` as its command, and one port (8080) in the range 3000–9999.

ChemSafeAgent's `qsar_modelling` skill finds this app at `car-antagonist`: the
model's name in lower case, with `-` for `_`. To use another pattern for all
models, set `QSAR_API_URL` in the agent's environment, e.g.
`https://ths-{subdomain}.serve.scilifelab.se`.

## Pinned versions

`requirements.txt` pins RDKit 2025.9.3 and NumPy 2.4.0, the versions in
chemsafe-agent's `requirements.txt`. Compare predictions before and after
changing either:

- **RDKit** computes the descriptors, and `validate.py` does not check them:
  it starts from fixed descriptors. RDKit 2026.03 changed the H-bond acceptor
  count for N-heterocycles
  ([#8997](https://github.com/rdkit/rdkit/releases/tag/Release_2026_03_1)),
  one of the 119 descriptors, which changes p-values for compounds such as
  1-vinylimidazole. The model was trained with the earlier definition.
- **NumPy** 2.4 already warns (`VisibleDeprecationWarning`) when it reads the
  2021 pickle; a later version may refuse it.

## Citation

> Dracheva, E.; Norinder, U.; Rydén, P.; Engelhardt, J.; Weiss, J. M.;
> Andersson, P. L. *In Silico* Identification of Potential Thyroid Hormone
> System Disruptors among Chemicals in Human Serum and Chemicals with a High
> Exposure Index. *Environ. Sci. Technol.* **2022**.
> DOI: [10.1021/acs.est.1c07762](https://doi.org/10.1021/acs.est.1c07762)
