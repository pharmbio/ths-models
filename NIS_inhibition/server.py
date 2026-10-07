'''HTTP API for one RiskMix conformal QSAR model. One model per container: on
SciLifeLab Serve each has its own URL, https://<subdomain>.serve.scilifelab.se.

The model is model.pkl beside this file, and model.json names it and sets what
it predicts and its conformal confidence. Set MODEL_DIR to read both from
another folder.

    GET  /                     what this model is and how to call it
    GET  /health               200 once the model is loaded
    GET  /predict?smiles=CCO   repeat `smiles`, or comma-separate, for several
    POST /predict              JSON {"smiles": "CCO"} or {"smiles": ["CCO", ...]},
                               a CSV/TSV upload (multipart field "file"), or a
                               CSV/TSV body (Content-Type text/csv or
                               text/tab-separated-values)

A table needs a column whose name contains "smiles". Add ?format=csv to
/predict for a CSV back instead of JSON. Interactive docs are at /docs.
'''

import csv
import io
import json
import logging
import os
from typing import List, Literal

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import JSONResponse, Response
from rdkit import RDLogger
from starlette.concurrency import run_in_threadpool

from conformal import descriptor_matrix, load_model, region

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.environ.get("MODEL_DIR", HERE)
# A 2000-compound request takes ~20 s on the largest model with 2 vCPU, well
# inside a proxy's timeout. Clients split larger inputs into batches.
MAX_COMPOUNDS = int(os.environ.get("MAX_COMPOUNDS", "2000"))
MAX_BODY_BYTES = int(os.environ.get("MAX_BODY_BYTES", str(50 * 1024 * 1024)))
SMOOTHING = os.environ.get("THS_SMOOTHING", "0") == "1"
SEED = int(os.environ.get("THS_SEED", "0"))

COLUMNS = ["smiles", "endpoint", "confidence", "p_inactive", "p_active", "prediction"]
CITATION = (
    "Dracheva, E.; Norinder, U.; Rydén, P.; Engelhardt, J.; Weiss, J. M.; Andersson, P. L. "
    "In Silico Identification of Potential Thyroid Hormone System Disruptors among Chemicals "
    "in Human Serum and Chemicals with a High Exposure Index. Environ. Sci. Technol. 2022. "
    "doi:10.1021/acs.est.1c07762"
)

logger = logging.getLogger("uvicorn.error")
# Unparsable SMILES are reported in the response; RDKit's own lines would only
# fill the log.
RDLogger.DisableLog("rdApp.*")

with open(os.path.join(MODEL_DIR, "model.json")) as fh:
    INFO = json.load(fh)
ENDPOINT = INFO["model"]
CONFIDENCE = float(INFO["confidence"])
# 1 - 0.8 is 0.19999999999999996 in binary floating point, and these p-values
# are rationals with small denominators (n_gt/(n_cal+1)), so an exact tie at
# the threshold is common rather than a corner case. Without the rounding,
# p == 0.2 counts as above 0.2 and the label is kept, widening the region on ties.
SIGNIFICANCE = round(1 - CONFIDENCE, 12)
# Loaded before uvicorn opens the port, so the app is reachable only when ready.
MODEL = load_model(os.path.join(MODEL_DIR, "model.pkl"))
logger.info("Loaded %s (confidence %s)", ENDPOINT, CONFIDENCE)

app = FastAPI(
    title="{} QSAR model".format(ENDPOINT),
    description="{}, as a conformal prediction region at confidence {}. {}".format(
        INFO["predicts"], CONFIDENCE, CITATION),
)


@app.exception_handler(Exception)
async def _unexpected(request: Request, exc: Exception):
    logger.exception("Prediction failed")
    return JSONResponse(status_code=500, content={"detail": "{}: {}".format(type(exc).__name__, exc)})


# Input handling

def _split(text: str) -> List[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def _smiles_from_json(payload) -> List[str]:
    smiles = payload.get("smiles") if isinstance(payload, dict) else None
    if isinstance(smiles, str):
        return _split(smiles)
    if isinstance(smiles, list) and all(isinstance(item, str) for item in smiles):
        return [item.strip() for item in smiles if item.strip()]
    raise HTTPException(400, 'Send JSON as {"smiles": "<SMILES>"} or {"smiles": ["<SMILES>", ...]}.')


def _smiles_from_table(data: bytes, tab: bool, source: str) -> List[str]:
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        raise HTTPException(400, "{} is not UTF-8 text.".format(source))
    rows = csv.reader(io.StringIO(text), delimiter="\t" if tab else ",")
    header = next(rows, None)
    if not header:
        raise HTTPException(400, "{} is empty.".format(source))
    column = next((i for i, name in enumerate(header) if "smiles" in name.lower()), None)
    if column is None:
        raise HTTPException(400, "No column containing 'smiles' found in {}. Columns present: {}.".format(
            source, ", ".join(header)))
    return [row[column].strip() for row in rows if len(row) > column and row[column].strip()]


# Prediction

def _predict(smiles_list: List[str]):
    '''One row per input SMILES, in input order, and the SMILES RDKit rejected.

    A SMILES RDKit cannot parse is a row with no prediction, listed in
    `unparsed`, not a failed request: a client sending a large input in
    batches gets every batch back, even one where nothing parses.
    '''

    if not smiles_list:
        raise HTTPException(400, "No valid SMILES strings were provided.")
    if len(smiles_list) > MAX_COMPOUNDS:
        raise HTTPException(413, "{} compounds in one request; the limit is {}. Split the input "
                                 "into batches.".format(len(smiles_list), MAX_COMPOUNDS))

    descriptors, bad_index, bad_smiles = descriptor_matrix(smiles_list)

    # Put predictions back on the rows they came from, so an unparsable
    # structure leaves a gap instead of shifting the rest up by one.
    full = np.full((len(smiles_list), 2), np.nan)
    if descriptors.shape[0]:
        rng = np.random.RandomState(SEED) if SMOOTHING else None
        bad = set(bad_index)
        full[[i for i in range(len(smiles_list)) if i not in bad]] = MODEL.predict_pvalues(
            descriptors, smoothing=SMOOTHING, rng=rng)

    rows = []
    for smiles, (p_inactive, p_active) in zip(smiles_list, full):
        parsed = not np.isnan(p_inactive)
        rows.append({
            "smiles": smiles,
            "endpoint": ENDPOINT,
            "confidence": CONFIDENCE,
            "p_inactive": float(p_inactive) if parsed else None,
            "p_active": float(p_active) if parsed else None,
            "prediction": region(p_inactive, p_active, SIGNIFICANCE) or None,
        })
    logger.info("Predicted %d compound(s), %d unparsable", len(rows), len(bad_smiles))
    return rows, bad_smiles


def _respond(rows, unparsed, fmt):
    if fmt == "csv":
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        return Response(buffer.getvalue(), media_type="text/csv", headers={
            "Content-Disposition": 'attachment; filename="{}_results.csv"'.format(ENDPOINT)})
    return {
        "model": ENDPOINT,
        "predicts": INFO["predicts"],
        "confidence": CONFIDENCE,
        "n_compounds": len(rows),
        "unparsed": unparsed,
        "predictions": rows,
    }


# Routes

@app.get("/")
def info():
    return {
        "model": ENDPOINT,
        "predicts": INFO["predicts"],
        "confidence": CONFIDENCE,
        "output": "A conformal prediction region per compound: active, inactive, both (undecided "
                  "at this confidence) or empty (outside the applicability domain), with the "
                  "p-values p_inactive and p_active. Not a probability.",
        "max_compounds_per_request": MAX_COMPOUNDS,
        "usage": {
            "GET /predict?smiles=<SMILES>": "one compound; repeat smiles= or comma-separate for more",
            "POST /predict (application/json)": '{"smiles": "<SMILES>"} or {"smiles": ["<SMILES>", ...]}',
            "POST /predict (multipart/form-data)": 'a CSV or TSV file in the field "file", with a column whose name contains "smiles"',
            "POST /predict (text/csv or text/tab-separated-values)": "the same table as the request body",
            "?format=csv": "return a CSV instead of JSON",
            "/docs": "interactive API documentation",
        },
        "citation": CITATION,
    }


@app.get("/health")
def health():
    return {"status": "ok", "model": ENDPOINT}


@app.get("/predict")
def predict_get(smiles: List[str] = Query(default=[]), format: Literal["json", "csv"] = "json"):
    rows, unparsed = _predict([item for value in smiles for item in _split(value)])
    return _respond(rows, unparsed, format)


@app.post("/predict", openapi_extra={"requestBody": {"required": True, "content": {
    "application/json": {"schema": {"type": "object", "required": ["smiles"], "properties": {
        "smiles": {"oneOf": [{"type": "string"}, {"type": "array", "items": {"type": "string"}}]}}}},
    "multipart/form-data": {"schema": {"type": "object", "required": ["file"], "properties": {
        "file": {"type": "string", "format": "binary"}}}},
    "text/csv": {"schema": {"type": "string"}},
    "text/tab-separated-values": {"schema": {"type": "string"}},
}}})
async def predict_post(request: Request, format: Literal["json", "csv"] = "json"):
    if int(request.headers.get("content-length") or 0) > MAX_BODY_BYTES:
        raise HTTPException(413, "The request body is larger than {} bytes.".format(MAX_BODY_BYTES))

    content_type = request.headers.get("content-type", "").split(";")[0].strip().lower()
    if content_type == "application/json":
        try:
            payload = await request.json()
        except ValueError:
            raise HTTPException(400, "The request body is not valid JSON.")
        smiles = _smiles_from_json(payload)
    elif content_type == "multipart/form-data":
        form = await request.form()
        upload = form.get("file")
        if upload is None or isinstance(upload, str):
            raise HTTPException(400, 'Upload the table in the form field "file".')
        extension = os.path.splitext(upload.filename or "")[1].lower()
        if extension not in (".csv", ".tsv"):
            raise HTTPException(400, "Only CSV or TSV files are supported for SMILES input.")
        smiles = _smiles_from_table(await upload.read(), extension == ".tsv", upload.filename)
    elif content_type in ("text/csv", "text/tab-separated-values"):
        smiles = _smiles_from_table(await request.body(), content_type != "text/csv", "The request body")
    else:
        raise HTTPException(415, "Send application/json, multipart/form-data (a file in the field "
                                 "\"file\"), text/csv or text/tab-separated-values.")

    rows, unparsed = await run_in_threadpool(_predict, smiles)
    return _respond(rows, unparsed, format)
