'''HTTP API for ADMET-AI, with the same routes, inputs and answer as the RiskMix
QSAR models: on SciLifeLab Serve it has its own URL,
https://admet-ai.serve.scilifelab.se.

The model is the admet_ai package's own (admet.py wraps it), and model.json
names it and says what it predicts. Set MODEL_DIR to read model.json from
another folder.

    GET  /                     what this model is, what each column means, how to call it
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

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import JSONResponse, Response
from starlette.concurrency import run_in_threadpool

from admet import Predictor, property_catalogue

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.environ.get("MODEL_DIR", HERE)
# The same limit as the QSAR models. With 2 vCPU, 1000 compounds (the agent's
# batch) take ~12 s on arm64 and ~27 s on amd64 under emulation, 2000 twice
# that. Clients split larger inputs into batches.
MAX_COMPOUNDS = int(os.environ.get("MAX_COMPOUNDS", "2000"))
MAX_BODY_BYTES = int(os.environ.get("MAX_BODY_BYTES", str(50 * 1024 * 1024)))
# 0 sizes torch's threads to the container's CPU limit.
NUM_THREADS = int(os.environ.get("ADMET_NUM_THREADS", "0"))

CITATION = (
    "Swanson, K.; Walther, P.; Leitz, J.; Mukherjee, S.; Wu, J. C.; Shivnaraine, R. V.; Zou, J. "
    "ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical "
    "libraries. Bioinformatics 2024, 40 (7), btae416. doi:10.1093/bioinformatics/btae416"
)

logger = logging.getLogger("uvicorn.error")

with open(os.path.join(MODEL_DIR, "model.json")) as fh:
    INFO = json.load(fh)
ENDPOINT = INFO["model"]
# Loaded before uvicorn opens the port, so the app is reachable only when ready.
PREDICTOR = Predictor(num_threads=NUM_THREADS)
PROPERTIES = property_catalogue()
COLUMNS = ["smiles", "endpoint"] + PREDICTOR.columns
logger.info("Loaded %s %s: %d columns, %d torch thread(s)",
            ENDPOINT, PREDICTOR.version, len(PREDICTOR.columns), PREDICTOR.num_threads)

N_ADMET = sum(p["category"] != "Physicochemical" for p in PROPERTIES)
OUTPUT = (
    "One row per compound: {} ADMET endpoints predicted by Chemprop-RDKit graph neural networks, "
    "{} physicochemical properties computed by RDKit. A classification endpoint is the predicted "
    "probability of the positive label (0-1), not a label; a regression endpoint is a value in the "
    "units under `properties`. There is no applicability domain: every structure RDKit parses gets "
    "a value, salts, mixtures and metals included. `confidence` is null: this is not a conformal model."
).format(N_ADMET, len(PROPERTIES) - N_ADMET)

app = FastAPI(
    title="{} ADMET model".format(ENDPOINT),
    description="{}. {}".format(INFO["predicts"], CITATION),
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

    A SMILES RDKit cannot parse is a row of None, listed in `unparsed`, not a
    failed request: a client sending a large input in batches gets every batch
    back, even one where nothing parses.
    '''

    if not smiles_list:
        raise HTTPException(400, "No valid SMILES strings were provided.")
    if len(smiles_list) > MAX_COMPOUNDS:
        raise HTTPException(413, "{} compounds in one request; the limit is {}. Split the input "
                                 "into batches.".format(len(smiles_list), MAX_COMPOUNDS))

    values, unparsed = PREDICTOR.predict(smiles_list)
    rows = [dict(zip(COLUMNS, [smiles, ENDPOINT] + row)) for smiles, row in zip(smiles_list, values)]
    logger.info("Predicted %d compound(s), %d unparsable", len(rows), len(unparsed))
    return rows, unparsed


def _respond(rows, unparsed, fmt):
    if fmt == "csv":
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        return Response(buffer.getvalue(), media_type="text/csv", headers={
            "Content-Disposition": 'attachment; filename="{}_results.csv"'.format(ENDPOINT)})
    # Plain numbers and strings already, so FastAPI's encoder has nothing to
    # convert; skipping it saves seconds on a 2000-row answer.
    return JSONResponse({
        "model": ENDPOINT,
        "predicts": INFO["predicts"],
        "confidence": None,
        "n_compounds": len(rows),
        "unparsed": unparsed,
        "predictions": rows,
    })


# Routes

@app.get("/")
def info():
    return {
        "model": ENDPOINT,
        "predicts": INFO["predicts"],
        "version": PREDICTOR.version,
        "confidence": None,
        "output": OUTPUT,
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
        "properties": PROPERTIES,
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
