#!/bin/bash
# Starts the PXR_agonist API on :8080 (server.py): the Dockerfile's CMD, and
# what Serve runs. The model loads before the port opens, so the app answers
# only once it can predict.
set -euo pipefail
cd "$(dirname "$0")"

echo "Starting PXR_agonist: Pregnane X receptor (PXR) agonism, at confidence 0.85"
exec uvicorn server:app --host 0.0.0.0 --port 8080
