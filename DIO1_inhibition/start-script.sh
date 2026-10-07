#!/bin/bash
# Starts the DIO1_inhibition API on :8080 (server.py): the Dockerfile's CMD, and
# what Serve runs. The model loads before the port opens, so the app answers
# only once it can predict.
set -euo pipefail
cd "$(dirname "$0")"

echo "Starting DIO1_inhibition: Type 1 iodothyronine deiodinase (DIO1) inhibition, at confidence 0.8"
exec uvicorn server:app --host 0.0.0.0 --port 8080
