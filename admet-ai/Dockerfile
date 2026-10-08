# ADMET-AI: absorption, distribution, metabolism, excretion and toxicity.
# The admet_ai package's Chemprop models behind an HTTP API on :8080
# (server.py), with the same routes and answer as the QSAR models; see README.md.
# Served at https://admet-ai.serve.scilifelab.se.
# Build from this folder (Serve runs linux/amd64):
#   docker build --platform linux/amd64 -t ths-admet-ai:1.0.0 .
FROM python:3.12.14-slim-trixie

ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1

ENV USER=admet
ENV HOME=/home/$USER
# Chemprop imports matplotlib, which needs a cache it can write. Built under
# emulation, root's steps leave $HOME/.cache root-owned; this one is the app's.
ENV MPLCONFIGDIR=$HOME/.matplotlib

RUN useradd -m -u 1000 $USER

WORKDIR $HOME/app
COPY requirements.txt ./
# torch first, CPU-only, at the version requirements.txt pins; pip then keeps
# it rather than fetching PyPI's CUDA build.
RUN pip install --index-url https://download.pytorch.org/whl/cpu "$(grep '^torch==' requirements.txt)" \
    && pip install -r requirements.txt

COPY admet.py server.py validate.py model.json ./
# Executable whatever its mode in the checkout.
COPY --chmod=0755 start-script.sh ./

USER $USER

# Refuse to build an image whose predictions differ from those the ADMET-AI
# authors ship. Run as the app's user, so nothing it needs is root-only.
RUN python validate.py

EXPOSE 8080
HEALTHCHECK --interval=10s --timeout=3s --start-period=60s --retries=3 \
    CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8080/health', timeout=2)"]
# Serve starts the container with the start-up script in WORKDIR.
CMD ["./start-script.sh"]
