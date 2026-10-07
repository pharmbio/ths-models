# DIO3_inhibition: type 3 iodothyronine deiodinase (DIO3) inhibition, at confidence 0.8.
# A conformal QSAR model behind an HTTP API on :8080 (server.py); see README.md.
# Served at https://dio3-inhibition.serve.scilifelab.se.
# Build from this folder (Serve runs linux/amd64):
#   docker build --platform linux/amd64 -t ths-dio3-inhibition:1.0.0 .
FROM python:3.12.14-slim-trixie

ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1

ENV USER=qsar
ENV HOME=/home/$USER

RUN useradd -m -u 1000 $USER

WORKDIR $HOME/app
COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY conformal.py server.py validate.py descriptors.npy model.json model.pkl reference.npy ./
# Refuse to build an image whose predictions differ from the original stack.
RUN python validate.py

# Executable whatever its mode in the checkout.
COPY --chmod=0755 start-script.sh ./

USER $USER

EXPOSE 8080
HEALTHCHECK --interval=10s --timeout=3s --start-period=60s --retries=3 \
    CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8080/health', timeout=2)"]
# Serve starts the container with the start-up script in WORKDIR.
CMD ["./start-script.sh"]
