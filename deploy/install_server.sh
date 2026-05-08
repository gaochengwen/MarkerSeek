#!/usr/bin/env bash
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HOST="${MARKERSEEK_WEB_HOST:-0.0.0.0}"
PORT="${MARKERSEEK_WEB_PORT:-80}"
DATA_DIR="${MARKERSEEK_WEB_DATA:-$APP_DIR/markerseek_jobs}"
MAFFT_BIN="${MARKERSEEK_MAFFT_BIN:-mafft}"
INSTALL_SERVICE="${MARKERSEEK_INSTALL_SERVICE:-1}"

need_cmd() {
  command -v "$1" >/dev/null 2>&1
}

echo "MarkerSeek deploy directory: $APP_DIR"

if ! need_cmd python3; then
  echo "python3 is required but was not found." >&2
  exit 1
fi

python3 - <<'PY'
import sys
if not ((3, 11) <= sys.version_info[:2] < (3, 14)):
    raise SystemExit("MarkerSeek requires Python >=3.11 and <3.14.")
PY

if need_cmd apt-get; then
  echo "Installing system dependencies with apt-get..."
  sudo apt-get update
  sudo apt-get install -y python3-venv python3-pip mafft
else
  echo "apt-get not found; skipping system package installation."
fi

if ! need_cmd "$MAFFT_BIN"; then
  echo "MAFFT executable '$MAFFT_BIN' was not found. Install MAFFT or set MARKERSEEK_MAFFT_BIN." >&2
  exit 1
fi

echo "Creating Python virtual environment..."
python3 -m venv "$APP_DIR/.venv"
"$APP_DIR/.venv/bin/python" -m pip install -U pip setuptools wheel
"$APP_DIR/.venv/bin/python" -m pip install "$APP_DIR"

mkdir -p "$DATA_DIR"

if [[ "$INSTALL_SERVICE" == "1" ]] && need_cmd systemctl; then
  echo "Installing systemd service markerseek-web..."
  sudo tee /etc/systemd/system/markerseek-web.service >/dev/null <<EOF
[Unit]
Description=MarkerSeek web server
After=network.target

[Service]
Type=simple
WorkingDirectory=$APP_DIR
Environment=MARKERSEEK_WEB_DATA=$DATA_DIR
Environment=MARKERSEEK_MAFFT_BIN=$MAFFT_BIN
ExecStart=$APP_DIR/.venv/bin/markerseek-web --host $HOST --port $PORT
Restart=always
RestartSec=5

[Install]
WantedBy=multi-user.target
EOF
  sudo systemctl daemon-reload
  sudo systemctl enable markerseek-web
  sudo systemctl restart markerseek-web
  echo "MarkerSeek web service started."
  echo "Open: http://SERVER_IP:${PORT}/markerseek"
else
  echo "Skipping systemd service installation."
  echo "Start manually with: bash $APP_DIR/deploy/start_server.sh"
fi
