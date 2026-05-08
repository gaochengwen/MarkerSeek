#!/usr/bin/env bash
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HOST="${MARKERSEEK_WEB_HOST:-0.0.0.0}"
PORT="${MARKERSEEK_WEB_PORT:-80}"
DATA_DIR="${MARKERSEEK_WEB_DATA:-$APP_DIR/markerseek_jobs}"
MAFFT_BIN="${MARKERSEEK_MAFFT_BIN:-mafft}"

mkdir -p "$DATA_DIR"
export MARKERSEEK_WEB_DATA="$DATA_DIR"
export MARKERSEEK_MAFFT_BIN="$MAFFT_BIN"

exec "$APP_DIR/.venv/bin/markerseek-web" --host "$HOST" --port "$PORT"
