#!/usr/bin/env bash
set -euo pipefail

exec /usr/bin/mafft --thread "${MAFFT_THREADS:-4}" "$@"
