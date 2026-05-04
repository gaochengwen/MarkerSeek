# MarkerSeek Web Deployment

This archive contains the MarkerSeek source code and a one-command deployment
script for a new Ubuntu server.

## Quick Deploy

```bash
tar -xzf MarkerSeek-web-deploy.tar.gz
cd MarkerSeek
bash deploy/install_server.sh
```

After installation, open:

```text
http://SERVER_IP/markerseek
```

or, after DNS is configured:

```text
http://www.bioseqhub.cn/markerseek
```

## What The Script Does

- Installs system packages when available through `apt-get`:
  `python3-venv`, `python3-pip`, and `mafft`.
- Creates `.venv` inside the project directory.
- Installs MarkerSeek and web dependencies from `pyproject.toml`.
- Creates the web data directory at `./markerseek_jobs`.
- Installs and starts a `systemd` service named `markerseek-web` on port `80`.

## Useful Commands

```bash
sudo systemctl status markerseek-web
sudo systemctl restart markerseek-web
sudo journalctl -u markerseek-web -f
```

## Environment Variables

Set these before running `deploy/install_server.sh` if needed:

- `MARKERSEEK_WEB_PORT`: default `80`
- `MARKERSEEK_WEB_HOST`: default `0.0.0.0`
- `MARKERSEEK_WEB_DATA`: default `<project>/markerseek_jobs`
- `MARKERSEEK_MAFFT_BIN`: default `mafft`
- `MARKERSEEK_INSTALL_SERVICE`: default `1`; set `0` to skip systemd service installation

If service installation is skipped, run manually:

```bash
bash deploy/start_server.sh
```
