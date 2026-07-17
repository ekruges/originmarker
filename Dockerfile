# OriginMarker - served as a subpage at https://ezrakruger.cc/originmarker/
# See deploy/README-deploy.md for the runbook.

# ---- stage 1: build the SPA -------------------------------------------------
FROM node:22-alpine AS web

WORKDIR /web
# Lockfile first: this layer only rebuilds when deps actually change.
COPY web/package.json web/package-lock.json ./
RUN npm ci
COPY web/ ./
# vite base is './' (relative), so the bundle works under /originmarker/ without
# an absolute origin baked in. Output: /web/dist
RUN npm run build


# ---- stage 2: runtime -------------------------------------------------------
FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# primer3-py is GPLv2 and this repo is Apache 2.0, so it is NOT in requirements.txt and the
# default image does not carry it: an image bundling both is a distributed combined work
# under terms this repo cannot meet. Installing it on a server you run is not distribution.
# So it is an operator's switch, exactly like the API keys:
#
#     docker compose build --build-arg WITH_PRIMERS=1
#
# Off, primers.py reports itself unavailable, panels build unchanged and carry no primers.
# If you redistribute an image built with this on, the combined work is GPLv2. Your call.
ARG WITH_PRIMERS=0
COPY requirements-primers.txt ./
RUN if [ "$WITH_PRIMERS" = "1" ]; then pip install --no-cache-dir -r requirements-primers.txt; fi

# Deploy-breaker guard. cloudflared has no path-rewrite, so the tunnel hands us the
# prefix intact ("/originmarker/api/health"). Only Starlette >=0.33 (FastAPI >=0.109)
# strips the app's root_path before routing. On older versions the LAN check
# (/api/health) still returns 200 while every public URL 404s - a silent, miserable
# bug. requirements.txt pins 0.139.0 today; this fails the build if that ever slips.
RUN python -c "import fastapi; v=tuple(int(p) for p in fastapi.__version__.split('.')[:3]); assert v>=(0,109,0), 'fastapi '+fastapi.__version__+' < 0.109.0: no root_path stripping, /originmarker/* would 404 through the tunnel'"

# panelbuilder.py + genetic_map.py sit at the root: genetic_map resolves its maps as
# Path(__file__).parent/"data"/"maps", so data/ must stay beside them. 23MB, bundled
# on purpose - never re-downloaded at runtime.
#
# primers.py ships whatever WITH_PRIMERS says: the module is ours and Apache 2.0, only the
# primer3 dependency is GPLv2. app/jobs.py imports it at module load, so leaving it out
# does not disable primers, it stops the app from starting.
COPY panelbuilder.py genetic_map.py build_info.py primers.py ./
COPY data/ ./data/
COPY app/ ./app/
# app/main.py resolves DIST as <repo>/web/dist -> /app/web/dist
COPY --from=web /web/dist ./web/dist

# Owned by 1000 so a fresh named volume inherits that ownership on first run
# (Docker seeds an empty named volume from the image dir, perms included).
RUN mkdir -p /cache && chown -R 1000:1000 /cache

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    ROOT_PATH=/originmarker \
    PANELBUILDER_CACHE=/cache

EXPOSE 8000

# Some NAS filesystems refuse a bind mount unless the container uid matches the host
# file owner, even at 777. Nothing here binds a host path: the cache is a named
# volume - but uid 1000 keeps that true if a bind is ever added.
USER 1000:1000

# /api/health is cheap: it touches no network. The Ensembl release number is warmed on a
# background thread at startup and read from memory, so health never blocks on upstream.
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/api/health', timeout=4)"

# NOTE: deliberately NO --root-path flag. uvicorn's --root-path *prepends* the prefix to
# an already-unstripped path, yielding /originmarker/originmarker/... -> 404. app/main.py
# reads ROOT_PATH into FastAPI(root_path=...), which strips the prefix only when present,
# so this serves both /api/health (LAN) and /originmarker/api/health (tunnel).
#
# --forwarded-allow-ips is required, not cosmetic. cloudflared reaches us from the docker
# bridge, not loopback, and uvicorn only honours X-Forwarded-Proto from addresses in this
# list (default: 127.0.0.1). Without it Starlette believes the request is http, and the
# trailing-slash redirect on /originmarker hands the user an http:// Location - a silent
# downgrade off TLS on a public URL. Trusting * is safe here because the only routes in
# are the tunnel and the LAN debug port; the worst a spoofed header can do is change a
# redirect's scheme.
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", \
     "--proxy-headers", "--forwarded-allow-ips", "*"]
