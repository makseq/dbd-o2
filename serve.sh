#!/usr/bin/env bash
# serve.sh — локальный статический сервер для dbd-o2 (2D-плеер + 1D-UI).
#
# ES-модули не грузятся с file:// (CORS), поэтому страницу надо открывать
# только через http. Зависимостей нет — нужен либо python3, либо node.
#
#   ./serve.sh            → http://127.0.0.1:8777/player/
#   ./serve.sh 9000       → другой порт
#   PORT=9000 ./serve.sh  → то же самое
#
# ВАЖНО: сервер поднимается из КОРНЯ репозитория, иначе плеер не найдёт ../data/.

set -euo pipefail
cd "$(dirname "$0")"

PORT="${1:-${PORT:-8777}}"
BASE="http://127.0.0.1:${PORT}"

echo "2D ДБР / O2 — плеер        →  ${BASE}/player/"
echo "  графики отдельно         →  ${BASE}/player/plots-demo.html"
echo "  1D-модель (живой солвер) →  ${BASE}/index.html"
echo "Ctrl+C — остановить"

if command -v python3 >/dev/null 2>&1; then
  exec python3 -m http.server "$PORT" --bind 127.0.0.1
elif command -v node >/dev/null 2>&1; then
  exec node --eval '
    const http = require("http"), fs = require("fs"), path = require("path");
    const TYPES = { ".html":"text/html; charset=utf-8", ".js":"text/javascript",
                    ".mjs":"text/javascript", ".css":"text/css", ".json":"application/json",
                    ".svg":"image/svg+xml", ".png":"image/png", ".md":"text/markdown; charset=utf-8" };
    const root = process.cwd(), port = Number(process.argv[1]);
    http.createServer((req, res) => {
      let rel = decodeURIComponent(req.url.split("?")[0]).replace(/^\/+/, "") || "index.html";
      let file = path.resolve(root, rel);
      if (!file.startsWith(root)) { res.writeHead(403).end(); return; }
      try { if (fs.statSync(file).isDirectory()) file = path.join(file, "index.html"); } catch (e) {}
      fs.readFile(file, (err, buf) => {
        if (err) { res.writeHead(404).end("not found"); return; }
        res.writeHead(200, { "content-type": TYPES[path.extname(file)] || "application/octet-stream" });
        res.end(buf);
      });
    }).listen(port, "127.0.0.1");
  ' "$PORT"
else
  echo "нужен python3 или node" >&2
  exit 1
fi
