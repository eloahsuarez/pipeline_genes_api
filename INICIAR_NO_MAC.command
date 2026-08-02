#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

show_error() {
  printf '\nERRO: %s\n\n' "$1"
  printf 'Pressione Enter para fechar...'
  read -r _
  exit 1
}

if ! command -v python3 >/dev/null 2>&1; then
  show_error "Python 3 não foi encontrado. Instale Python 3.11 ou superior pelo instalador oficial do python.org e execute este arquivo novamente."
fi

if ! python3 - <<'PY'
import sys
raise SystemExit(0 if sys.version_info >= (3, 11) else 1)
PY
then
  show_error "É necessário Python 3.11 ou superior."
fi

if ! python3 - <<'PY'
try:
    import tkinter
except Exception:
    raise SystemExit(1)
PY
then
  show_error "O módulo Tkinter não está disponível. No macOS, use preferencialmente o instalador oficial do Python em python.org, que inclui o suporte gráfico necessário."
fi

if [ ! -x ".venv/bin/python" ]; then
  echo "Preparando o ambiente Python..."
  python3 -m venv .venv
  .venv/bin/python -m pip install --upgrade pip
  .venv/bin/python -m pip install -r requirements.txt
elif ! .venv/bin/python -c 'import keyring' >/dev/null 2>&1; then
  echo "Instalando os componentes da atualização..."
  .venv/bin/python -m pip install -r requirements.txt
fi

echo "Abrindo o Gene Conservado..."
exec .venv/bin/python app.py
