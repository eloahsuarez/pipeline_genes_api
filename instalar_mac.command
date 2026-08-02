#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

pause_on_exit() {
  printf '\nPressione Enter para fechar...'
  read -r _
}
trap pause_on_exit EXIT

if ! command -v python3 >/dev/null 2>&1; then
  echo "ERRO: Python 3 não foi encontrado."
  echo "Instale Python 3.11 ou superior pelo instalador oficial do python.org."
  exit 1
fi

python3 - <<'PY'
import sys
if sys.version_info < (3, 11):
    raise SystemExit("ERRO: é necessário Python 3.11 ou superior.")
try:
    import tkinter
except Exception as exc:
    raise SystemExit(
        "ERRO: Tkinter não está disponível. Use preferencialmente o instalador oficial "
        "do Python em python.org. Detalhe: " + str(exc)
    )
print("Python e Tkinter encontrados.")
PY

rm -rf .venv
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt

echo
 echo "Instalação concluída. Agora abra executar_mac.command."
