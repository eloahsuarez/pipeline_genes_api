#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

if [ ! -x ".venv/bin/python" ]; then
  echo "O ambiente ainda não foi instalado."
  echo "Execute instalar_mac.command ou INICIAR_NO_MAC.command primeiro."
  printf '\nPressione Enter para fechar...'
  read -r _
  exit 1
fi

exec .venv/bin/python app.py
