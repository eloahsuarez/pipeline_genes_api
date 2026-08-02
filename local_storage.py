from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import keyring
from keyring.errors import KeyringError, PasswordDeleteError


KEYRING_SERVICE = "GeneConservado"
KEYRING_ACCOUNT = "credenciais-padrao"


def default_settings_path() -> Path:
    if sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    elif os.name == "nt":
        base = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming"))
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return base / "GeneConservado" / "configuracoes.json"


class LocalStateStore:
    def __init__(self, settings_path: Path | None = None, keyring_backend=None) -> None:
        self.settings_path = settings_path or default_settings_path()
        self.keyring = keyring_backend or keyring

    def load(self) -> tuple[dict[str, Any], dict[str, str]]:
        settings: dict[str, Any] = {}
        if self.settings_path.exists():
            loaded = json.loads(self.settings_path.read_text(encoding="utf-8"))
            if not isinstance(loaded, dict):
                raise ValueError("O arquivo local de configurações não contém um objeto JSON.")
            settings = loaded

        try:
            raw_credentials = self.keyring.get_password(KEYRING_SERVICE, KEYRING_ACCOUNT)
        except KeyringError as exc:
            raise RuntimeError("Não foi possível acessar o cofre seguro do sistema.") from exc

        credentials: dict[str, str] = {}
        if raw_credentials:
            loaded_credentials = json.loads(raw_credentials)
            if not isinstance(loaded_credentials, dict):
                raise ValueError("As credenciais locais estão em formato inválido.")
            credentials = {
                str(name): str(value)
                for name, value in loaded_credentials.items()
                if isinstance(name, str) and isinstance(value, (str, int, float, bool))
            }
        return settings, credentials

    def save(self, settings: dict[str, Any], credentials: dict[str, str]) -> None:
        credential_payload = json.dumps(credentials, ensure_ascii=False)
        try:
            self.keyring.set_password(KEYRING_SERVICE, KEYRING_ACCOUNT, credential_payload)
        except KeyringError as exc:
            raise RuntimeError("Não foi possível salvar as credenciais no cofre seguro.") from exc

        self.settings_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.settings_path.with_name(f".{self.settings_path.name}.tmp")
        temporary.write_text(
            json.dumps(settings, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        if os.name != "nt":
            temporary.chmod(0o600)
        os.replace(temporary, self.settings_path)

    def delete(self) -> None:
        self.settings_path.unlink(missing_ok=True)
        try:
            self.keyring.delete_password(KEYRING_SERVICE, KEYRING_ACCOUNT)
        except PasswordDeleteError:
            pass
        except KeyringError as exc:
            raise RuntimeError("Não foi possível apagar as credenciais do cofre seguro.") from exc
