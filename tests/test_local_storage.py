import json
import os

import pytest

from local_storage import KEYRING_ACCOUNT, KEYRING_SERVICE, LocalStateStore


class MemoryKeyring:
    def __init__(self):
        self.values = {}

    def get_password(self, service, account):
        return self.values.get((service, account))

    def set_password(self, service, account, password):
        self.values[(service, account)] = password

    def delete_password(self, service, account):
        self.values.pop((service, account), None)


def test_local_state_round_trip_separates_settings_and_credentials(tmp_path):
    backend = MemoryKeyring()
    settings_path = tmp_path / "GeneConservado" / "configuracoes.json"
    store = LocalStateStore(settings_path=settings_path, keyring_backend=backend)

    settings = {"gene": "TP53", "max_records": 25}
    credentials = {"ncbi_api_key": "segredo", "idt_password": "senha"}
    store.save(settings, credentials)

    assert json.loads(settings_path.read_text(encoding="utf-8")) == settings
    assert "segredo" not in settings_path.read_text(encoding="utf-8")
    assert "senha" not in settings_path.read_text(encoding="utf-8")
    assert json.loads(backend.values[(KEYRING_SERVICE, KEYRING_ACCOUNT)]) == credentials
    assert store.load() == (settings, credentials)

    if os.name != "nt":
        assert settings_path.stat().st_mode & 0o777 == 0o600


def test_delete_removes_local_settings_and_keyring_entry(tmp_path):
    backend = MemoryKeyring()
    settings_path = tmp_path / "configuracoes.json"
    store = LocalStateStore(settings_path=settings_path, keyring_backend=backend)
    store.save({"gene": "BRCA1"}, {"idt_password": "senha"})

    store.delete()

    assert not settings_path.exists()
    assert backend.values == {}


def test_invalid_settings_file_is_rejected(tmp_path):
    settings_path = tmp_path / "configuracoes.json"
    settings_path.write_text("[]", encoding="utf-8")
    store = LocalStateStore(settings_path=settings_path, keyring_backend=MemoryKeyring())

    with pytest.raises(ValueError, match="objeto JSON"):
        store.load()
