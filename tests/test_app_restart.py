from pathlib import Path
from types import SimpleNamespace

import app


def test_restart_application_reopens_source_with_current_python(monkeypatch, tmp_path):
    executable = "/ambiente/python"
    target = tmp_path / "app.py"
    calls = []

    monkeypatch.setattr(app.sys, "executable", executable)
    monkeypatch.delattr(app.sys, "frozen", raising=False)
    monkeypatch.setattr(app.os, "execv", lambda path, args: calls.append((path, args)))

    app.restart_application(target)

    assert calls == [(executable, [executable, str(target.resolve())])]


def test_restart_application_reopens_packaged_executable(monkeypatch):
    executable = "/aplicativos/gene-conservado"
    calls = []

    monkeypatch.setattr(app.sys, "executable", executable)
    monkeypatch.setattr(app.sys, "frozen", True, raising=False)
    monkeypatch.setattr(app.os, "execv", lambda path, args: calls.append((path, args)))

    app.restart_application(Path("ignorado.py"))

    assert calls == [(executable, [executable])]


def test_exported_configuration_excludes_secrets():
    class Value:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

    fake_app = SimpleNamespace(
        vars={
            "gene": Value("TP53"),
            "ncbi_email": Value("pesquisador@example.com"),
            "ncbi_api_key": Value("chave-ncbi"),
            "idt_client_secret": Value("segredo-idt"),
            "idt_password": Value("senha-idt"),
        }
    )

    exported = app.GenePipelineApp.config_data(fake_app)

    assert exported == {"gene": "TP53", "ncbi_email": "pesquisador@example.com"}


def test_local_state_sends_credentials_only_to_keyring_payload():
    class Value:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

    fake_app = SimpleNamespace(
        vars={
            "gene": Value("TP53"),
            "ncbi_email": Value("pesquisador@example.com"),
            "ncbi_api_key": Value("chave-ncbi"),
            "idt_password": Value("senha-idt"),
        }
    )

    settings = app.GenePipelineApp._local_settings_data(fake_app)
    credentials = app.GenePipelineApp._credential_data(fake_app)

    assert settings == {"gene": "TP53"}
    assert credentials == {
        "ncbi_email": "pesquisador@example.com",
        "ncbi_api_key": "chave-ncbi",
        "idt_password": "senha-idt",
    }
