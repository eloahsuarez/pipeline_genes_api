from pathlib import Path

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
