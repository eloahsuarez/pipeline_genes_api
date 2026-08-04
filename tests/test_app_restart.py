from pathlib import Path
from types import SimpleNamespace

import app
import pytest
from models import SequenceRecord
from ncbi_blast import REFSEQ_MRNA_DATABASE, REFSEQ_SELECT_DATABASE


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


def test_selected_idt_pair_uses_the_selected_ranking_row():
    fake_app = SimpleNamespace(
        idt_tree=SimpleNamespace(selection=lambda: ("2",)),
        pairs=[
            SimpleNamespace(rank=1, idt={"forward": {}}),
            SimpleNamespace(rank=2, idt={"forward": {}}),
        ],
    )

    selected = app.GenePipelineApp._selected_idt_pair(fake_app)

    assert selected.rank == 2


def test_ncbi_target_accessions_include_selected_records_and_pair_reference():
    fake_app = SimpleNamespace(
        records=[
            SimpleNamespace(accession="NM_000546.6", selected=True),
            SimpleNamespace(accession="NM_001126112.3", selected=False),
        ],
        pairs=[
            SimpleNamespace(reference_accession="NM_000546.6"),
            SimpleNamespace(reference_accession="NM_001276760.2"),
        ],
    )

    accessions = app.GenePipelineApp._ncbi_target_accessions(fake_app)

    assert accessions == ("NM_000546.6", "NM_001276760.2")


def test_ncbi_classification_labels_are_user_facing():
    assert app.GenePipelineApp._ncbi_classification_label("target") == "Alvo"
    assert app.GenePipelineApp._ncbi_classification_label("off_target") == "Outro gene"


def test_ncbi_remote_database_options_keep_exact_aliases_and_friendly_labels():
    assert app.DEFAULT_NCBI_SPECIFICITY_DATABASE == REFSEQ_MRNA_DATABASE
    assert app.DEFAULT_NCBI_SPECIFICITY_TOP_PAIRS == 1
    assert app.DEFAULT_NCBI_SPECIFICITY_HITLIST_SIZE == 50_000
    assert app.DEFAULT_NCBI_SPECIFICITY_MAX_AMPLICON == 4_000
    assert app.DEFAULT_NCBI_SPECIFICITY_MAX_ESTIMATED_MISMATCHES == 6
    assert app.NCBI_SPECIFICITY_THREE_PRIME_WINDOW == 5
    assert app.NCBI_SPECIFICITY_WORD_SIZE == 7
    assert app.NCBI_SPECIFICITY_DATABASE_LABELS == {
        REFSEQ_MRNA_DATABASE: "RefSeq mRNA (refseq_mrna)",
        REFSEQ_SELECT_DATABASE: "RefSeq Select RNA (refseq_select_rna)",
    }


def test_old_remote_profile_settings_are_migrated_as_one_profile():
    old_settings = {
        "ncbi_specificity_top_n": 10,
        "ncbi_specificity_max_hits": 100,
        "ncbi_specificity_max_amplicon": 1000,
        "ncbi_specificity_identity_pct": 80.0,
    }

    migrated = app.migrate_ncbi_specificity_profile_settings(old_settings)

    assert migrated["ncbi_specificity_top_n"] == 1
    assert migrated["ncbi_specificity_max_hits"] == 50_000
    assert migrated["ncbi_specificity_max_amplicon"] == 4_000
    assert migrated["ncbi_specificity_max_estimated_mismatches"] == 6
    assert old_settings["ncbi_specificity_top_n"] == 10


def test_new_remote_profile_settings_keep_user_values():
    settings = {
        "ncbi_specificity_top_n": 3,
        "ncbi_specificity_max_hits": 12_000,
        "ncbi_specificity_max_amplicon": 2_500,
        "ncbi_specificity_max_estimated_mismatches": 4,
    }

    assert app.migrate_ncbi_specificity_profile_settings(settings) == settings


def test_load_local_state_migrates_old_remote_profile():
    class Variable:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

        def set(self, value):
            self.value = value

    variables = {
        "ncbi_specificity_top_n": Variable(1),
        "ncbi_specificity_max_hits": Variable(50_000),
        "ncbi_specificity_max_amplicon": Variable(4_000),
        "ncbi_specificity_max_estimated_mismatches": Variable(6),
    }
    fake_app = SimpleNamespace(
        local_state_store=SimpleNamespace(
            load=lambda: (
                {
                    "ncbi_specificity_top_n": 10,
                    "ncbi_specificity_max_hits": 100,
                    "ncbi_specificity_max_amplicon": 1000,
                },
                {},
            )
        ),
        vars=variables,
        _refresh_junction_references=lambda: None,
        log=lambda _message: None,
    )

    app.GenePipelineApp._load_local_state(fake_app)

    assert variables["ncbi_specificity_top_n"].get() == 1
    assert variables["ncbi_specificity_max_hits"].get() == 50_000
    assert variables["ncbi_specificity_max_amplicon"].get() == 4_000
    assert variables["ncbi_specificity_max_estimated_mismatches"].get() == 6


def test_load_exported_config_migrates_old_remote_profile(monkeypatch, tmp_path):
    class Variable:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

        def set(self, value):
            self.value = value

    path = tmp_path / "configuracao-v1.7.json"
    path.write_text(
        '{"ncbi_specificity_top_n": 10, "ncbi_specificity_max_hits": 100, '
        '"ncbi_specificity_max_amplicon": 1000}',
        encoding="utf-8",
    )
    variables = {
        "ncbi_specificity_top_n": Variable(1),
        "ncbi_specificity_max_hits": Variable(50_000),
        "ncbi_specificity_max_amplicon": Variable(4_000),
        "ncbi_specificity_max_estimated_mismatches": Variable(6),
    }
    fake_app = SimpleNamespace(
        vars=variables,
        _refresh_junction_references=lambda: None,
    )
    monkeypatch.setattr(
        app.filedialog,
        "askopenfilename",
        lambda **_kwargs: str(path),
    )

    app.GenePipelineApp.load_config(fake_app)

    assert variables["ncbi_specificity_top_n"].get() == 1
    assert variables["ncbi_specificity_max_hits"].get() == 50_000
    assert variables["ncbi_specificity_max_amplicon"].get() == 4_000
    assert variables["ncbi_specificity_max_estimated_mismatches"].get() == 6


def test_exported_configuration_excludes_secrets():
    class Value:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

    fake_app = SimpleNamespace(
        vars={
            "gene": Value("TP53"),
            "ncbi_specificity_database": Value(REFSEQ_MRNA_DATABASE),
            "ncbi_specificity_max_estimated_mismatches": Value(6),
            "ncbi_email": Value("pesquisador@example.com"),
            "ncbi_api_key": Value("chave-ncbi"),
            "idt_client_secret": Value("segredo-idt"),
            "idt_password": Value("senha-idt"),
        }
    )

    exported = app.GenePipelineApp.config_data(fake_app)

    assert exported == {
        "gene": "TP53",
        "ncbi_specificity_database": REFSEQ_MRNA_DATABASE,
        "ncbi_specificity_max_estimated_mismatches": 6,
        "ncbi_email": "pesquisador@example.com",
    }


def test_local_state_sends_credentials_only_to_keyring_payload():
    class Value:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

    fake_app = SimpleNamespace(
        vars={
            "gene": Value("TP53"),
            "ncbi_specificity_database": Value(REFSEQ_MRNA_DATABASE),
            "ncbi_specificity_max_estimated_mismatches": Value(6),
            "ncbi_email": Value("pesquisador@example.com"),
            "ncbi_api_key": Value("chave-ncbi"),
            "idt_password": Value("senha-idt"),
        }
    )

    settings = app.GenePipelineApp._local_settings_data(fake_app)
    credentials = app.GenePipelineApp._credential_data(fake_app)

    assert settings == {
        "gene": "TP53",
        "ncbi_specificity_database": REFSEQ_MRNA_DATABASE,
        "ncbi_specificity_max_estimated_mismatches": 6,
    }
    assert credentials == {
        "ncbi_email": "pesquisador@example.com",
        "ncbi_api_key": "chave-ncbi",
        "idt_password": "senha-idt",
    }


def test_alignment_strategy_uses_pairwise_for_two_and_clustal_for_three():
    two = [SimpleNamespace(selected=True), SimpleNamespace(selected=True)]
    three = [*two, SimpleNamespace(selected=True)]

    assert app.alignment_strategy(two) == "pairwise"
    assert app.alignment_strategy(three) == "clustal"


def test_alignment_strategy_counts_only_selected_records():
    records = [
        SimpleNamespace(selected=True),
        SimpleNamespace(selected=False),
        SimpleNamespace(selected=False),
    ]

    with pytest.raises(ValueError, match="pelo menos duas"):
        app.alignment_strategy(records)


def test_primer_design_source_accepts_one_sequence_without_alignment():
    record = SimpleNamespace(selected=True, sequence="AACCGG", accession="SEQ1")

    mode, payload, selected = app.primer_design_source("", [record])

    assert mode == "single"
    assert payload == "AACCGG"
    assert selected is record


def test_primer_design_source_prefers_an_existing_alignment():
    records = [SimpleNamespace(selected=True), SimpleNamespace(selected=True)]

    assert app.primer_design_source(">a\nAA\n>b\nAA\n", records) == (
        "alignment",
        ">a\nAA\n>b\nAA\n",
        None,
    )


@pytest.mark.parametrize(
    "records, message",
    [
        ([], "Selecione uma sequência"),
        (
            [SimpleNamespace(selected=True), SimpleNamespace(selected=True)],
            "Alinhe as sequências",
        ),
    ],
)
def test_primer_design_source_rejects_ambiguous_selection(records, message):
    with pytest.raises(ValueError, match=message):
        app.primer_design_source("", records)


def _sequence_record(accession: str, selected: bool = True) -> SequenceRecord:
    sequence = "ACGTACGT"
    return SequenceRecord(
        uid=accession,
        accession=accession,
        definition=f"Registro {accession}",
        organism="Organismo teste",
        sequence=sequence,
        length=len(sequence),
        selected=selected,
    )


def test_run_clustal_routes_exactly_two_sequences_to_local_pairwise(monkeypatch):
    scheduled = {}
    fake_app = SimpleNamespace(
        records=[_sequence_record("A"), _sequence_record("B")],
        vars={},
        log=lambda _message: None,
        _run_background=lambda label, task, done: scheduled.update(
            label=label, task=task, done=done
        ),
    )
    monkeypatch.setattr(app, "pairwise_align_fasta", lambda _fasta: ">A\nACGT\n>B\nACGT\n")

    app.GenePipelineApp.run_clustal(fake_app)

    assert scheduled["label"] == "Alinhando duas sequências localmente…"
    assert scheduled["task"]() == (
        "",
        "pareado-local-fasta",
        ">A\nACGT\n>B\nACGT\n",
    )


def test_run_clustal_keeps_three_sequences_on_ebi_client(monkeypatch):
    class Value:
        def __init__(self, value):
            self.value = value

        def get(self):
            return self.value

    scheduled = {}
    calls = []

    class FakeClient:
        def __init__(self, log):
            self.log = log

        def run(self, fasta, params):
            calls.append((fasta, params))
            return "job-1", "aln-fasta", "alinhamento"

    fake_app = SimpleNamespace(
        records=[
            _sequence_record("A"),
            _sequence_record("B"),
            _sequence_record("C"),
        ],
        vars={
            "ebi_email": Value("pesquisador@example.com"),
            "clustal_title": Value("Teste"),
            "iterations": Value(0),
            "dealign": Value(False),
            "mbed": Value(True),
            "poll_seconds": Value(3),
            "timeout_minutes": Value(20),
        },
        log=lambda _message: None,
        _run_background=lambda label, task, done: scheduled.update(
            label=label, task=task, done=done
        ),
    )
    monkeypatch.setattr(app, "EbiClustalClient", FakeClient)

    app.GenePipelineApp.run_clustal(fake_app)

    assert scheduled["label"] == "Executando Clustal Omega…"
    assert scheduled["task"]() == ("job-1", "aln-fasta", "alinhamento")
    assert len(calls) == 1
    assert calls[0][0].count(">") == 3


def test_close_app_waits_for_background_process(monkeypatch):
    warnings = []
    destroyed = []
    fake_app = SimpleNamespace(
        busy=True,
        destroy=lambda: destroyed.append(True),
    )
    monkeypatch.setattr(
        app.messagebox,
        "showwarning",
        lambda title, message: warnings.append((title, message)),
    )

    app.GenePipelineApp.close_app(fake_app)

    assert destroyed == []
    assert warnings and warnings[0][0] == "Processo em andamento"
