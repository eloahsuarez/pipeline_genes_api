from __future__ import annotations

import json
import os
import queue
import sys
import threading
import traceback
from datetime import datetime
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText

from bioinformatics import (
    PrimerDesignParams,
    find_conserved_regions,
    generate_exon_junction_primer_pairs,
    generate_primer_pairs,
    records_to_fasta,
)
from ebi_client import ClustalParams, EbiClustalClient
from exporter import (
    export_primers_xlsx,
    export_regions_csv,
    export_sequences_csv,
    write_json,
    write_text,
)
from idt_client import IdtClient, IdtConditions, IdtCredentials
from local_storage import LocalStateStore
from models import ConservedRegion, PrimerPair, SequenceRecord
from ncbi_client import NcbiClient, NcbiSearchParams
from version import APP_VERSION


APP_TITLE = "Gene Conservado — NCBI + Clustal Omega + IDT"
CREDENTIAL_VARIABLES = {
    "ncbi_email",
    "ncbi_api_key",
    "ebi_email",
    "idt_client_id",
    "idt_client_secret",
    "idt_username",
    "idt_password",
}
SENSITIVE_CONFIG_VARIABLES = {"ncbi_api_key", "idt_client_secret", "idt_password"}


def restart_application(app_path: Path | None = None) -> None:
    """Replace the current process with a fresh instance of the application."""
    if getattr(sys, "frozen", False):
        arguments = [sys.executable]
    else:
        target = (app_path or Path(__file__)).resolve()
        arguments = [sys.executable, str(target)]
    os.execv(sys.executable, arguments)


class GenePipelineApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title(f"{APP_TITLE} — v{APP_VERSION}")
        self.geometry("1260x820")
        self.minsize(1080, 700)

        self.records: list[SequenceRecord] = []
        self.query_used = ""
        self.ncbi_total = 0
        self.alignment = ""
        self.clustal_job_id = ""
        self.clustal_result_type = ""
        self.regions: list[ConservedRegion] = []
        self.pairs: list[PrimerPair] = []
        self.worker_queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self.busy = False
        self.local_state_store = LocalStateStore()
        self.local_persistence_enabled = True

        self.vars: dict[str, tk.Variable] = {}
        self._create_menu()
        self._create_ui()
        self._load_local_state()
        self.protocol("WM_DELETE_WINDOW", self.close_app)
        self.after(150, self._poll_worker_queue)

    def _var(self, name: str, value, kind: str = "str"):
        cls = {"str": tk.StringVar, "int": tk.IntVar, "float": tk.DoubleVar, "bool": tk.BooleanVar}[kind]
        variable = cls(value=value)
        self.vars[name] = variable
        return variable

    def _create_menu(self) -> None:
        menu = tk.Menu(self)
        file_menu = tk.Menu(menu, tearoff=False)
        file_menu.add_command(label="Salvar configuração…", command=self.save_config)
        file_menu.add_command(label="Carregar configuração…", command=self.load_config)
        file_menu.add_separator()
        file_menu.add_command(label="Exportar projeto…", command=self.export_project)
        file_menu.add_separator()
        file_menu.add_command(label="Sair", command=self.close_app)
        menu.add_cascade(label="Arquivo", menu=file_menu)

        local_menu = tk.Menu(menu, tearoff=False)
        local_menu.add_command(
            label="Salvar configurações e credenciais agora",
            command=lambda: self.save_local_state(show_confirmation=True),
        )
        local_menu.add_command(
            label="Apagar dados salvos deste computador…",
            command=self.delete_local_state,
        )
        menu.add_cascade(label="Dados locais", menu=local_menu)
        self.config(menu=menu)

    def _create_ui(self) -> None:
        root = ttk.Frame(self, padding=8)
        root.pack(fill="both", expand=True)
        header = ttk.Frame(root)
        ttk.Label(
            header,
            text=f"Versão {APP_VERSION}",
            font=("TkDefaultFont", 10, "bold"),
        ).pack(side="left")
        self.notebook = ttk.Notebook(root)

        self.tab_ncbi = ttk.Frame(self.notebook, padding=10)
        self.tab_sequences = ttk.Frame(self.notebook, padding=10)
        self.tab_clustal = ttk.Frame(self.notebook, padding=10)
        self.tab_design = ttk.Frame(self.notebook, padding=10)
        self.tab_idt = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(self.tab_ncbi, text="1. NCBI")
        self.notebook.add(self.tab_sequences, text="2. Sequências")
        self.notebook.add(self.tab_clustal, text="3. Clustal Omega")
        self.notebook.add(self.tab_design, text="4. Conservação e primers")
        self.notebook.add(self.tab_idt, text="5. IDT")

        self._build_ncbi_tab()
        self._build_sequences_tab()
        self._build_clustal_tab()
        self._build_design_tab()
        self._build_idt_tab()

        bottom = ttk.Frame(root)
        self.status_var = tk.StringVar(value="Pronto.")
        ttk.Label(bottom, textvariable=self.status_var).pack(side="left")
        ttk.Button(bottom, text="Exportar projeto", command=self.export_project).pack(side="right")
        ttk.Button(
            bottom,
            text="Atualizar ferramenta",
            command=self.restart_app,
        ).pack(side="right", padx=(0, 6))

        log_frame = ttk.LabelFrame(root, text="Log", padding=4)
        self.log_text = ScrolledText(log_frame, height=8, wrap="word")
        self.log_text.pack(fill="both", expand=True)

        # Reserve o rodapé antes da área expansível para que os controles não
        # sejam empurrados para fora da janela em telas menores.
        log_frame.pack(side="bottom", fill="x", expand=False, pady=(8, 0))
        bottom.pack(side="bottom", fill="x", pady=(8, 0))
        header.pack(side="top", fill="x", pady=(0, 6))
        self.notebook.pack(side="top", fill="both", expand=True)
        self.log(
            f"Aplicativo v{APP_VERSION} iniciado. Credenciais não são enviadas ao ChatGPT."
        )

    @staticmethod
    def _entry(parent, row, label, variable, width=34, show=None, column=0):
        ttk.Label(parent, text=label).grid(row=row, column=column, sticky="w", padx=4, pady=3)
        entry = ttk.Entry(parent, textvariable=variable, width=width, show=show)
        entry.grid(row=row, column=column + 1, sticky="ew", padx=4, pady=3)
        return entry

    def _build_ncbi_tab(self) -> None:
        frame = self.tab_ncbi
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(3, weight=1)

        self._entry(frame, 0, "E-mail NCBI *", self._var("ncbi_email", ""), column=0)
        self._entry(frame, 0, "API key NCBI", self._var("ncbi_api_key", ""), column=2)
        self._entry(frame, 1, "Gene *", self._var("gene", "TP53"), column=0)
        self._entry(frame, 1, "Organismo *", self._var("organism", "Homo sapiens"), column=2)
        self._entry(frame, 2, "Consulta Entrez adicional", self._var("extra_query", ""), width=55, column=0)

        ttk.Label(frame, text="Banco").grid(row=3, column=0, sticky="w", padx=4, pady=3)
        ttk.Combobox(frame, textvariable=self._var("database", "nuccore"), values=["nuccore"], state="readonly").grid(row=3, column=1, sticky="ew", padx=4, pady=3)
        ttk.Label(frame, text="Tipo de sequência").grid(row=3, column=2, sticky="w", padx=4, pady=3)
        ttk.Combobox(frame, textvariable=self._var("sequence_type", "mRNA"), values=["mRNA", "CDS", "genômico", "RNA", "qualquer"], state="readonly").grid(row=3, column=3, sticky="ew", padx=4, pady=3)

        self._entry(frame, 4, "Máximo de registros", self._var("max_records", 100, "int"), column=0)
        self._entry(frame, 4, "Comprimento mínimo", self._var("min_length", 100, "int"), column=2)
        self._entry(frame, 5, "Comprimento máximo", self._var("max_length", 100000, "int"), column=0)

        checks = ttk.LabelFrame(frame, text="Filtros", padding=8)
        checks.grid(row=6, column=0, columnspan=4, sticky="ew", padx=4, pady=8)
        ttk.Checkbutton(checks, text="Somente RefSeq", variable=self._var("refseq_only", True, "bool")).pack(side="left", padx=8)
        ttk.Checkbutton(checks, text="Excluir PREDICTED", variable=self._var("exclude_predicted", True, "bool")).pack(side="left", padx=8)
        ttk.Checkbutton(checks, text="Excluir partial", variable=self._var("exclude_partial", True, "bool")).pack(side="left", padx=8)
        ttk.Checkbutton(checks, text="Exigir anotação exata do gene", variable=self._var("require_gene_feature", True, "bool")).pack(side="left", padx=8)

        info = (
            "A busca inicial é automática, mas você poderá ativar ou desativar cada registro na aba Sequências. "
            "Para comparar isoformas, desmarque 'Exigir anotação exata' apenas quando necessário."
        )
        ttk.Label(frame, text=info, wraplength=1000).grid(row=7, column=0, columnspan=4, sticky="w", padx=4, pady=8)
        ttk.Button(frame, text="Buscar e baixar do NCBI", command=self.search_ncbi).grid(row=8, column=0, columnspan=4, pady=12)

    def _build_sequences_tab(self) -> None:
        top = ttk.Frame(self.tab_sequences)
        top.pack(fill="x")
        ttk.Button(top, text="Selecionar todas", command=lambda: self._set_all_records(True)).pack(side="left", padx=3)
        ttk.Button(top, text="Desmarcar todas", command=lambda: self._set_all_records(False)).pack(side="left", padx=3)
        ttk.Button(top, text="Exportar FASTA…", command=self.export_selected_fasta).pack(side="left", padx=3)
        ttk.Label(top, text="Dê duplo clique em uma linha para alternar Usar = Sim/Não.").pack(side="right")

        columns = ("use", "accession", "organism", "length", "genes", "definition")
        self.seq_tree = ttk.Treeview(self.tab_sequences, columns=columns, show="headings", selectmode="browse")
        headings = {"use": "Usar", "accession": "Acesso", "organism": "Organismo", "length": "bp", "genes": "Genes", "definition": "Definição"}
        widths = {"use": 55, "accession": 140, "organism": 160, "length": 75, "genes": 160, "definition": 600}
        for col in columns:
            self.seq_tree.heading(col, text=headings[col])
            self.seq_tree.column(col, width=widths[col], anchor="w")
        self.seq_tree.pack(fill="both", expand=True, pady=(8, 0))
        self.seq_tree.bind("<Double-1>", self.toggle_record)

    def _build_clustal_tab(self) -> None:
        frame = self.tab_clustal
        frame.columnconfigure(1, weight=1)
        self._entry(frame, 0, "E-mail EMBL-EBI *", self._var("ebi_email", ""), column=0)
        self._entry(frame, 1, "Título da tarefa", self._var("clustal_title", "Gene conserved regions"), column=0)
        self._entry(frame, 2, "Iterações", self._var("iterations", 0, "int"), column=0)
        self._entry(frame, 3, "Intervalo de consulta (s)", self._var("poll_seconds", 3, "int"), column=0)
        self._entry(frame, 4, "Tempo máximo (min)", self._var("timeout_minutes", 20, "int"), column=0)
        ttk.Checkbutton(frame, text="Remover gaps prévios (dealign)", variable=self._var("dealign", False, "bool")).grid(row=5, column=0, columnspan=2, sticky="w", padx=4, pady=3)
        ttk.Checkbutton(frame, text="Usar mBed", variable=self._var("mbed", True, "bool")).grid(row=6, column=0, columnspan=2, sticky="w", padx=4, pady=3)
        ttk.Button(frame, text="Alinhar sequências selecionadas", command=self.run_clustal).grid(row=7, column=0, columnspan=2, pady=12)
        self.clustal_info = tk.StringVar(value="Nenhum alinhamento executado.")
        ttk.Label(frame, textvariable=self.clustal_info, wraplength=900).grid(row=8, column=0, columnspan=2, sticky="w", padx=4, pady=8)
        self.alignment_preview = ScrolledText(frame, height=18, wrap="none")
        self.alignment_preview.grid(row=9, column=0, columnspan=2, sticky="nsew", padx=4, pady=4)
        frame.rowconfigure(9, weight=1)

    def _build_design_tab(self) -> None:
        outer = self.tab_design
        params = ttk.LabelFrame(outer, text="Parâmetros", padding=8)
        params.pack(fill="x")
        for col in (1, 3, 5):
            params.columnconfigure(col, weight=1)

        self._entry(params, 0, "Identidade mínima (%)", self._var("identity_pct", 100.0, "float"), column=0)
        self._entry(params, 0, "Cobertura mínima (%)", self._var("coverage_pct", 100.0, "float"), column=2)
        self._entry(params, 0, "Região mínima (bp)", self._var("min_region", 18, "int"), column=4)
        self._entry(params, 1, "Primer mínimo", self._var("primer_min_len", 18, "int"), column=0)
        self._entry(params, 1, "Primer máximo", self._var("primer_max_len", 25, "int"), column=2)
        self._entry(params, 1, "Top pares", self._var("top_pairs", 50, "int"), column=4)
        self._entry(params, 2, "GC mínimo (%)", self._var("gc_min", 35.0, "float"), column=0)
        self._entry(params, 2, "GC máximo (%)", self._var("gc_max", 65.0, "float"), column=2)
        self._entry(params, 2, "Tm alvo (°C)", self._var("tm_target", 60.0, "float"), column=4)
        self._entry(params, 3, "Tm mínimo (°C)", self._var("tm_min", 55.0, "float"), column=0)
        self._entry(params, 3, "Tm máximo (°C)", self._var("tm_max", 65.0, "float"), column=2)
        self._entry(params, 3, "Primer conc. local (nM)", self._var("local_primer_nm", 250.0, "float"), column=4)
        self._entry(params, 4, "Amplicon mínimo", self._var("amp_min", 80, "int"), column=0)
        self._entry(params, 4, "Amplicon máximo", self._var("amp_max", 250, "int"), column=2)
        self._entry(params, 4, "Amplicon alvo", self._var("amp_target", 150, "int"), column=4)

        junction = ttk.LabelFrame(outer, text="Junção éxon–éxon", padding=8)
        junction.pack(fill="x", pady=(8, 0))
        junction.columnconfigure(1, weight=1)
        ttk.Checkbutton(
            junction,
            text="Exigir que pelo menos um primer do par atravesse uma junção",
            variable=self._var("require_exon_junction", False, "bool"),
            command=self._update_junction_controls,
        ).grid(row=0, column=0, columnspan=6, sticky="w", padx=4, pady=3)
        ttk.Label(junction, text="Transcrito de referência").grid(
            row=1, column=0, sticky="w", padx=4, pady=3
        )
        self.junction_reference_combo = ttk.Combobox(
            junction,
            textvariable=self._var("junction_reference", ""),
            values=[],
            state="disabled",
        )
        self.junction_reference_combo.grid(row=1, column=1, sticky="ew", padx=4, pady=3)
        ttk.Label(junction, text="Ancoragem mínima 5′").grid(
            row=1, column=2, sticky="w", padx=4, pady=3
        )
        self.junction_min_5p_entry = ttk.Entry(
            junction, textvariable=self._var("junction_min_5p", 7, "int"), width=7
        )
        self.junction_min_5p_entry.grid(row=1, column=3, sticky="w", padx=4, pady=3)
        ttk.Label(junction, text="Ancoragem mínima 3′").grid(
            row=1, column=4, sticky="w", padx=4, pady=3
        )
        self.junction_min_3p_entry = ttk.Entry(
            junction, textvariable=self._var("junction_min_3p", 4, "int"), width=7
        )
        self.junction_min_3p_entry.grid(row=1, column=5, sticky="w", padx=4, pady=3)
        self.junction_info = tk.StringVar(
            value="Busque transcritos mRNA/RNA com anotação de éxons no NCBI."
        )
        ttk.Label(junction, textvariable=self.junction_info, wraplength=1000).grid(
            row=2, column=0, columnspan=6, sticky="w", padx=4, pady=3
        )
        self._update_junction_controls()

        actions = ttk.Frame(outer)
        actions.pack(fill="x", pady=8)
        ttk.Button(actions, text="Encontrar regiões conservadas", command=self.analyze_conservation).pack(side="left", padx=4)
        ttk.Button(actions, text="Gerar pares de primers", command=self.generate_primers).pack(side="left", padx=4)

        panes = ttk.Panedwindow(outer, orient="vertical")
        panes.pack(fill="both", expand=True)
        region_frame = ttk.LabelFrame(panes, text="Regiões conservadas")
        pair_frame = ttk.LabelFrame(panes, text="Pares candidatos")
        panes.add(region_frame, weight=1)
        panes.add(pair_frame, weight=2)

        rcols = ("start", "end", "length", "identity", "coverage", "sequence")
        self.region_tree = ttk.Treeview(region_frame, columns=rcols, show="headings")
        for col, text, width in [
            ("start", "Início", 80), ("end", "Fim", 80), ("length", "bp", 70),
            ("identity", "Identidade", 90), ("coverage", "Cobertura", 90), ("sequence", "Sequência", 650),
        ]:
            self.region_tree.heading(col, text=text)
            self.region_tree.column(col, width=width, anchor="w")
        self.region_tree.pack(fill="both", expand=True)

        pcols = ("rank", "score", "amp", "junction", "fwd", "ftm", "fgc", "rev", "rtm", "rgc")
        self.pair_tree = ttk.Treeview(pair_frame, columns=pcols, show="headings")
        for col, text, width in [
            ("rank", "Rank", 55), ("score", "Score", 70), ("amp", "Amplicon", 80),
            ("junction", "Junção", 180),
            ("fwd", "Forward", 210), ("ftm", "F Tm", 65), ("fgc", "F GC", 65),
            ("rev", "Reverse", 210), ("rtm", "R Tm", 65), ("rgc", "R GC", 65),
        ]:
            self.pair_tree.heading(col, text=text)
            self.pair_tree.column(col, width=width, anchor="w")
        self.pair_tree.pack(fill="both", expand=True)

    def _build_idt_tab(self) -> None:
        frame = self.tab_idt
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(3, weight=1)
        self._entry(frame, 0, "Client ID", self._var("idt_client_id", ""), column=0)
        self._entry(frame, 0, "Client Secret", self._var("idt_client_secret", ""), show="*", column=2)
        self._entry(frame, 1, "Usuário IDT", self._var("idt_username", ""), column=0)
        self._entry(frame, 1, "Senha IDT", self._var("idt_password", ""), show="*", column=2)
        self._entry(frame, 2, "Na+ (mM)", self._var("idt_na", 50.0, "float"), column=0)
        self._entry(frame, 2, "Mg2+ (mM)", self._var("idt_mg", 3.0, "float"), column=2)
        self._entry(frame, 3, "dNTP (mM)", self._var("idt_dntp", 0.8, "float"), column=0)
        self._entry(frame, 3, "Oligo (µM)", self._var("idt_oligo", 0.25, "float"), column=2)
        self._entry(frame, 4, "Temperatura de folding (°C)", self._var("idt_fold", 37.0, "float"), column=0)
        self._entry(frame, 4, "Quantidade de pares", self._var("idt_top_n", 10, "int"), column=2)

        actions = ttk.Frame(frame)
        actions.grid(row=5, column=0, columnspan=4, pady=10)
        ttk.Button(actions, text="Testar autenticação", command=self.test_idt).pack(side="left", padx=4)
        ttk.Button(actions, text="Analisar melhores pares na IDT", command=self.run_idt).pack(side="left", padx=4)
        ttk.Label(frame, text="A senha não é gravada no arquivo de configuração. Não cole credenciais nesta conversa.", wraplength=900).grid(row=6, column=0, columnspan=4, sticky="w", padx=4, pady=4)

        columns = ("rank", "f_tm", "r_tm", "f_hp", "r_hp", "f_sd", "r_sd", "hetero")
        self.idt_tree = ttk.Treeview(frame, columns=columns, show="headings")
        for col, text, width in [
            ("rank", "Rank", 55), ("f_tm", "F Tm IDT", 90), ("r_tm", "R Tm IDT", 90),
            ("f_hp", "F hairpin ΔG", 110), ("r_hp", "R hairpin ΔG", 110),
            ("f_sd", "F self ΔG", 100), ("r_sd", "R self ΔG", 100), ("hetero", "Hetero ΔG", 100),
        ]:
            self.idt_tree.heading(col, text=text)
            self.idt_tree.column(col, width=width, anchor="center")
        self.idt_tree.grid(row=7, column=0, columnspan=4, sticky="nsew", padx=4, pady=8)
        frame.rowconfigure(7, weight=1)

    def log(self, message: str) -> None:
        timestamp = datetime.now().strftime("%H:%M:%S")
        def append():
            self.log_text.insert("end", f"[{timestamp}] {message}\n")
            self.log_text.see("end")
        if threading.current_thread() is threading.main_thread():
            append()
        else:
            self.after(0, append)

    def _run_background(self, label: str, function, on_success) -> None:
        if self.busy:
            messagebox.showwarning("Processo em andamento", "Aguarde a operação atual terminar.")
            return
        self.busy = True
        self.status_var.set(label)

        def worker():
            try:
                result = function()
                self.worker_queue.put(("success", (on_success, result)))
            except Exception as exc:
                self.worker_queue.put(("error", (exc, traceback.format_exc())))
        threading.Thread(target=worker, daemon=True).start()

    def _poll_worker_queue(self) -> None:
        try:
            while True:
                kind, payload = self.worker_queue.get_nowait()
                self.busy = False
                self.status_var.set("Pronto.")
                if kind == "success":
                    callback, result = payload
                    callback(result)
                else:
                    exc, trace = payload
                    self.log(trace)
                    messagebox.showerror("Erro", str(exc))
        except queue.Empty:
            pass
        self.after(150, self._poll_worker_queue)

    def restart_app(self) -> None:
        if self.busy:
            messagebox.showwarning(
                "Processo em andamento",
                "Aguarde a operação atual terminar antes de atualizar a ferramenta.",
            )
            return
        if self.local_persistence_enabled and not self.save_local_state():
            return
        self.status_var.set("Atualizando e reiniciando…")
        self.log("Reiniciando para carregar as atualizações da pasta do projeto.")
        self.update_idletasks()
        self.after(100, self._restart_process)

    def _restart_process(self) -> None:
        self.destroy()
        restart_application()

    def _local_settings_data(self) -> dict:
        return {
            name: variable.get()
            for name, variable in self.vars.items()
            if name not in CREDENTIAL_VARIABLES
        }

    def _credential_data(self) -> dict[str, str]:
        return {
            name: str(self.vars[name].get())
            for name in CREDENTIAL_VARIABLES
            if name in self.vars
        }

    def _load_local_state(self) -> None:
        try:
            settings, credentials = self.local_state_store.load()
        except Exception as exc:
            self.log(f"Não foi possível carregar os dados locais: {exc}")
            messagebox.showwarning("Dados locais", str(exc))
            return

        for name, value in {**settings, **credentials}.items():
            if name in self.vars:
                self.vars[name].set(value)
        self._refresh_junction_references()
        if settings or credentials:
            self.log("Configurações e credenciais locais carregadas com segurança.")

    def save_local_state(self, show_confirmation: bool = False) -> bool:
        try:
            self.local_state_store.save(
                self._local_settings_data(),
                self._credential_data(),
            )
        except Exception as exc:
            self.log(f"Não foi possível salvar os dados locais: {exc}")
            messagebox.showerror("Dados locais", str(exc))
            return False

        self.local_persistence_enabled = True
        self.log("Configurações salvas localmente e credenciais guardadas no cofre do sistema.")
        if show_confirmation:
            messagebox.showinfo(
                "Dados locais",
                "Configurações e credenciais foram salvas somente neste computador.",
            )
        return True

    def delete_local_state(self) -> None:
        confirmed = messagebox.askyesno(
            "Apagar dados locais",
            "Apagar configurações e credenciais salvas neste computador?",
        )
        if not confirmed:
            return
        try:
            self.local_state_store.delete()
        except Exception as exc:
            self.log(f"Não foi possível apagar os dados locais: {exc}")
            messagebox.showerror("Dados locais", str(exc))
            return
        for name in CREDENTIAL_VARIABLES:
            if name in self.vars:
                self.vars[name].set("")
        self.local_persistence_enabled = False
        self.log("Configurações e credenciais locais apagadas.")
        messagebox.showinfo("Dados locais", "Os dados salvos neste computador foram apagados.")

    def close_app(self) -> None:
        if self.local_persistence_enabled and not self.save_local_state():
            close_without_saving = messagebox.askyesno(
                "Fechar sem salvar",
                "Não foi possível salvar os dados locais. Deseja fechar mesmo assim?",
            )
            if not close_without_saving:
                return
        self.destroy()

    def _invalidate_alignment_results(self, announce: bool = False) -> None:
        had_results = bool(self.alignment or self.regions or self.pairs)
        self.alignment = ""
        self.clustal_job_id = ""
        self.clustal_result_type = ""
        self.regions = []
        self.pairs = []
        self.clustal_info.set("Nenhum alinhamento executado.")
        self.alignment_preview.delete("1.0", "end")
        self.region_tree.delete(*self.region_tree.get_children())
        self.pair_tree.delete(*self.pair_tree.get_children())
        self.idt_tree.delete(*self.idt_tree.get_children())
        if announce and had_results:
            self.log("A seleção de sequências mudou; execute o Clustal Omega novamente.")

    def search_ncbi(self) -> None:
        params = NcbiSearchParams(
            gene=self.vars["gene"].get(), organism=self.vars["organism"].get(),
            database=self.vars["database"].get(), max_records=self.vars["max_records"].get(),
            sequence_type=self.vars["sequence_type"].get(), refseq_only=self.vars["refseq_only"].get(),
            exclude_predicted=self.vars["exclude_predicted"].get(), exclude_partial=self.vars["exclude_partial"].get(),
            require_gene_feature=self.vars["require_gene_feature"].get(), min_length=self.vars["min_length"].get(),
            max_length=self.vars["max_length"].get(), extra_query=self.vars["extra_query"].get(),
        )
        ncbi_email = self.vars["ncbi_email"].get()
        ncbi_api_key = self.vars["ncbi_api_key"].get()

        def task():
            client = NcbiClient(ncbi_email, ncbi_api_key, log=self.log)
            return client.search_and_fetch(params)

        def done(result):
            self.records, self.query_used, self.ncbi_total = result
            self._invalidate_alignment_results()
            self._refresh_sequences()
            self.notebook.select(self.tab_sequences)
            self.log(f"Busca concluída: {len(self.records)} registros utilizáveis de {self.ncbi_total} encontrados.")

        self._run_background("Consultando o NCBI…", task, done)

    def _refresh_sequences(self) -> None:
        self.seq_tree.delete(*self.seq_tree.get_children())
        for index, record in enumerate(self.records):
            self.seq_tree.insert("", "end", iid=str(index), values=(
                "Sim" if record.selected else "Não", record.accession, record.organism,
                record.length, "; ".join(record.genes), record.definition,
            ))
        self._refresh_junction_references()

    @staticmethod
    def _record_has_exon_junctions(record: SequenceRecord) -> bool:
        exons = record.exons
        molecule_type = record.molecule_type.casefold()
        if "rna" not in molecule_type or "genomic" in molecule_type:
            return False
        return (
            len(exons) >= 2
            and exons == sorted(exons, key=lambda exon: (exon.start, exon.end))
            and all(1 <= exon.start <= exon.end <= record.length for exon in exons)
            and all(
                left.end + 1 == right.start
                for left, right in zip(exons, exons[1:])
            )
        )

    def _refresh_junction_references(self) -> None:
        eligible = [
            record.accession
            for record in self.records
            if record.selected and self._record_has_exon_junctions(record)
        ]
        self.junction_reference_combo.configure(values=eligible)
        if eligible and self.vars["junction_reference"].get() not in eligible:
            self.vars["junction_reference"].set(eligible[0] if eligible else "")
        elif self.records and not eligible:
            self.vars["junction_reference"].set("")
        if eligible:
            self.junction_info.set(
                f"{len(eligible)} transcrito(s) selecionado(s) possuem junções anotadas; "
                "escolha a isoforma de referência."
            )
        else:
            self.junction_info.set(
                "Nenhum transcrito selecionado possui pelo menos dois éxons contíguos anotados."
            )
        self._update_junction_controls()

    def _update_junction_controls(self) -> None:
        enabled = self.vars.get("require_exon_junction") is not None and self.vars[
            "require_exon_junction"
        ].get()
        self.junction_reference_combo.configure(state="readonly" if enabled else "disabled")
        entry_state = "normal" if enabled else "disabled"
        self.junction_min_5p_entry.configure(state=entry_state)
        self.junction_min_3p_entry.configure(state=entry_state)

    def toggle_record(self, _event=None) -> None:
        if self.busy:
            messagebox.showwarning(
                "Processo em andamento", "Aguarde a operação atual antes de alterar a seleção."
            )
            return
        selected = self.seq_tree.selection()
        if not selected:
            return
        index = int(selected[0])
        self.records[index].selected = not self.records[index].selected
        self._invalidate_alignment_results(announce=True)
        self._refresh_sequences()

    def _set_all_records(self, value: bool) -> None:
        if self.busy:
            messagebox.showwarning(
                "Processo em andamento", "Aguarde a operação atual antes de alterar a seleção."
            )
            return
        changed = any(record.selected != value for record in self.records)
        for record in self.records:
            record.selected = value
        if changed:
            self._invalidate_alignment_results(announce=True)
        self._refresh_sequences()

    def export_selected_fasta(self) -> None:
        fasta = records_to_fasta(self.records)
        if not fasta:
            messagebox.showwarning("Sem sequências", "Nenhuma sequência está selecionada.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".fasta", filetypes=[("FASTA", "*.fasta *.fa"), ("Todos", "*.*")])
        if path:
            Path(path).write_text(fasta, encoding="utf-8")

    def run_clustal(self) -> None:
        fasta = records_to_fasta(self.records)
        if fasta.count(">") < 3:
            messagebox.showwarning("Poucas sequências", "Selecione pelo menos três sequências.")
            return

        params = ClustalParams(
            email=self.vars["ebi_email"].get(), title=self.vars["clustal_title"].get(),
            iterations=self.vars["iterations"].get(), dealign=self.vars["dealign"].get(),
            mbed=self.vars["mbed"].get(), poll_seconds=self.vars["poll_seconds"].get(),
            timeout_minutes=self.vars["timeout_minutes"].get(),
        )

        def task():
            return EbiClustalClient(log=self.log).run(fasta, params)

        def done(result):
            job_id, result_type, alignment = result
            self._invalidate_alignment_results()
            self.clustal_job_id = job_id
            self.clustal_result_type = result_type
            self.alignment = alignment
            self.clustal_info.set(f"Job: {self.clustal_job_id} | Resultado: {self.clustal_result_type} | {len(self.alignment):,} caracteres")
            self.alignment_preview.delete("1.0", "end")
            self.alignment_preview.insert("1.0", self.alignment[:100000])
            self.notebook.select(self.tab_design)

        self._run_background("Executando Clustal Omega…", task, done)

    def _design_params(self) -> PrimerDesignParams:
        return PrimerDesignParams(
            identity_threshold=self.vars["identity_pct"].get() / 100.0,
            coverage_threshold=self.vars["coverage_pct"].get() / 100.0,
            min_primer_len=self.vars["primer_min_len"].get(), max_primer_len=self.vars["primer_max_len"].get(),
            min_gc=self.vars["gc_min"].get(), max_gc=self.vars["gc_max"].get(),
            min_tm=self.vars["tm_min"].get(), max_tm=self.vars["tm_max"].get(), target_tm=self.vars["tm_target"].get(),
            min_amplicon=self.vars["amp_min"].get(), max_amplicon=self.vars["amp_max"].get(),
            target_amplicon=self.vars["amp_target"].get(), top_pairs=self.vars["top_pairs"].get(),
            na_mm=self.vars["idt_na"].get(), mg_mm=self.vars["idt_mg"].get(), dntp_mm=self.vars["idt_dntp"].get(),
            primer_conc_nm=self.vars["local_primer_nm"].get(),
        )

    def analyze_conservation(self) -> None:
        if not self.alignment:
            messagebox.showwarning("Sem alinhamento", "Execute o Clustal Omega primeiro.")
            return
        self.regions = find_conserved_regions(
            self.alignment, self.vars["identity_pct"].get() / 100.0,
            self.vars["coverage_pct"].get() / 100.0, self.vars["min_region"].get(),
        )
        self.region_tree.delete(*self.region_tree.get_children())
        for region in self.regions:
            self.region_tree.insert("", "end", values=(
                region.consensus_start, region.consensus_end, region.length,
                f"{region.mean_identity * 100:.1f}%", f"{region.mean_coverage * 100:.1f}%", region.sequence,
            ))
        self.log(f"Foram encontradas {len(self.regions)} regiões conservadas.")

    def generate_primers(self) -> None:
        if not self.alignment:
            messagebox.showwarning("Sem alinhamento", "Execute o Clustal Omega primeiro.")
            return

        alignment = self.alignment
        design_params = self._design_params()
        require_junction = self.vars["require_exon_junction"].get()
        reference_record: SequenceRecord | None = None
        min_5_prime_match = self.vars["junction_min_5p"].get()
        min_3_prime_match = self.vars["junction_min_3p"].get()
        if require_junction:
            reference_id = self.vars["junction_reference"].get()
            reference_record = next(
                (
                    record
                    for record in self.records
                    if record.selected and record.accession == reference_id
                ),
                None,
            )
            if reference_record is None or not self._record_has_exon_junctions(reference_record):
                messagebox.showwarning(
                    "Sem referência anotada",
                    "Escolha um transcrito mRNA/RNA selecionado com pelo menos dois éxons anotados.",
                )
                return
            if min_5_prime_match < 1 or min_3_prime_match < 1:
                messagebox.showwarning(
                    "Ancoragem inválida", "As ancoragens mínimas 5′ e 3′ devem ser maiores que zero."
                )
                return

        def task():
            if require_junction and reference_record is not None:
                return generate_exon_junction_primer_pairs(
                    alignment,
                    design_params,
                    reference_record.accession,
                    reference_record.exons,
                    min_5_prime_match=min_5_prime_match,
                    min_3_prime_match=min_3_prime_match,
                )
            return generate_primer_pairs(alignment, design_params)

        def done(result):
            self.pairs = result
            self._refresh_pairs()
            self.idt_tree.delete(*self.idt_tree.get_children())
            if require_junction and reference_record is not None:
                self.log(
                    f"Foram gerados {len(self.pairs)} pares com primer em junção "
                    f"para {reference_record.accession}."
                )
            else:
                self.log(f"Foram gerados {len(self.pairs)} pares candidatos.")
            self.notebook.select(self.tab_idt)

        self._run_background("Gerando candidatos de primers…", task, done)

    def _refresh_pairs(self) -> None:
        self.pair_tree.delete(*self.pair_tree.get_children())
        for pair in self.pairs:
            junctions: list[str] = []
            for candidate in (pair.forward, pair.reverse):
                for match in candidate.junctions:
                    exon_names = (
                        f"E{match.left_exon_number}–E{match.right_exon_number}"
                        if match.left_exon_number and match.right_exon_number
                        else f"posição {match.junction_position}"
                    )
                    junctions.append(
                        f"{candidate.orientation}: {exon_names} "
                        f"(5′ {match.primer_5_prime_bases}/3′ {match.primer_3_prime_bases})"
                    )
            self.pair_tree.insert("", "end", values=(
                pair.rank, f"{pair.score:.2f}", pair.amplicon_length,
                "; ".join(junctions) or "—",
                pair.forward.sequence, f"{pair.forward.tm_c:.2f}", f"{pair.forward.gc_percent:.1f}%",
                pair.reverse.sequence, f"{pair.reverse.tm_c:.2f}", f"{pair.reverse.gc_percent:.1f}%",
            ))

    def _idt_values(self):
        credentials = IdtCredentials(
            client_id=self.vars["idt_client_id"].get(), client_secret=self.vars["idt_client_secret"].get(),
            username=self.vars["idt_username"].get(), password=self.vars["idt_password"].get(),
        )
        conditions = IdtConditions(
            na_mm=self.vars["idt_na"].get(), mg_mm=self.vars["idt_mg"].get(),
            dntp_mm=self.vars["idt_dntp"].get(), oligo_um=self.vars["idt_oligo"].get(),
            folding_temp_c=self.vars["idt_fold"].get(),
        )
        return credentials, conditions

    def test_idt(self) -> None:
        credentials, _ = self._idt_values()

        def task():
            return IdtClient(credentials, log=self.log).authenticate()

        def done(_):
            messagebox.showinfo("IDT", "Autenticação concluída.")
        self._run_background("Autenticando na IDT…", task, done)

    def run_idt(self) -> None:
        if not self.pairs:
            messagebox.showwarning("Sem primers", "Gere os pares de primers primeiro.")
            return
        count = min(self.vars["idt_top_n"].get(), len(self.pairs))
        credentials, conditions = self._idt_values()
        pairs_to_analyze = self.pairs[:count]

        def task():
            client = IdtClient(credentials, log=self.log)
            client.authenticate()
            for pair in pairs_to_analyze:
                pair.idt = client.analyze_pair(pair.forward.sequence, pair.reverse.sequence, conditions)
            return count

        def done(processed):
            self._refresh_idt()
            self.log(f"IDT: {processed} pares analisados.")

        self._run_background("Analisando primers na IDT…", task, done)

    @staticmethod
    def _value(data, *keys):
        current = data
        for key in keys:
            if not isinstance(current, dict):
                return ""
            current = current.get(key, "")
        return current

    def _refresh_idt(self) -> None:
        self.idt_tree.delete(*self.idt_tree.get_children())
        for pair in self.pairs:
            if not pair.idt:
                continue
            f = pair.idt.get("forward", {})
            r = pair.idt.get("reverse", {})
            self.idt_tree.insert("", "end", values=(
                pair.rank,
                self._value(f, "analysis", "MeltTemp"), self._value(r, "analysis", "MeltTemp"),
                self._value(f, "strongest_hairpin", "deltaG"), self._value(r, "strongest_hairpin", "deltaG"),
                self._value(f, "strongest_self_dimer", "DeltaG"), self._value(r, "strongest_self_dimer", "DeltaG"),
                self._value(pair.idt, "strongest_hetero_dimer", "DeltaG"),
            ))

    def config_data(self) -> dict:
        return {
            name: variable.get()
            for name, variable in self.vars.items()
            if name not in SENSITIVE_CONFIG_VARIABLES
        }

    def save_config(self) -> None:
        path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON", "*.json")])
        if path:
            Path(path).write_text(json.dumps(self.config_data(), ensure_ascii=False, indent=2), encoding="utf-8")

    def load_config(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("JSON", "*.json"), ("Todos", "*.*")])
        if not path:
            return
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        for name, value in data.items():
            if name in self.vars:
                self.vars[name].set(value)
        self._refresh_junction_references()

    def export_project(self) -> None:
        folder = filedialog.askdirectory(title="Escolha a pasta de saída")
        if not folder:
            return
        target = Path(folder) / f"gene_pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        target.mkdir(parents=True, exist_ok=True)
        write_json(target / "configuracao.json", self.config_data())
        write_json(target / "resumo.json", {
            "consulta_ncbi": self.query_used, "total_ncbi": self.ncbi_total,
            "clustal_job_id": self.clustal_job_id, "clustal_result_type": self.clustal_result_type,
        })
        if self.records:
            export_sequences_csv(target / "sequencias.csv", self.records)
            write_text(target / "sequencias_selecionadas.fasta", records_to_fasta(self.records))
        if self.alignment:
            write_text(target / "alinhamento.fasta", self.alignment)
        if self.regions:
            export_regions_csv(target / "regioes_conservadas.csv", self.regions)
        if self.pairs:
            export_primers_xlsx(target / "primers.xlsx", self.pairs)
            write_json(target / "primers.json", self.pairs)
        self.log(f"Projeto exportado para {target}")
        messagebox.showinfo("Exportação", f"Arquivos salvos em:\n{target}")


if __name__ == "__main__":
    app = GenePipelineApp()
    app.mainloop()
