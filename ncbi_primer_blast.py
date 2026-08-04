import webbrowser
import urllib.parse
import requests
import re
import time
from typing import Callable, Optional
from bs4 import BeautifulSoup

PRIMER_BLAST_URL = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi"

def open_primer_blast(
    forward_primer: str,
    reverse_primer: str,
    exon_junction_span: str = "0",
    database: str = "refseq_mrna",
    organism: str = "Homo sapiens (taxid:9606)"
) -> None:
    """
    Constrói a requisição para o NCBI Primer-BLAST e abre o navegador
    na página de resultados. (Usado pelo botão Abrir no Navegador)
    """
    params = _build_params(forward_primer, reverse_primer, exon_junction_span, database, organism)

    try:
        response = requests.post(PRIMER_BLAST_URL, data=params, timeout=15)
        response.raise_for_status()
        
        match = re.search(r'RID=([A-Z0-9]+)', response.text)
        if not match:
            match = re.search(r'name="RID"\s+value="([A-Z0-9]+)"', response.text, re.IGNORECASE)
            
        if match:
            rid = match.group(1)
            result_url = f"{PRIMER_BLAST_URL}?RID={rid}"
            webbrowser.open(result_url)
        else:
            match_job = re.search(r'name="job_key"\s+value="([^"]+)"', response.text, re.IGNORECASE)
            match_ctg = re.search(r'name="ctg_time"\s+value="([^"]+)"', response.text, re.IGNORECASE)
            if match_job:
                # Polling parameters check page
                query_string = urllib.parse.urlencode({
                    "CMD": "request",
                    "job_key": match_job.group(1),
                    "ctg_time": match_ctg.group(1) if match_ctg else "",
                    "CheckStatus": "Check"
                })
                webbrowser.open(f"{PRIMER_BLAST_URL}?{query_string}")
            else:
                query_string = urllib.parse.urlencode(params)
                full_url = f"{PRIMER_BLAST_URL}?{query_string}"
                webbrowser.open(full_url)
    except Exception as e:
        print(f"Erro ao submeter para o Primer-BLAST: {e}")
        query_string = urllib.parse.urlencode(params)
        full_url = f"{PRIMER_BLAST_URL}?{query_string}"
        webbrowser.open(full_url)


def run_primer_blast_sync(
    forward_primer: str,
    reverse_primer: str,
    exon_junction_span: str = "0",
    database: str = "refseq_mrna",
    organism: str = "Homo sapiens (taxid:9606)",
    status_callback: Optional[Callable[[str], None]] = None
) -> str:
    """
    Executa a pesquisa de forma síncrona/blocking e extrai os resultados textuais
    do HTML de resposta do NCBI.
    """
    def log(msg):
        if status_callback:
            status_callback(msg)
            
    params = _build_params(forward_primer, reverse_primer, exon_junction_span, database, organism)
    sess = requests.Session()
    
    log("Enviando requisição ao NCBI...")
    try:
        resp = sess.post(PRIMER_BLAST_URL, data=params, timeout=15)
        resp.raise_for_status()
    except Exception as e:
        return f"Erro ao acessar NCBI: {e}"

    match_job = re.search(r'name="job_key"\s+value="([^"]+)"', resp.text, re.IGNORECASE)
    match_ctg = re.search(r'name="ctg_time"\s+value="([^"]+)"', resp.text, re.IGNORECASE)
    
    if not match_job:
        # Tenta RID caso venha direto
        match = re.search(r'RID=([A-Z0-9]+)', resp.text)
        if not match:
            # Pode ser que os resultados já estejam prontos no HTML atual (ex: nenhum erro)
            if "Primer-BLAST can only find non-specific primers" in resp.text or "alignInfo" in resp.text:
                return _parse_html_results(resp.text)
            return "Erro: 'job_key' ou 'RID' não encontrado na resposta do NCBI. O formato da página mudou."

    job_key = match_job.group(1)
    ctg_time = match_ctg.group(1) if match_ctg else ""
    
    poll_params = {
        "CMD": "request",
        "job_key": job_key,
        "ctg_time": ctg_time,
        "CheckStatus": "Check"
    }
    
    attempts = 0
    while attempts < 30: # Limite de 150 segundos
        attempts += 1
        log(f"Aguardando NCBI processar ({attempts * 5}s)...")
        time.sleep(5)
        try:
            poll_resp = sess.get(PRIMER_BLAST_URL, params=poll_params, timeout=15)
            poll_resp.raise_for_status()
        except Exception as e:
            return f"Erro durante a verificação (polling): {e}"
            
        if "alignInfo" in poll_resp.text or 'class="error"' in poll_resp.text or "pr_Result" in poll_resp.text:
            log("Processamento concluído. Extraindo dados...")
            return _parse_html_results(poll_resp.text)
            
        if "Status" in poll_resp.text or "job_key" in poll_resp.text:
            continue
            
        return "Erro desconhecido: O NCBI retornou uma página inesperada."
        
    return "Timeout: A pesquisa demorou muito para responder."


def _parse_html_results(html: str) -> str:
    """Extrai os resultados das tabelas do HTML"""
    soup = BeautifulSoup(html, "html.parser")
    
    align_info = soup.find(id="alignInfo")
    if not align_info:
        # Se não há alignInfo, tenta buscar div pr_Result inteira
        align_info = soup.find("div", class_="pr_Result")
        
    if not align_info:
        # Se houve erro tipo "No primers found"
        err = soup.find("p", class_="error")
        if err:
            return f"Erro retornado pelo NCBI:\n{err.get_text().strip()}"
        return "Nenhum resultado de alinhamento encontrado (div 'alignInfo' ausente)."
        
    # Limpar tags script/style e links, mas preservar table e form
    for tag in align_info(["script", "style", "a", "button", "i"]):
        tag.decompose()
        
    lines = align_info.get_text(separator="\n").split("\n")
    cleaned_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped and stripped not in ["Help", "Download primer pairs", "Text", "CSV", "Tabular"]:
            cleaned_lines.append(stripped)
            
    if not cleaned_lines:
        return "Resultados processados, mas vazios."
        
    return "\n".join(cleaned_lines)


def _build_params(forward_primer, reverse_primer, exon_junction_span, database, organism):
    return {
        "CMD": "request",
        "PRIMER_LEFT_INPUT": forward_primer,
        "PRIMER_RIGHT_INPUT": reverse_primer,
        "PRIMER_ON_SPLICE_SITE": exon_junction_span,
        "PRIMER_SPECIFICITY_DATABASE": database,
        "ORGANISM": organism,
        "SEARCH_SPECIFIC_PRIMER": "on",
        "SPLICE_SITE_OVERLAP_5END": "7",
        "SPLICE_SITE_OVERLAP_3END": "4",
        "SPLICE_SITE_OVERLAP_3END_MAX": "8",
        "MIN_INTRON_SIZE": "1000",
        "MAX_INTRON_SIZE": "1000000",
        "PRIMER_PRODUCT_MIN": "70",
        "PRIMER_PRODUCT_MAX": "1000",
        "PRIMER_NUM_RETURN": "10",
        "PRIMER_MIN_TM": "57.0",
        "PRIMER_OPT_TM": "60.0",
        "PRIMER_MAX_TM": "63.0",
        "PRIMER_MAX_DIFF_TM": "3"
    }
