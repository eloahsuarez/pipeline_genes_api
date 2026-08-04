import webbrowser
import urllib.parse
import requests

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
    na página de resultados.
    """
    
    # Parâmetros padrão necessários para submeter um formulário de Primer-BLAST
    params = {
        "CMD": "request",
        "PRIMER_LEFT_INPUT": forward_primer,
        "PRIMER_RIGHT_INPUT": reverse_primer,
        "PRIMER_ON_SPLICE_SITE": exon_junction_span,
        "PRIMER_SPECIFICITY_DATABASE": database,
        "ORGANISM": organism,
        "SEARCH_SPECIFIC_PRIMER": "on",
        # Parâmetros de compatibilidade
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

    import re
    # NCBI Primer-BLAST responde a requests POST com uma página de "Format Request"
    # que contém um meta refresh ou um Job ID (RID).
    try:
        response = requests.post(PRIMER_BLAST_URL, data=params, timeout=15)
        response.raise_for_status()
        
        # Procura pelo RID na resposta
        # Exemplo: name="RID" value="XYZ123" ou URL com RID=XYZ123
        match = re.search(r'RID=([A-Z0-9]+)', response.text)
        if not match:
            match = re.search(r'name="RID"\s+value="([A-Z0-9]+)"', response.text, re.IGNORECASE)
            
        if match:
            rid = match.group(1)
            result_url = f"{PRIMER_BLAST_URL}?RID={rid}"
            webbrowser.open(result_url)
        else:
            # Se falhar em achar o RID, cai de volta para abrir o formulário preenchido
            query_string = urllib.parse.urlencode(params)
            full_url = f"{PRIMER_BLAST_URL}?{query_string}"
            webbrowser.open(full_url)
    except Exception as e:
        print(f"Erro ao submeter para o Primer-BLAST: {e}")
        # Fallback para GET preenchido
        query_string = urllib.parse.urlencode(params)
        full_url = f"{PRIMER_BLAST_URL}?{query_string}"
        webbrowser.open(full_url)
