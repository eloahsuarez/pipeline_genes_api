# Validação da aplicação

Executado em 2 de agosto de 2026:

- compilação de todos os módulos Python: aprovada;
- importação do aplicativo: aprovada;
- teste de região 100% conservada: aprovado;
- teste de gap interrompendo região conservada: aprovado;
- teste de geração de pares de primers: aprovado;
- alinhamento global pareado de duas sequências, preservando cabeçalhos e bases: aprovado;
- roteamento de três ou mais sequências para o cliente Clustal Omega: aprovado;
- desenho direto a partir de uma sequência sem Clustal: aprovado;
- bloqueio de primers sobre ou através de bases `N`, sem deslocar coordenadas: aprovado;
- leitura de intervalos e números de éxons do GenBank XML: aprovada;
- preservação de coordenadas na presença de bases IUPAC ambíguas: aprovada;
- mapeamento do transcrito de referência através de gaps do alinhamento: aprovado;
- primers forward e reverse atravessando junções, com ancoragens 5′/3′: aprovados;
- filtro por junção antes do limite de ranking: aprovado;
- rejeição de intervalos inválidos, sobrepostos, não contíguos ou fora do transcrito: aprovada;
- regressão do modo convencional de primers: aprovada;
- exportação dos metadados de éxons e junções em CSV/XLSX: aprovada;
- descoberta segura de BLAST+ e MFEprimer por caminho local ou `PATH`: aprovada;
- cache privado e transacional dos índices, sem alterar o FASTA original: aprovado;
- parsing de hits BLAST e produtos MFEprimer 4.5 em JSON/TSV: aprovado;
- exportação do relatório de especificidade em JSON e CSV: aprovada;
- geração automática do FASTA F/R dos melhores pares, com validação de ranks e bases: aprovada;
- submissão BLAST simulada com `refseq_mrna` e `refseq_select_rna`, taxid 9606 e perfil
  sensível para oligos (`blastn`, sem megablast, `WORD_SIZE=7`, `EXPECT=30000`,
  `HITLIST_SIZE=50000`): aprovada;
- validação fail-closed dos identificadores de banco, sem fallback silencioso para um
  banco diferente: aprovada;
- parsing XML2_S com namespace, múltiplas queries, taxonomia, sequências alinhadas e
  HSPs parciais: aprovado;
- rejeição fail-closed de XML incompleto, consulta ausente, erro embutido, banco vazio e
  ausência de `db-num`/`db-len` positivos em qualquer consulta: aprovada;
- correlação F-R, F-F e R-R no mesmo acesso, em fitas opostas, 3′ voltados para dentro e
  limite inclusivo de amplicon: aprovada;
- extrapolação conservadora das extremidades ausentes de HSPs parciais, respeitando a
  fita: aprovada;
- estimativa de diferenças no oligo inteiro e na janela 3′, incluindo bases fora do HSP,
  substituições e gaps: aprovada;
- exclusão de produtos quando um sítio excede o máximo de diferenças estimadas: aprovada;
- classificação de produtos no gene alvo e em outros genes: aprovada;
- prioridade do símbolo de gene inferido sobre menções textuais ao alvo, evitando falso
  parecer de especificidade: aprovada;
- exportação de combinação F-R/F-F/R-R e métricas dos sítios físicos esquerdo/direito no
  CSV remoto, preservando as colunas F/R legadas: aprovada;
- polling com RID/RTOE, XML2_S, intervalo seguro, timeout e erros HTTP simulados: aprovado;
- smoke real do NCBI BLAST com os dois identificadores de banco: aprovado; RIDs retornados
  `6ZKHZN34014` e `6ZKHYN1E014`;
- resposta real do perfil sensível em `refseq_mrna`, com 10.000 hits por primer:
  aprovada; o RID `6ZNRAXA8016` recuperou os três transcritos CCR5 de 99 bp e os
  off-targets ATP11B e CNTNAP3B, este último como R-R de 2.969 bp;
- smoke real com dois pares, BLAST+ 2.17 e MFEprimer 4.5: aprovado;
- reinício local usa o mesmo interpretador e reabre o código atualizado: aprovado;
- reinício de uma futura versão empacotada: aprovado;
- versão semântica centralizada e exibida na interface: aprovada;
- rodapé e botão de atualização visíveis em janela de 1260 × 820: aprovado;
- configuração e resultados de especificidade separados nas duas modalidades: aprovado;
- tabela e detalhes de especificidade visíveis em 1260 × 820 (198 px/106 px) e
  1080 × 700 (118 px/66 px): aprovados;
- validação amigável de campos numéricos remotos e limite de 1.000 produtos renderizados,
  sem truncar o relatório exportável: aprovados;
- rolagem das configurações de especificidade no tamanho mínimo: aprovada;
- pares candidatos em aba própria, com área útil de 427 px e rolagem horizontal/vertical: aprovado;
- mapa dos melhores pares em subaba própria, redimensionável e com rolagem vertical: aprovado;
- coordenadas inclusivas 1-based convertidas sem erro de uma base: aprovadas;
- orientação Forward para a direita e Reverse para a esquerda: aprovada;
- régua do consenso e do transcrito de referência, incluindo extremos: aprovada;
- abertura real do mapa no macOS, com canvas de 1093 × 414 px e 32 elementos: aprovada;
- separação entre configurações comuns e credenciais: aprovada;
- credenciais ausentes do JSON local e dos arquivos exportados: aprovada;
- permissões privadas do arquivo local no macOS/Linux: aprovadas;
- cofre nativo detectado nesta máquina: `keyring.backends.macOS.Keyring`;
- exclusão dos dados locais e da entrada no cofre: aprovada;
- suíte automatizada: aprovada; os smokes externos opcionais permanecem separados da
  execução padrão;
- validação sintática dos scripts do macOS com `bash -n`: aprovada;
- permissões de execução dos arquivos `.command`: configuradas no ZIP.

## Compatibilidade macOS

O código usa Python, Tkinter e bibliotecas multiplataforma. Foram adicionados scripts próprios para a estrutura de ambientes virtuais do macOS (`.venv/bin/python`).

Nesta atualização, a interface Tkinter foi aberta fisicamente no macOS e o mapa foi
renderizado com um par de validação sobre um alvo de 3.358 nt. A validação em Windows 11
continua dependendo de um teste manual nesse sistema.

## Serviços de rede neste ciclo

O NCBI BLAST remoto foi executado com dados reais em `refseq_mrna` e
`refseq_select_rna`, retornando os RIDs registrados acima. Além dos smokes de banco, o
perfil sensível da versão 1.8.0 foi submetido com os primers CCR5 do caso de regressão e
limite de 10.000 hits por primer. Essa resposta real confirmou a correlação F-R dos três
transcritos alvo e a correlação R-R do produto CNTNAP3B de 2.969 bp. O limite padrão de
50.000 hits, bem como os casos CLK3 e NKD1 mais profundos na lista, permanecem cobertos
por testes determinísticos; não foram submetidos em outro job externo neste ciclo. As
integrações locais BLAST+/MFEprimer também foram executadas.

As chamadas externas abaixo não foram submetidas com dados reais neste ciclo:

- NCBI ESearch/EFetch;
- EMBL-EBI Clustal Omega Job Dispatcher;
- IDT SciTools Plus/OligoAnalyzer.

O ambiente usado para empacotar o projeto não tinha as credenciais da IDT. Os smokes do
BLAST confirmam a aceitação dos dois identificadores pelo serviço. O aplicativo não
executa o alinhamento global Needleman–Wunsch usado pelo Primer-BLAST; portanto, a
extrapolação conservadora de HSPs e as diferenças estimadas não constituem uma reprodução
integral do serviço oficial.
