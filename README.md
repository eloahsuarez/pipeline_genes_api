# Gene Conservado — NCBI + Clustal Omega + BLAST+/MFEprimer + IDT

Versão atual: **1.10.2**

Aplicativo de desktop compatível com **macOS e Windows** que automatiza este fluxo:

1. pesquisa registros nucleotídicos no NCBI/GenBank;
2. filtra os registros por organismo, gene, tipo, tamanho e qualidade;
3. permite marcar manualmente quais sequências devem ser usadas;
4. alinha localmente quando há duas sequências ou envia três ou mais ao Clustal Omega do EMBL-EBI;
5. encontra regiões conservadas pelo limiar de identidade e cobertura;
6. também permite desenhar diretamente a partir de uma única sequência, sem alinhamento;
7. gera pares candidatos de primers, inclusive com um primer atravessando uma junção éxon–éxon;
8. envia o FASTA dos melhores primers ao NCBI e prevê produtos F-R, F-F e R-R;
9. opcionalmente verifica os mesmos pares em um banco local com BLAST+ e MFEprimer;
10. mostra a localização dos melhores pares ao longo da sequência alvo;
11. opcionalmente consulta OligoAnalyzer, hairpin, self-dimer e heterodimer na API da IDT;
12. exporta FASTA, CSV, JSON e Excel.

## Instalação no macOS

### Forma mais simples

1. Extraia completamente o ZIP.
2. Abra a pasta extraída.
3. Clique duas vezes em `INICIAR_NO_MAC.command`.

Na primeira execução, esse arquivo cria o ambiente virtual, instala as dependências e abre o programa. Nas execuções seguintes, ele apenas inicia o aplicativo.

### Requisitos no Mac

- macOS com processador Intel ou Apple Silicon;
- Python 3.11 ou superior;
- Tkinter disponível.

É recomendável usar o instalador oficial do Python, pois ele normalmente inclui o suporte ao Tkinter necessário para a interface gráfica.

Caso o macOS bloqueie o arquivo por ele ter sido baixado da internet, clique com o botão direito em `INICIAR_NO_MAC.command` e selecione **Abrir**. Alternativamente, no Terminal aberto na pasta do programa, execute:

```bash
chmod +x *.command
./INICIAR_NO_MAC.command
```

Também estão disponíveis os scripts separados:

- `instalar_mac.command`: recria o ambiente e instala as dependências;
- `executar_mac.command`: abre o aplicativo depois da instalação.

## Aplicar atualizações locais

Quando os arquivos do programa forem atualizados diretamente nesta pasta, clique em
**Atualizar ferramenta**, na barra inferior. O aplicativo encerra a instância atual e
abre novamente com o mesmo Python, carregando o código novo.

O código da versão aparece no canto superior esquerdo da janela e também no título.
Cada nova atualização deve alterar `APP_VERSION` no arquivo `version.py`.

O botão não baixa arquivos da internet e não altera a `.venv`: ele apenas reinicia o
programa para aplicar o que já foi substituído na pasta do projeto. Aguarde qualquer
consulta ou análise em andamento terminar. Resultados que ainda estiverem apenas na
memória devem ser exportados antes do reinício.

## Configurações e credenciais locais

O aplicativo carrega os dados locais ao abrir e salva automaticamente ao fechar ou
reiniciar. Também é possível usar o menu **Dados locais** para salvar imediatamente ou
apagar tudo que ficou armazenado neste computador.

- configurações comuns ficam fora da pasta do projeto, no diretório privado de dados do
  usuário;
- credenciais ficam no cofre seguro do sistema: **Chaves do macOS** ou **Gerenciador de
  Credenciais do Windows**;
- API key do NCBI, Client Secret e senha da IDT não são incluídos em configurações
  exportadas, logs, Git ou GitHub.

## Instalação no Windows 11

1. Extraia a pasta.
2. Execute `instalar.bat` uma vez.
3. Execute `executar.bat` para abrir o programa.

Requer Python 3.11 ou superior instalado com o comando `py` disponível.

## Credenciais

### NCBI

- O e-mail é necessário.
- A API key é opcional, mas aumenta o limite padrão de requisições.

### EMBL-EBI

- O Job Dispatcher exige um e-mail válido.
- Não exige chave de API.

### IDT

Na conta IDT, abra **My Account → API access** e gere Client ID e Client Secret. O programa solicita um token OAuth quando a análise é iniciada.

Depois da análise, selecione uma linha no ranking e clique em **Ver sequências do par
selecionado** (ou dê duplo clique na linha) para abrir as sequências Forward e Reverse
do par analisado.

Os arquivos de configuração exportados **não gravam API keys, senha nem Client Secret**.
O salvamento automático local usa o cofre seguro do sistema. Não envie credenciais pelo
ChatGPT.

### NCBI remoto (especificidade dos pares)

Na aba **6. Especificidade → NCBI remoto**, escolha o banco e o aplicativo monta
automaticamente um FASTA com os primeiros pares do ranking. Cada par gera duas entradas,
por exemplo `pair_0001_F` e `pair_0001_R`. Use **Visualizar FASTA dos melhores pares**
para conferir a entrada ou **Exportar FASTA** para salvá-la.

Estão disponíveis duas opções:

- **RefSeq mRNA** (`refseq_mrna`), padrão mais abrangente, que inclui múltiplos mRNAs e
  isoformas RefSeq e se aproxima melhor do escopo de transcritos usado pelo Primer-BLAST;
- **RefSeq Select RNA** (`refseq_select_rna`), banco compacto e mais rápido, com um
  transcrito representativo selecionado por gene codificante, mas sem cobertura de todas
  as isoformas.

Nas duas modalidades de especificidade, **Configuração** e **Resultados** ficam em
subabas separadas para manter a tabela e os detalhes completamente visíveis. Ao concluir
uma análise, o aplicativo abre **Resultados** automaticamente. A configuração possui
rolagem vertical quando a janela está no tamanho mínimo.

A submissão usa a
[Common URL API do BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html)
com um perfil sensível fixo para oligos:

- banco **RefSeq mRNA** (`refseq_mrna`) ou **RefSeq Select RNA**
  (`refseq_select_rna`);
- organismo **Homo sapiens** (`taxid:9606`);
- programa `blastn`, sem megablast, com `WORD_SIZE=7` e `EXPECT=30000`;
- ajuste automático para queries curtas desativado, pois os parâmetros sensíveis são
  enviados explicitamente;
- `HITLIST_SIZE=50000`, para não limitar a triagem aos primeiros poucos transcritos;
- por padrão, somente o melhor par do ranking, amplicons de até 4.000 bp e no máximo
  seis diferenças estimadas em cada sítio de ligação.

`WORD_SIZE=7`, `EXPECT=30000` e o limite de 50.000 sequências acompanham os parâmetros
atualmente expostos pelo
[formulário do Primer-BLAST no NCBI](https://www.ncbi.nlm.nih.gov/tools/primer-blast/).
O
[artigo do Primer-BLAST hospedado pelo NCBI](https://pmc.ncbi.nlm.nih.gov/articles/PMC3412702/)
explica por que a busca sensível precisa ser refinada além do alinhamento local. Esses
valores ampliam a recuperação de alvos divergentes, mas podem aumentar bastante o tempo
de fila e de processamento. Todos os oligos escolhidos são enviados em um único job.

O app não trata dois hits isolados como prova de amplificação. Para cada par candidato,
ele testa as combinações **F-R**, **F-F** e **R-R** e só forma um produto quando dois
sítios atingem o mesmo acesso, estão em fitas opostas, têm as extremidades 3′ voltadas
uma para a outra e respeitam o **Amplicon máximo**. Assim, também detecta situações nas
quais um único oligo pode atuar nos dois lados de um produto.

O BLAST produz alinhamentos locais (HSPs), que às vezes cobrem somente parte do oligo.
Nesses casos, o aplicativo extrapola conservadoramente as coordenadas das bases não
alinhadas conforme a fita e conta cada base fora do HSP como diferença. Também estima as
diferenças totais e as diferenças na janela dos cinco nucleotídeos da extremidade 3′.
Essas contagens são estimativas de triagem: um sítio com mais de seis diferenças
estimadas não participa da formação do produto padrão.

O gene informado na aba NCBI é o alvo. Acesso do transcrito selecionado e símbolo no
título do resultado são usados para separar **Alvo** de **Outro gene**. Para cada produto,
a tabela mostra a combinação de primers e os dois sítios físicos. O CSV inclui
`combinacao_primers` e, para `sitio_esquerdo` e `sitio_direito`, orientação F/R, query
cover, identidade, diferenças estimadas total/3′, coordenadas e fita. As colunas legadas
Forward/Reverse foram mantidas para compatibilidade. O JSON conserva ainda E-value, bit
score e os dados completos dos hits.

Um acesso alvo conhecido tem prioridade na classificação. Fora disso, quando o título
permite inferir um símbolo de gene, ele precisa coincidir exatamente com o gene alvo; uma
simples menção ao nome em outra parte do título não transforma um off-target em alvo. A
busca textual é usada somente quando o título não oferece símbolo estruturado.

Para evitar travar o Tkinter em resultados muito grandes, a tabela e o painel de detalhes
renderizam no máximo 1.000 produtos de uma vez e informam quantos foram omitidos da
visualização. O relatório em memória e as exportações JSON/CSV continuam contendo todos
os produtos encontrados.

As chaves antigas `ncbi_specificity_identity_pct` e
`ncbi_specificity_coverage_pct` deixaram de filtrar a busca remota, porque um corte
prévio poderia eliminar HSPs parciais que precisam da extrapolação conservadora.
Configurações locais antigas podem continuar contendo essas chaves, mas elas são inertes
no perfil sensível e podem ser removidas. O exemplo de configuração já usa os novos
padrões e `ncbi_specificity_max_estimated_mismatches=6`.

O NCBI pode enfileirar uma consulta por vários minutos. A ferramenta respeita o polling
do serviço e impede reinício enquanto a análise está ativa. A API key configurada para
E-utilities não é enviada ao BLAST; somente o e-mail de identificação do cliente.

O resultado remoto é aceito somente quando contém todas as consultas submetidas e
estatísticas positivas `db-num` e `db-len` para cada uma. XML incompleto, banco vazio ou
erro embutido interrompem a análise. Se uma consulta atingir o limite configurado de hits
sem revelar off-target, o parecer é **Inconclusivo: limite de hits atingido**, e não
**Específico no banco**.

Mesmo no RefSeq mRNA, o resultado continua limitado à coleção e aos parâmetros
selecionados: não exclui produtos ausentes do banco, pseudogenes ou DNA genômico. Embora
o modo remoto agora verifique F-R/F-F/R-R, ele **não executa o alinhamento global
Needleman–Wunsch** que o
[Primer-BLAST do NCBI](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) usa para completar
os alinhamentos locais do BLAST. A extrapolação e as diferenças estimadas do aplicativo
não fornecem equivalência integral ao Primer-BLAST nem substituem suas regras adicionais
de especificidade. Para uma checagem complementar, siga o
[guia oficial do Primer-BLAST](https://www.ncbi.nlm.nih.gov/guide/howto/design-pcr-primers/)
ou use a subaba local com uma referência adequada.

### BLAST+ e MFEprimer (especificidade local)

Essa etapa é opcional e usa executáveis instalados no computador; eles não são pacotes
Python e não são baixados pelo botão **Atualizar ferramenta**. Instale o
[BLAST+ distribuído pelo NCBI](https://www.ncbi.nlm.nih.gov/books/NBK569861/) e o
[MFEprimer 4.x](https://github.com/quwubin/MFEprimer-3.0/releases) na versão adequada ao
seu sistema. Na aba **6. Especificidade → Banco local — BLAST+ e MFEprimer**, informe
os caminhos completos dos executáveis ou deixe os campos
em branco para procurar em `.tools/bin` e no `PATH`.

Selecione um FASTA de referência representativo — por exemplo, o genoma ou transcriptoma
do organismo — e clique em **Preparar banco**. A ferramenta copia esse FASTA para o
diretório privado de dados do aplicativo e cria ali os índices do BLAST+ e do MFEprimer;
o arquivo original não é modificado. Bancos genômicos grandes podem exigir bastante
tempo e espaço na primeira preparação.

O BLAST+ é usado como triagem de sítios individuais de ligação, enquanto o MFEprimer
avalia os produtos formados pelo par. Um resultado só descreve a especificidade em
relação ao banco FASTA escolhido e não substitui validação experimental.

## O que significa “sequência aplicável ao gene”

Não existe um filtro universal. O programa separa os critérios para que você possa decidir:

- mRNA, CDS, RNA ou sequência genômica;
- RefSeq ou todos os registros;
- incluir ou excluir `PREDICTED` e `partial`;
- exigir que a feature `/gene=` corresponda exatamente ao gene informado;
- tamanho mínimo e máximo;
- seleção manual final de cada registro.

Isso reduz o risco de alinhar misturas biologicamente incompatíveis.

## Regiões conservadas

Com identidade e cobertura em 100%, uma posição só é marcada como conservada quando todas as sequências possuem a mesma base canônica naquela coluna e nenhuma possui gap. Reduzir os limiares permite tolerar variantes ou registros incompletos.

## Primers

A geração local é uma triagem computacional. Ela usa Tm por nearest-neighbor do Biopython, GC, tamanho, diferença de Tm, tamanho do amplicon e penalidades simples. A etapa IDT acrescenta análise físico-química. Ainda é recomendável verificar especificidade contra sequências não alvo, idealmente por Primer-BLAST ou BLAST.

Depois de gerar os primers, a ferramenta abre automaticamente a aba **5. Resultados**,
na subaba **Mapa dos primers**. O mapa exibe os X melhores pares definidos pelo campo
**Top pares**, em ordem de ranking. A régua representa o consenso sem gaps no desenho
convencional ou o transcrito de referência no modo de junção éxon–éxon. Cada linha mostra
o amplicon e as posições dos primers: Forward aponta para a direita e Reverse para a
esquerda. Use a rolagem vertical quando houver muitos pares.

A subaba **Pares candidatos** mantém a tabela detalhada, com barras de rolagem horizontal
e vertical para visualizar todas as colunas e todos os pares.

Quando apenas um registro estiver marcado na aba **Sequências**, use **Desenhar com 1
sequência** ou abra diretamente a aba **Conservação e primers** e clique em **Gerar pares
de primers**. Nesse modo o programa preserva as coordenadas originais e impede que um
primer contenha ou atravesse bases ambíguas `N`; nenhum alinhamento é criado ou exportado.

Com exatamente duas sequências, o aplicativo faz um alinhamento global pareado local.
Com três ou mais, continua usando o Clustal Omega do EMBL-EBI. Essa separação segue a
[orientação do próprio serviço](https://www.ebi.ac.uk/jdispatcher/docs/faqs/clustal/),
pois o Clustal Omega é destinado a alinhamentos múltiplos.

### Primers em junções éxon–éxon

Na aba **Conservação e primers**, marque **Exigir que pelo menos um primer do par atravesse uma junção**. Em seguida:

1. escolha um transcrito de referência entre os acessos selecionados que possuem anotação de éxons;
2. ajuste a quantidade mínima de bases que deve permanecer nos lados 5′ e 3′ do primer;
3. gere os pares normalmente.

O programa lê as coordenadas das features `exon` do registro GenBank, mapeia o transcrito de referência para o alinhamento e mantém apenas pares em que ao menos um oligo contém bases dos dois éxons. Os padrões de 7 bases no lado 5′ e 4 bases no lado 3′ evitam aceitar uma sobreposição de apenas uma base.

Esse modo requer um registro de mRNA/RNA processado com pelo menos dois éxons contíguos anotados. A disponibilidade dessas features varia entre organismos e bancos; os limites não são inferidos a partir de gaps do Clustal. A isoforma de referência deve ser escolhida conscientemente, pois isoformas diferentes podem ter junções diferentes.

## Limitações

- o modo remoto correlaciona F-R/F-F/R-R, mas usa extrapolação conservadora de HSPs e
  diferenças estimadas; sem alinhamento global Needleman–Wunsch, não reproduz todos os
  algoritmos do Primer-BLAST;
- o RefSeq Select RNA não representa todas as isoformas; o RefSeq mRNA é mais abrangente,
  mas nenhum dos dois representa o genoma humano completo;
- o modo local depende da abrangência do FASTA escolhido pelo usuário;
- não interpreta clinicamente variantes;
- não escolhe sozinho entre isoformas biologicamente relevantes;
- o modo de junção depende da anotação de éxons do transcrito de referência e não monta cDNA a partir de uma sequência genômica;
- os endpoints da IDT dependem da disponibilidade e das permissões da conta;
- serviços públicos podem impor limites, indisponibilidade temporária ou mudanças de API;
- a interface foi validada no macOS; ainda é recomendado um teste manual final no Windows 11.

## Arquivos principais

- `app.py`: interface gráfica e coordenação;
- `ncbi_client.py`: ESearch e EFetch;
- `ncbi_blast.py`: BLAST remoto, parsing XML2 e correlação dos hits por par;
- `ebi_client.py`: Job Dispatcher/Clustal Omega;
- `idt_client.py`: OAuth e OligoAnalyzer;
- `bioinformatics.py`: conservação, consenso e primers;
- `specificity.py`: integração local e isolada de BLAST+ e MFEprimer;
- `primer_plot.py`: escala e geometria do mapa de primers;
- `exporter.py`: exportação dos resultados;
- `INICIAR_NO_MAC.command`: instalação automática e inicialização no macOS.
