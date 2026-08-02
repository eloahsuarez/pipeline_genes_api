# Gene Conservado — NCBI + Clustal Omega + IDT

Versão atual: **1.1.0**

Aplicativo de desktop compatível com **macOS e Windows** que automatiza este fluxo:

1. pesquisa registros nucleotídicos no NCBI/GenBank;
2. filtra os registros por organismo, gene, tipo, tamanho e qualidade;
3. permite marcar manualmente quais sequências devem ser usadas;
4. envia o FASTA ao Clustal Omega do EMBL-EBI;
5. encontra regiões conservadas pelo limiar de identidade e cobertura;
6. gera pares candidatos de primers, inclusive com um primer atravessando uma junção éxon–éxon;
7. opcionalmente consulta OligoAnalyzer, hairpin, self-dimer e heterodimer na API da IDT;
8. exporta FASTA, CSV, JSON e Excel.

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

O arquivo de configuração **não grava a senha nem o Client Secret**. Não envie credenciais pelo ChatGPT.

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

### Primers em junções éxon–éxon

Na aba **Conservação e primers**, marque **Exigir que pelo menos um primer do par atravesse uma junção**. Em seguida:

1. escolha um transcrito de referência entre os acessos selecionados que possuem anotação de éxons;
2. ajuste a quantidade mínima de bases que deve permanecer nos lados 5′ e 3′ do primer;
3. gere os pares normalmente.

O programa lê as coordenadas das features `exon` do registro GenBank, mapeia o transcrito de referência para o alinhamento e mantém apenas pares em que ao menos um oligo contém bases dos dois éxons. Os padrões de 7 bases no lado 5′ e 4 bases no lado 3′ evitam aceitar uma sobreposição de apenas uma base.

Esse modo requer um registro de mRNA/RNA processado com pelo menos dois éxons contíguos anotados. A disponibilidade dessas features varia entre organismos e bancos; os limites não são inferidos a partir de gaps do Clustal. A isoforma de referência deve ser escolhida conscientemente, pois isoformas diferentes podem ter junções diferentes.

## Limitações

- não executa Primer-BLAST automaticamente;
- não interpreta clinicamente variantes;
- não escolhe sozinho entre isoformas biologicamente relevantes;
- o modo de junção depende da anotação de éxons do transcrito de referência e não monta cDNA a partir de uma sequência genômica;
- os endpoints da IDT dependem da disponibilidade e das permissões da conta;
- serviços públicos podem impor limites, indisponibilidade temporária ou mudanças de API;
- o pacote foi validado em ambiente Linux, mas não foi executado fisicamente em um Mac durante a geração do ZIP.

## Arquivos principais

- `app.py`: interface gráfica e coordenação;
- `ncbi_client.py`: ESearch e EFetch;
- `ebi_client.py`: Job Dispatcher/Clustal Omega;
- `idt_client.py`: OAuth e OligoAnalyzer;
- `bioinformatics.py`: conservação, consenso e primers;
- `exporter.py`: exportação dos resultados;
- `INICIAR_NO_MAC.command`: instalação automática e inicialização no macOS.
