# Validação da função de primers em junções éxon–éxon

Executado em 2 de agosto de 2026:

- compilação de todos os módulos Python: aprovada;
- importação do aplicativo: aprovada;
- teste de região 100% conservada: aprovado;
- teste de gap interrompendo região conservada: aprovado;
- teste de geração de pares de primers: aprovado;
- leitura de intervalos e números de éxons do GenBank XML: aprovada;
- preservação de coordenadas na presença de bases IUPAC ambíguas: aprovada;
- mapeamento do transcrito de referência através de gaps do alinhamento: aprovado;
- primers forward e reverse atravessando junções, com ancoragens 5′/3′: aprovados;
- filtro por junção antes do limite de ranking: aprovado;
- rejeição de intervalos inválidos, sobrepostos, não contíguos ou fora do transcrito: aprovada;
- regressão do modo convencional de primers: aprovada;
- exportação dos metadados de éxons e junções em CSV/XLSX: aprovada;
- reinício local usa o mesmo interpretador e reabre o código atualizado: aprovado;
- reinício de uma futura versão empacotada: aprovado;
- versão semântica centralizada e exibida na interface: aprovada;
- rodapé e botão de atualização visíveis em janela de 1260 × 820: aprovado;
- suíte automatizada: **26 testes aprovados**;
- validação sintática dos scripts do macOS com `bash -n`: aprovada;
- permissões de execução dos arquivos `.command`: configuradas no ZIP.

## Compatibilidade macOS

O código usa Python, Tkinter e bibliotecas multiplataforma. Foram adicionados scripts próprios para a estrutura de ambientes virtuais do macOS (`.venv/bin/python`).

O pacote não foi executado fisicamente em um computador Mac durante a geração. Por isso, a validação confirma a portabilidade do código e a sintaxe dos scripts, mas não substitui um teste final em macOS real.

## Não executado com dados reais

As chamadas externas não foram executadas neste ambiente de construção:

- NCBI ESearch/EFetch;
- EMBL-EBI Clustal Omega Job Dispatcher;
- IDT SciTools Plus/OligoAnalyzer.

Motivos: o ambiente usado para empacotar o projeto não tinha as credenciais da IDT e não foi usado para submeter trabalhos reais aos serviços externos.
