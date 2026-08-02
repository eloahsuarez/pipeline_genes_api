# Memória do projeto — Gene Conservado

Estas regras devem ser preservadas em todas as alterações futuras deste repositório.

## Interface

- O aplicativo continua sendo uma interface desktop em Python/Tkinter.
- Exibir o código da versão no canto superior esquerdo e no título da janela.
- Manter o botão **Atualizar ferramenta** sempre visível na barra inferior, à esquerda
  de **Exportar projeto**.
- Empacotar o cabeçalho e o rodapé antes do `Notebook` expansível para evitar que a
  barra inferior seja empurrada para fora da janela em telas menores.

## Versionamento

- `version.py` é a fonte única da versão, por meio de `APP_VERSION`.
- Usar versionamento semântico `MAJOR.MINOR.PATCH`.
- Incrementar `PATCH` em correções, `MINOR` em funcionalidades compatíveis e `MAJOR`
  em alterações incompatíveis.
- Toda atualização visível ao usuário deve atualizar `APP_VERSION` e a versão indicada
  no `README.md`.

## Atualização local

- O botão **Atualizar ferramenta** não baixa arquivos da internet.
- Ele deve encerrar a instância atual e reabrir o aplicativo com o mesmo interpretador,
  carregando os arquivos já atualizados na pasta do projeto.
- Não reiniciar enquanto uma consulta ou análise estiver em andamento.
- Preservar `.venv`, configurações, credenciais e arquivos exportados.

## Validação e publicação

- Antes de concluir uma alteração, executar `python -m pytest -q` e verificar a
  compilação dos módulos modificados.
- Não versionar `.venv`, caches, segredos ou configurações locais; respeitar o
  `.gitignore`.
- O repositório remoto oficial é `https://github.com/eloahsuarez/pipeline_genes_api`.
