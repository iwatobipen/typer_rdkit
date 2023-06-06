# typer_rdkit

## Description
 This package is CLI tool of rdkit

## Requirements
- rdkit
- numpy
- typer
- typer-cli
- useful_rdkit_utils

## Install
```
gh repo clone iwatobipen/typer_rdkit
cd typer_rdkit
pip install -e .
```

## Basic Usage
- if you would like to use CLI completion
```
rdkit_cli --install-completion
```

```
 iwatobipen@penguin:~/dev/sandbox/typertest$ rdkit_cli --help
                                                                                                                                      
 Usage: rdkit_cli [OPTIONS] COMMAND [ARGS]...                                                                                         
                                                                                                                                      
 RDKit CLI tools                                                                                                                      
                                                                                                                                      
╭─ Options ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --install-completion          Install completion for the current shell.                                                            │
│ --show-completion             Show completion for the current shell, to copy it or customize the installation.                     │
│ --help                        Show this message and exit.                                                                          │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ calcfp                                                                                                                             │
│ standardize                                                                                                                        │
╰────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

```

- standrdize molecules
```
rdkit_cli standardize --help
                                                                                       
 Usage: rdkit_cli standardize [OPTIONS] FILE_PATH                                      
                                                                                       
╭─ Arguments ─────────────────────────────────────────────────────────────────────────╮
│ *    file_path      TEXT  [default: None] [required]                                │
╰─────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ───────────────────────────────────────────────────────────────────────────╮
│ --remove-salt    --no-remove-salt          [default: remove-salt]                   │
│ --out-file                           TEXT  [default: None]                          │
│ --help                                     Show this message and exit.              │
╰─────────────────────────────────────────────────────────────────────────────────────╯
```

- calc fp and save results as npy

```
$ rdkit_cli calcfp --help
                                                                                       
 Usage: rdkit_cli calcfp [OPTIONS] FILE_PATH                                           
                                                                                       
╭─ Arguments ─────────────────────────────────────────────────────────────────────────╮
│ *    file_path      TEXT  [default: None] [required]                                │
╰─────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ───────────────────────────────────────────────────────────────────────────╮
│ --radius          INTEGER  [default: 2]                                             │
│ --nbits           INTEGER  [default: 2048]                                          │
│ --out-file        TEXT     [default: None]                                          │
│ --help                     Show this message and exit.                              │
╰─────────────────────────────────────────────────────────────────────────────────────╯
```