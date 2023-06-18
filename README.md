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

╭─ Options ──────────────────────────────────────────────────────────────╮
│ --install-completion          Install completion for the current       │
│                               shell.                                   │
│ --show-completion             Show completion for the current shell,   │
│                               to copy it or customize the              │
│                               installation.                            │
│ --help                        Show this message and exit.              │
╰────────────────────────────────────────────────────────────────────────╯
╭─ Commands ─────────────────────────────────────────────────────────────╮
│ calcfp            Calculate Morgan fingerprint                         │
│ checkprops        Read SDF and return prop names                       │
│ confgen           Generate 3D conformers                               │
│ freewilson        CLI based Free Wilson Analysis                       │
│ nibrssfilter      Substructure filters for hit triaginge               │
│ sdf2smi           Converter of SDF to SMILES                           │
│ standardize       Standardize molecules which came from SDF            │
╰────────────────────────────────────────────────────────────────────────╯
```

- standrdize molecules
```
rdkit_cli standardize --help
                                                                                       
 Usage: rdkit_cli standardize [OPTIONS] FILE_PATH                                      

Usage: rdkit_cli standardize [OPTIONS] FILE_PATH

 Standardize molecules which came from SDF
 Molecule standardization command with chemblstructure pipeline

╭─ Arguments ────────────────────────────────────────────────────────────╮
│ *    file_path      TEXT  str path to input SDF [default: None]        │
│                           [required]                                   │
╰────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────╮
│ --remove-salt    --no-remove-salt          remove salt if the compound │
│                                            has salt                    │
│                                            [default: remove-salt]      │
│ --out-file                           TEXT  path of output file, output │
│                                            file is saved in same       │
│                                            driectry of inputfile       │
│                                            [default: None]             │
│ --help                                     Show this message and exit. │
╰────────────────────────────────────────────────────────────────────────╯                                                                                       
```

- calc fp and save results as npy

```
$ rdkit_cli calcfp --help

 Usage: rdkit_cli calcfp [OPTIONS] FILE_PATH

 Calculate Morgan fingerprint
 Calculate Morgan fingerprint from given SDF, and save the FP as numpy
 array.

╭─ Arguments ────────────────────────────────────────────────────────────╮
│ *    file_path      TEXT  str path to input SDF [default: None]        │
│                           [required]                                   │
╰────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────╮
│ --radius          INTEGER  radius of the fingerprint [default: 2]      │
│ --nbits           INTEGER  number of bits of MorganFP [default: 2048]  │
│ --out-file        TEXT     path of output file, output file is saved   │
│                            in same driectry of inputfile and extension │
│                            of the file should be npy                   │
│                            [default: None]                             │
│ --help                     Show this message and exit.                 │
╰────────────────────────────────────────────────────────────────────────╯                                                                                       
```

- generate conformer from 2D sdf
```
$ rdkit_cli confgen --help

 Usage: rdkit_cli confgen [OPTIONS] FILE_PATH [NUMCONFS]

 Generate 3D conformers
 Generate 3D conformers from SDF and save them in new SDF. Compounds
 should be standardize and de-salted.

╭─ Arguments ────────────────────────────────────────────────────────────╮
│ *    file_path      TEXT        str path to input SDF [default: None]  │
│                                 [required]                             │
│      numConfs       [NUMCONFS]  the number of conformers to generate   │
│                                 [default: 10]                          │
╰────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────╮
│ --randomseed        INTEGER  provide a seed for the random number      │
│                              generator so that the same coordinates    │
│                              can be obtained for a molecule on         │
│                              multiple runs. If -1, the RNG will not be │
│                              seeded.                                   │
│                              [default: None]                           │
│ --out-file          TEXT     path of output file, output file is saved │
│                              in same driectry of inputfile and         │
│                              extension of the file should be npy       │
│                              [default: None]                           │
│ --help                       Show this message and exit.               │
╰────────────────────────────────────────────────────────────────────────╯
```


- chek props in SDF
```
$ rdkit_cli checkprops --help

 Usage: rdkit_cli checkprops [OPTIONS] FILE_PATH

 Read SDF and return prop names

╭─ Arguments ────────────────────────────────────────────────────────────╮
│ *    file_path      TEXT  str path to input SDF [default: None]        │
│                           [required]                                   │
╰────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────╮
│ --help          Show this message and exit.                            │
╰────────────────────────────────────────────────────────────────────────╯
```

- sdf2smi
```
$ rdkit_cli sdf2smi --help

 Usage: rdkit_cli sdf2smi [OPTIONS] FILE_PATH

 Converter of SDF to SMILES

╭─ Arguments ────────────────────────────────────────────────────────────╮
│ *    file_path      TEXT  str path to input SDF [default: None]        │
│                           [required]                                   │
╰────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────╮
│ --separator        TEXT  separator of output smi file                  │
│                          commma|tab|whitespace                         │
│                          [default: comma]                              │
│ --out-file         TEXT  path of output file, output file is saved in  │
│                          same driectry of inputfile                    │
│                          [default: None]                               │
│ --help                   Show this message and exit.                   │
╰────────────────────────────────────────────────────────────────────────
```

- NIBRSS Filter
```
 rdkit_cli nibrssfilter --help

 Usage: rdkit_cli nibrssfilter [OPTIONS] FILE_PATH SMILESCOL

 Substructure filters for hit triaginge
 More details on this work can be found in a recent publication:
 Schuffenhauer, A. et al. Evolution of Novartis' small molecule screening
 deck design, J. Med. Chem. (2020), DOI.
 https://dx.doi.org/10.1021/acs.jmedchem.0c01332

╭─ Arguments ────────────────────────────────────────────────────────────╮
│ *    file_path      TEXT  str path to input csv the files should have  │
│                           header line                                  │
│                           [default: None]                              │
│                           [required]                                   │
│ *    smilescol      TEXT  name of SMILES column [default: None]        │
│                           [required]                                   │
╰────────────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────────────╮
│ --separator        TEXT  separator of csv, comma|tab|white-space       │
│                          [default: comma]                              │
│ --out-file         TEXT  path of output file csv [default: None]       │
│ --help                   Show this message and exit.                   │
╰────────────────────────────────────────────────────────────────────────╯
```
