import os
import sys
import typer
import pathlib
from typer import Typer
from useful_rdkit_utils import mol2numpy_fp
from chembl_structure_pipeline import standardize_mol
from chembl_structure_pipeline import get_parent_mol
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem import AllChem
from rdkit.Chem import RDConfig
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
import numpy as np
import pandas as pd
from typing import Optional
from rich.progress import track
fw_path = os.path.join(RDConfig.RDContribDir, 'FreeWilson')
nibrssf_path = os.path.join(RDConfig.RDContribDir, 'NIBRSubstructureFilters')
sys.path.append(fw_path)
sys.path.append(nibrssf_path)
from freewilson import FWDecompose, FWBuild, predictions_to_csv
import assignSubstructureFilters
app = Typer(help="RDKit CLI tools")

def mol2smi(mol):
    if mol == None:
        return "*"
    else:
        return Chem.MolToSmiles(mol)


@app.command()
def checkprops(file_path: str=typer.Argument(..., help="str path to input SDF")):
    """Read SDF and return prop names
    """
    df = PandasTools.LoadSDF(file_path)
    print(df.columns.to_list())               

@app.command()
def sdf2smi(file_path: str=typer.Argument(..., help="str path to input SDF"),
            separator: str=typer.Option("comma", help="separator of output smi file commma|tab|whitespace"),
            out_file: Optional[str]=typer.Option(None, help="path of output file, output file is saved in same driectry of inputfile")
            ):
    """Converter of SDF to SMILES
    """
    if separator=="comma":
        s=","
    elif separator=="tab":
        s="\t"
    else:
        s=" "
    df = PandasTools.LoadSDF(file_path)
    if out_file == None:
        inputpath = pathlib.Path(file_path)
        from_suffix = inputpath.suffix
        out_file = f"{file_path}".replace(from_suffix,'.smi')
    df['smiles'] = df.ROMol.apply(mol2smi)
    df.to_csv(out_file, sep=s, index=False)


@app.command()
def standardize(file_path: str=typer.Argument(..., help="str path to input SDF"),
                remove_salt: bool=typer.Option(True, help="remove salt if the compound has salt"),
                out_file: Optional[str]=typer.Option(None, help="path of output file, output file is saved in same driectry of inputfile")):
    """Standardize molecules which came from SDF

    Molecule standardization command with chemblstructure pipeline

    """
    mols = [m for m in Chem.ForwardSDMolSupplier(file_path)]
    if out_file == None:
        out_file = f"{file_path}".replace('.sdf','_std.sdf')
    w = Chem.SDWriter(out_file)
    for m in track(mols):
        if m == None:
            pass
        m.RemoveAllConformers()
        AllChem.Compute2DCoords(m)
        stdm = standardize_mol(m)
        if stdm:
            if remove_salt:
                stdm, _ = get_parent_mol(stdm)
            w.write(stdm)

@app.command()
def calcfp(file_path: str=typer.Argument(..., help="str path to input SDF"),
           radius: int=typer.Option(2, help="radius of the fingerprint"),
           nBits: int=typer.Option(2048, help="number of bits of MorganFP"),
           out_file: Optional[str]=typer.Option(None, help="path of output file, output file is saved in same driectry of inputfile and extension of the file should be npy")):
    """Calculate Morgan fingerprint

    Calculate Morgan fingerprint from given SDF, and save the FP as numpy array.

    """
 
    mols = [m for m in Chem.ForwardSDMolSupplier(file_path)]
    fps = []
    for m in track(mols):
        fps.append(mol2numpy_fp(m, radius=radius, nBits=nBits))
    if out_file == None:
        out_file = f"{file_path}".replace('.sdf','.npy')
    with open(out_file, 'wb') as outf:
        np.save(outf, np.array(fps))

@app.command()
def confgen(file_path: str=typer.Argument(..., help="str path to input SDF"),
            numConfs: int=typer.Argument(10, help="the number of conformers to generate"),
            randomSeed: Optional[int]=typer.Option(None, help="provide a seed for the random number generator so that the same coordinates can be obtained for a molecule on multiple runs. If -1, the RNG will not be seeded."),
            out_file: Optional[str]=typer.Option(None, help="path of output file, output file is saved in same driectry of inputfile and extension of the file should be npy")):
    """Generate 3D conformers

    Generate 3D conformers from SDF and save them in new SDF. Compounds should be standardize and de-salted.

    """
 
    mols = [m for m in Chem.ForwardSDMolSupplier(file_path)]
    if out_file == None:
        out_file = f"{file_path}".replace('.sdf','_3d.sdf')
    w = Chem.SDWriter(out_file)
    for m in mols:
        m = Chem.AddHs(m)
        cids = rdDistGeom.EmbedMultipleConfs(m)
        for cid in cids:
            m.SetProp('confId', str(cid))
            w.write(m, confId=cid)
    w.close()

@app.command()
def freewilson(file_path: str=typer.Argument(..., help="str path to input SDF"),
               target_val: str=typer.Argument(..., help="property name of prediction target such as pIC50"),

               usefmcs: bool=typer.Option(False, help="use FMCS for build FW model. If the option is set True, scaffold will be ignored"),
               threshold: float=typer.Option(0.8, help="threshold of FMCS"),
               scaffold: Optional[str]=typer.Option(None, help="path of scaffold shoud be mol format"),
               out_file: Optional[str]=typer.Option(None, help="path of output file")
               ):
    """CLI based Free Wilson Analysis
    
    """
    
    mols = [m for m in Chem.SDMolSupplier(file_path) if m != None]
    scores = [float(m.GetProp(target_val)) for m in mols]
    train_data = {}
    for idx, m in enumerate(mols):
        train_data[Chem.MolToSmiles(m)] = scores[idx]
    if scaffold != None:
        scaf = Chem.MolFromMolFile(scaffold)
    else:
        scaf = None
    if usefmcs or scaf==None:
        print('Will find MCS, it will take long time if you path thousands or more compounds')
        mcs = rdFMCS.FindMCS(mols, threshold=threshold, 
                             atomCompare=rdFMCS.AtomCompare.CompareAny, 
                             completeRingsOnly=True)
        scaf = mcs.queryMol
    decomp = FWDecompose(scaf, mols, scores)
    print(f"Training R^2 is {decomp.r2:0.2f}")
    preds = FWBuild(decomp)
    if out_file == None:
        out_file = f"{file_path}".replace('.sdf','_fw.csv')

    out_file = open(out_file, 'w')
    out_file.write("idx\tsmi\tpred_val\treal_val\n")
    for idx, row in enumerate(train_data.items()):
        out_file.write(f"train_{idx}\t{row[0]}\t*\t{row[1]}\n")
    for idx, pred in enumerate(preds):
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(pred.smiles))
        pred_val = pred.prediction
        if smi in list(train_data.keys()):
            real_val = train_data[smi]
        else:
            real_val = ""
        out_file.write(f"enum_{idx}\t{smi}\t{pred_val}\t{real_val}\n")
    out_file.close()

@app.command()
def nibrssfilter(file_path: str=typer.Argument(..., help="str path to input csv\nthe files should have header line"),
                 smilescol: str=typer.Argument(..., help='name of SMILES column'),
                 separator: str=typer.Option('comma', help="separator of csv, comma|tab|white-space"),
                 out_file: Optional[str]=typer.Option(None, help="path of output file csv")
                 ):
    """Substructure filters for hit triaginge

    More details on this work can be found in a recent publication: Schuffenhauer, A. et al. Evolution of Novartis' small molecule screening deck design, J. Med. Chem. (2020), DOI. https://dx.doi.org/10.1021/acs.jmedchem.0c01332
    """
    if separator=="comma":
        s=","
    elif separator=="tab":
        s="\t"
    else:
        s=" "
    df = pd.read_csv(file_path, sep=s)
    res = assignSubstructureFilters.assignFilters(df, nameSmilesColumn=smilescol)
    tmpdf = pd.DataFrame.from_records(res, columns=assignSubstructureFilters.FilterMatch._fields)
    df = df.merge(tmpdf, how='left', left_index=True, right_index=True)
    if out_file == None:
        inputpath = pathlib.Path(file_path)
        from_suffix = inputpath.suffix
        out_file = f"{file_path}".replace(from_suffix,'_filtered.csv')

    df.to_csv(out_file,sep=s, index=False)
    print('Done')


def main():
    app()
