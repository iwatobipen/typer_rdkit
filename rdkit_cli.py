import typer
from typer import Typer
from useful_rdkit_utils import mol2numpy_fp
from chembl_structure_pipeline import standardize_mol
from chembl_structure_pipeline import get_parent_mol
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdForceFieldHelpers

import numpy as np
from typing import Optional
from rich.progress import track
from rdkit.Chem import AllChem
app = Typer(help="RDKit CLI tools")

@app.command()
def standardize(file_path: str=typer.Argument(..., help="str path to input SDF"),
                remove_salt: bool=typer.Argument(True, help="remove salt if the compound has salt"),
                out_file: Optional[str]=typer.Argument(None, help="path of output file, output file is saved in same driectry of inputfile")):
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
           radius: int=typer.Argument(2, help="radius of the fingerprint"),
           nBits: int=typer.Argument(2048, help="number of bits of MorganFP"),
           out_file: Optional[str]=typer.Argument(None, help="path of output file, output file is saved in same driectry of inputfile and extension of the file should be npy")):
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
            randomSeed: Optional[int]=typer.Argument(None, help="provide a seed for the random number generator so that the same coordinates can be obtained for a molecule on multiple runs. If -1, the RNG will not be seeded."),
            out_file: Optional[str]=typer.Argument(None, help="path of output file, output file is saved in same driectry of inputfile and extension of the file should be npy")):
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

   
def main():
    app()
