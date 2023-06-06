from typer import Typer
from useful_rdkit_utils import mol2numpy_fp
from chembl_structure_pipeline import standardize_mol
from chembl_structure_pipeline import get_parent_mol
from rdkit import Chem
import numpy as np
from typing import Optional
from rich.progress import track
from rdkit.Chem import AllChem
app = Typer(help="RDKit CLI tools")

@app.command()
def standardize(file_path: str,
                remove_salt: bool=True,
                out_file: Optional[str]=None):
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
def calcfp(file_path: str,
           radius: int=2,
           nBits: int=2048,
           out_file: Optional[str]=None):
    mols = [m for m in Chem.ForwardSDMolSupplier(file_path)]
    fps = []
    for m in track(mols):
        fps.append(mol2numpy_fp(m, radius=radius, nBits=nBits))
    if out_file == None:
        out_file = f"{file_path}".replace('.sdf','.npy')
    with open(out_file, 'wb') as outf:
        np.save(outf, np.array(fps))


def main():
    app()

"""
if __name__=="__main__":
    app()
"""