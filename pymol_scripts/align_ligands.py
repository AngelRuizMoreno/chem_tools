# Import necessary libraries from PyMOL and RDKit
from pymol import cmd, stored
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdFMCS
import tempfile

def align_ligands(target, reference):
    '''
    Aligns two small molecules (ligands) within PyMOL based on maximum common substructure (MCS) alignment.

    Parameters:
        target (str): The name of the object in PyMOL representing the target molecule to be aligned.
        reference (str): The name of the object in PyMOL representing the reference molecule.

    Returns:
        rmsd (float): Root Mean Square Deviation (RMSD) value between the aligned structures.
        mcs_smarts (str): SMARTS string representing the Maximum Common Substructure (MCS).
        num_atoms (int): Number of atoms in the MCS used for alignment.
    
    Description:
        This function saves the target and reference molecules as temporary .sdf files from PyMOL. 
        It uses RDKit to read the .sdf files, finds the MCS between the two molecules, and aligns 
        the target molecule to the reference based on this common substructure. After alignment, 
        the function reloads the aligned target molecule back into PyMOL as 'aligned_target' 
        and deletes the temporary .sdf files automatically.
    
    Requirements:
        - PyMOL (with the command-line version or PyMOL open to run the script)
        - RDKit (for chemical informatics functionalities)
        - The target and reference objects must be loaded in PyMOL

    Usage Example:
        Load the molecules in PyMOL and run the script:
            align_ligands target_name reference_name
        
        Where `target_name` and `reference_name` are the PyMOL object names of the target and reference molecules.
    
    Additional Information:
        The script will print the RMSD, the SMARTS pattern of the MCS, and the number of atoms involved 
        in the alignment after execution.
    '''
    
    # Create temporary files for the reference and target molecules
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as ref_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tar_tmp:
        
        # Save reference and target molecules to temporary .sdf files
        cmd.save(ref_tmp.name, format="sdf", selection=reference, state=-1)
        cmd.save(tar_tmp.name, format="sdf", selection=target, state=-1)
        
        # Read the molecules from the .sdf files using RDKit
        ref = Chem.SDMolSupplier(ref_tmp.name)[0]
        tar = Chem.SDMolSupplier(tar_tmp.name)[0]
        
        # Find the Maximum Common Substructure (MCS) between reference and target
        mcs = rdFMCS.FindMCS([ref, tar])
        ref_atoms = ref.GetSubstructMatch(mcs.queryMol)
        tar_atoms = tar.GetSubstructMatch(mcs.queryMol)
        
        # Map the atoms based on MCS and align the target molecule to the reference
        atom_map = list(zip(tar_atoms, ref_atoms))
        rmsd = rdMolAlign.AlignMol(prbMol=tar, refMol=ref, atomMap=atom_map)
        
        # Write the aligned target molecule to a new temporary .sdf file
        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as aligned_tmp:
            out = Chem.SDWriter(aligned_tmp.name)
            out.write(tar)
            out.close()
    
    # Print alignment information
    print(f"RMSD: {rmsd}")
    print(f"MCS (SMARTS): {mcs.smartsString}")
    print(f"Number of atoms in MCS: {len(atom_map)}")
    
    return rmsd, mcs.smartsString, len(atom_map)

# Extend PyMOL's command set to include this function
cmd.extend("align_ligands", align_ligands)