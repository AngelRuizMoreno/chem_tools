# Import necessary libraries from PyMOL and RDKit
from pymol import cmd, stored
import os, tempfile, subprocess, shutil
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdFMCS

def __init__(self):
    
    self.menuBar.addcascademenu('Plugin', 'ChemTools', 'ChemTools', label='ChemTools')
    self.menuBar.addmenuitem('ChemTools', 'command', 'align_with_rdkit_help', label='Align with RDKit Help', command=lambda: align_with_rdkit_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'align_with_fkcombu_help', label='Align with FKcombu Help', command=lambda: align_with_fkcombu_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'docking_with_smina_help', label='Docking with Smina Help', command=lambda: docking_with_smina_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'pockets_with_fpocket_help', label='Pockets with Fpocket Help', command=lambda: pockets_with_fpocket_Help())


def align_with_rdkit_Help():
    Usage="""

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
    """
    print(Usage)
    return

def align_with_fkcombu_Help():
    Usage="""
    
    Aligns the ligand to the receptor in PyMOL using fkcombu, with the receptor surface as a constraint.

    Parameters:
    - receptor (str): Name of the receptor object in PyMOL.
    - reference (str): Name of the reference ligand object in PyMOL.
    - target (str): Name of the target ligand object in PyMOL.

    Returns:
    - None
    """
    print(Usage)
    return
    
def docking_with_smina_Help():
    Usage="""

    docking_with_smina with specified parameters, loading receptor and ligand directly from PyMOL session.

    Parameters:
    - receptor (str): Name of the receptor object in PyMOL.
    - ligand (str): Name of the ligand object in PyMOL.
    - center_x (float): X-coordinate for the center of the search box.
    - center_y (float): Y-coordinate for the center of the search box.
    - center_z (float): Z-coordinate for the center of the search box.
    - size_x (float): Size of the search box in the X dimension.
    - size_y (float): Size of the search box in the Y dimension.
    - size_z (float): Size of the search box in the Z dimension.
    - exhaustiveness (int): Exhaustiveness of the search (default is 8).

    Returns:
    - None
    """
    print(Usage)
    return

def pockets_with_fpocket_Help():
    Usage="""
    Identifies pockets on a protein surface using fpocket and displays them in PyMOL.

    Parameters:
    - receptor (str): Name of the receptor object in PyMOL.

    Returns:
    - None

    Description:
    This function uses fpocket to identify pockets on the surface of a protein loaded in PyMOL. It saves the
    receptor structure as a temporary PDB file, runs fpocket on this file, and loads the identified pockets back
    into PyMOL as separate objects. Pockets are then displayed as semi-transparent colored spheres, with each
    pocket numbered sequentially.
    
    Requirements:
    - PyMOL
    - fpocket (accessible from the command line)
    
    Usage Example:
        pockets_with_fpocket receptor_name

    The function will print fpocket’s output to the console and display each identified pocket as a separate
    PyMOL object with a unique color.
    """
    print(Usage)
    return
    
def align_with_rdkit(reference:str,target:str):
    
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
    os.remove(ref_tmp.name)
    os.remove(tar_tmp.name)
    
    return rmsd, mcs.smartsString, len(atom_map)


def align_with_fkcombu(receptor: str, reference: str, target: str):
    
    # Create temporary files for receptor and ligands in PDB and SDF formats
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as receptor_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as reference_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as target_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as out_tmp:

        # Save the receptor, reference, and target molecules from PyMOL to temporary files
        cmd.save(receptor_tmp.name, receptor, format="pdb")
        cmd.save(reference_tmp.name, reference, format="sdf")
        cmd.save(target_tmp.name, target, format="sdf")

        # Prepare the fkcombu command
        command = [
            "fkcombu",
            "-T", target_tmp.name,    # Target ligand file path
            "-R", reference_tmp.name,  # Reference ligand file path
            "-P", receptor_tmp.name,   # Receptor protein file path
            "-osdfT", out_tmp.name     # Output file path for the aligned ligand
        ]

        # Run the fkcombu command
        try:
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Output any messages from fkcombu
            print("fkcombu output:")
            print(result.stdout)
            if result.stderr:
                print("fkcombu errors:")
                print(result.stderr)
    
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running fkcombu: {e}")
            print("Standard output:", e.stdout)
            print("Standard error:", e.stderr)

        # Check if the output file was created and has content
        if os.path.exists(out_tmp.name) and os.path.getsize(out_tmp.name) > 0:
            # Load the aligned ligand back into PyMOL
            cmd.load(out_tmp.name, object="aligned_ligand", format="sdf")
            print("Aligned ligand loaded as 'aligned_ligand'")
        else:
            print("Error: The output file was not created or is empty.")

        # Clean up temporary files
        os.remove(receptor_tmp.name)
        os.remove(reference_tmp.name)
        os.remove(target_tmp.name)
        os.remove(out_tmp.name)


def docking_with_smina(receptor: str, ligand: str, center_x: float, center_y: float, center_z: float,
              size_x: float, size_y: float, size_z: float, exhaustiveness: int = 8):

    
    # Create temporary files for receptor and ligand in PDBQT format
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as receptor_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as ligand_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as out_tmp:
        
        # Export receptor and ligand from PyMOL to temporary PDBQT files
        cmd.save(receptor_tmp.name, receptor, format="pdb")
        cmd.save(ligand_tmp.name, ligand, format="sdf")

        # Prepare the command to run smina
        command = [
            "smina",
            "--receptor", receptor_tmp.name,
            "--ligand", ligand_tmp.name,
            "--out", out_tmp.name,
            "--center_x", str(center_x),
            "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--exhaustiveness", str(exhaustiveness)
        ]

        # Run the command and capture the output and errors
        try:
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Output the results to the console
            print("Smina output:")
            print(result.stdout)
            
            if result.stderr:
                print("Smina errors:")
                print(result.stderr)

            # Load the docking result back into PyMOL as a new object
            cmd.load(out_tmp.name, object="docked_ligand")
            print("Docking results loaded as 'docked_ligand'")

        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running smina: {e}")
            print("Standard output:", e.stdout)
            print("Standard error:", e.stderr)
    

def pockets_with_fpocket(receptor: str):
    # Create a temporary file for the receptor PDB
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as receptor_tmp:
        # Save the receptor from PyMOL to a temporary PDB file
        cmd.save(receptor_tmp.name, receptor, format="pdb")
        
        # Prepare the fpocket command
        command = ["fpocket", "-f", receptor_tmp.name]
        
        try:
            # Run the fpocket command
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Output the results to the console
            print("Fpocket output:")
            print(result.stdout)
            
            if result.stderr:
                print("Fpocket errors:")
                print(result.stderr)
                
            # Construct paths for the fpocket output directory and files
            pocket_output_dir = receptor_tmp.name.replace(".pdb", "_out")
            pocket_pdb_file = os.path.join(pocket_output_dir, os.path.basename(receptor_tmp.name).replace(".pdb", "_out.pdb"))
            info_file = os.path.join(pocket_output_dir, os.path.basename(receptor_tmp.name).replace(".pdb", "_info.txt"))

            # Load the fpocket output PDB file into PyMOL
            if os.path.exists(pocket_pdb_file):
                cmd.load(pocket_pdb_file,format="pdb",object="pockets_found")
                print("Pockets loaded")

                # Display fpocket info file content
                if os.path.exists(info_file):
                    with open(info_file) as pocket_data:
                        print(pocket_data.read())
                else:
                    print("Warning: Info file not found.")
            else:
                print("Error: Pocket PDB file not found.")
        
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running fpocket: {e}")
            print("Fpocket output:", e.stdout)
            print("Fpocket error:", e.stderr)
            return
        
        # Selecting and displaying pockets
        stored.list = []
        cmd.iterate("(resn STP)", "stored.list.append(resi)")
        lastSTP = stored.list[-1] if stored.list else None
        
        if lastSTP:
            for my_index in range(1, int(lastSTP) + 1):
                pocket_selection = f"pocket{my_index}"
                cmd.select(pocket_selection, f"resn STP and resi {my_index}")
                cmd.color(my_index+1, pocket_selection)
                cmd.show("spheres", pocket_selection)
                cmd.set("sphere_scale", 1, pocket_selection)
                cmd.set("sphere_transparency", 0.4, pocket_selection)
        cmd.remove("pockets_found and polymer")
        
        # Clean up temporary and output files
        os.remove(receptor_tmp.name)
        # Remove the output directory and its contents
        if os.path.exists(pocket_output_dir):
            shutil.rmtree(pocket_output_dir)


# Extend PyMOL's command set to include this function
cmd.extend("align_with_rdkit",align_with_rdkit)
cmd.extend("align_with_fkcombu", align_with_fkcombu)
cmd.extend("docking_with_smina", docking_with_smina)
cmd.extend("pockets_with_fpocket", pockets_with_fpocket)

