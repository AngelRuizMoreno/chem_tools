from pymol import cmd, stored
import subprocess
import tempfile


def __init__(self):
	self.menuBar.addcascademenu('Plugin','Smina','Smina',label = 'Smina')
	self.menuBar.addmenuitem('Smina', 'command','Help',label = 'Usage',command = lambda s=self : Help())

def Help():
    Usage="""
    Run smina with specified parameters, loading receptor and ligand directly from PyMOL session.

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
    
def run_smina(receptor: str, ligand: str, center_x: float, center_y: float, center_z: float,
              size_x: float, size_y: float, size_z: float, exhaustiveness: int = 8):
    """
    Run smina with specified parameters, loading receptor and ligand directly from PyMOL session.

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

# Extend PyMOL's command set to include this function
cmd.extend("smina_pymol", run_smina)
