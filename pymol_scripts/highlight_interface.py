from pymol import cmd, stored

def highlight_interface(protein, ligand, distance=4.0):
    """
    Highlights residues on a protein that are in close contact with a ligand.
    
    Parameters:
        protein (str): PyMOL object name for the protein.
        ligand (str): PyMOL object name for the ligand.
        distance (float): Maximum distance to consider an interaction.
    """
    cmd.select("interface", f"({protein}) within {distance} of ({ligand})")
    cmd.show("sticks", "interface")
    cmd.color("orange", "interface")

# Extend PyMOL's command set to include this function
cmd.extend("highlight_interface", highlight_interface)