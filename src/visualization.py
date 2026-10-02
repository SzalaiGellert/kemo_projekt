import os
import webbrowser
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_and_show_3d(smiles_code: str, filename: str = "molecule.pdb"):
    """SMILES-ből 3D modellt generál, elmenti és AZONNAL megjeleníti egy ablakban."""
    print(f"1. SMILES kód feldolgozása: {smiles_code}")
    mol = Chem.MolFromSmiles(smiles_code)
    
    if mol is None:
        print("Hiba: Érvénytelen SMILES kód!")
        return

    # 3D koordináták generálása
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    
    # Mentés PDB fájlba
    Chem.MolToPDBFile(mol, filename)
    print(f"2. 3D modell generálva: {filename}")

    # --- AUTOMATIKUS 3D MEGJELENÍTÉS ---
    # Létrehozunk egy apró interaktív 3D nézegető fájlt
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    </head>
    <body style="margin: 0; padding: 0; overflow: hidden; background-color: #111;">
        <div id="container" style="width: 100vw; height: 100vh;"></div>
        <script>
            let viewer = $3Dmol.createViewer("container", {{backgroundColor: "black"}});
            let pdbData = `{Chem.MolToPDBBlock(mol)}`;
            viewer.addModel(pdbData, "pdb");
            viewer.setStyle({{}}, {{stick: {{radius: 0.15}}, sphere: {{scale: 0.25}}}});
            viewer.zoomTo();
            viewer.render();
        </script>
    </body>
    </html>
    """
    
    view_file = "view_3d.html"
    with open(view_file, "w") as f:
        f.write(html_content)
        
    # Azonnal megnyitja a rendszer saját ablakában!
    webbrowser.open('file://' + os.path.realpath(view_file))
    print("3. A 3D ablak automatikusan megnyílt a képernyőn!\n")

if __name__ == "__main__":
    # Futtasd le ezt:
    generate_and_show_3d("CC(=O)OC1=CC=CC=C1C(=O)O") # Aszpirin