import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor

# 1. MOLEKULA-ADATBÁZIS (Ezt a "könyvtárat" fogjuk átfuttatni a szűrőn)
compound_library = [
    {'name': 'Ecetsav', 'smiles': 'CC(=O)O'},
    {'name': 'Benzol', 'smiles': 'C1=CC=CC=C1'},
    {'name': 'Paracetamol', 'smiles': 'CC(=O)NC1=CC=C(O)C=C1'},
    {'name': 'Methanol', 'smiles': 'CO'},
    {'name': 'Hexán', 'smiles': 'CCCCCC'},
    {'name': 'Koffein', 'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'},
    {'name': 'Ibuprofén', 'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'}
]

# Betanító adatsor a QSAR modellhez
train_data = {
    'smiles': [
        'CC(=O)OC1=CC=CC=C1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'C1=CC=C(C=C1)O', 'CCO', 'CCCCCC', 'C1=CC=CC=C1'
    ],
    'solubility': [-2.24, -1.18, -3.32, -0.82, 0.5, -3.9, -2.1]
}

def smiles_to_fp(smiles: str, n_bits: int = 2048):
    """Morgan Fingerprint generálás."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits))

def run_virtual_screening(min_solubility: float = -2.5):
    """
    Virtuális szűrést végző függvény:
    Kiszűri azokat a molekulákat, amelyek várható oldékonysága jobb a megadott küszöbértéknél.
    """
    print("--- 1. QSAR Modell betanítása ---")
    df_train = pd.DataFrame(train_data)
    X_train = np.array([smiles_to_fp(s) for s in df_train['smiles']])
    y_train = df_train['solubility'].values

    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    print(" -> A modell készen áll a szűrésre.\n")

    print(f"--- 2. Virtuális szűrés indítása (Küszöbérték: logS >= {min_solubility}) ---")
    results = []

    for comp in compound_library:
        fp = smiles_to_fp(comp['smiles'])
        if fp is not None:
            pred_sol = model.predict(fp.reshape(1, -1))[0]
            
            # Elfogadási kritérium ellenőrzése
            passed = pred_sol >= min_solubility
            
            results.append({
                'Név': comp['name'],
                'SMILES': comp['smiles'],
                'Becsült logS': round(pred_sol, 2),
                'Státusz': 'ELFOGADVA' if passed else 'ELUTASÍTVA'
            })

    results_df = pd.DataFrame(results)
    print(results_df.to_string(index=False))

    # Elfogadott molekulák kimentése
    passed_compounds = results_df[results_df['Státusz'] == 'ELFOGADVA']
    print(f"\n--- 3. Eredmény ---")
    print(f"Összes szűrt molekula: {len(compound_library)}")
    print(f"Kritériumoknak megfelelt: {len(passed_compounds)}")

if __name__ == "__main__":
    # Futtassuk a szűrést pl. -2.5-ös küszöbértékre
    run_virtual_screening(min_solubility=-2.5)