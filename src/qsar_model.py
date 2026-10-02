import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error

# 1. KÉMIAI ADATSOR LÉTREHOZÁSA (Minta adatsor oldékonyságra - logS)
# A logS a vízoldékonyságot jelöli (minél kisebb/negatívabb, annál rosszabbul oldódik)
data = {
    'smiles': [
        'CC(=O)OC1=CC=CC=C1C(=O)O',  # Aszpirin
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', # Koffein
        'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', # Ibuprofén
        'C1=CC=C(C=C1)O',            # Fenol
        'CCO',                       # Etanol
        'CCCCCC',                    # Hexán
        'C1=CC=CC=C1'                # Benzol
    ],
    'solubility': [-2.24, -1.18, -3.32, -0.82, 0.5, -3.9, -2.1] # logS értékek
}

def smiles_to_fingerprint(smiles: str, n_bits: int = 2048):
    """SMILES kódból Morgan Fingerprint (ujjlenyomat) számsort készít."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Morgan Fingerprint generálása (radius=2 a standard választás)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
    return np.array(fp)

def build_and_train_qsar():
    print("1. Adatok előkészítése...")
    df = pd.DataFrame(data)
    
    # Fingerprintek kiszámítása minden molekulára
    X_list = []
    y_list = []
    
    for idx, row in df.iterrows():
        fp = smiles_to_fingerprint(row['smiles'])
        if fp is not None:
            X_list.append(fp)
            y_list.append(row['solubility'])
            
    X = np.array(X_list)
    y = np.array(y_list)
    
    print(f" -> Adatmátrix mérete: {X.shape} (Molekulák száma: {X.shape[0]}, Jellemzők: {X.shape[1]})")
    
    # 2. QSAR Gépi tanulási modell tanítása (Random Forest)
    print("2. QSAR (Random Forest) modell tanítása...")
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X, y)
    
    # 3. Modell értékelése a meglévő adatokon
    predictions = model.predict(X)
    r2 = r2_score(y, predictions)
    print(f" -> Modell pontossága (R² score): {r2:.2f}")
    
    # 4. ÚJ, ISMERETLEN MOLEKULA TULAJDONSÁGÁNAK BECSLÉSE
    new_smiles = "CC(=O)O" # Ecetsav
    new_fp = smiles_to_fingerprint(new_smiles).reshape(1, -1)
    predicted_solubility = model.predict(new_fp)[0]
    
    print("\n--- TESZT BECSLÉS ---")
    print(f"Új molekula (Ecetsav) SMILES: {new_smiles}")
    print(f"Becsült oldékonyság (logS): {predicted_solubility:.2f}")

if __name__ == "__main__":
    build_and_train_qsar()