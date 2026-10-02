# Automatizált Kemoinformatikai Keretrendszer

**Gépészmérnöki Projektfeladat**  
**Projekt témája:** 3D molekulamodellezés, QSAR gépi tanulási becslés és virtuális szűrés Python alapon  

---

## A Projekt Leírása

A projekt célja egy olyan moduláris, Python-alapú szoftveres adatfeldolgozó lánc létrehozása, amely a kémiai analízis és szerkezetvizsgálat három alapvető lépését automatizálja:

1. **3D Molekulamodellezés (`src/visualization.py`):**  
   Szöveges kémiai reprezentációkból (SMILES kódokból) állít elő molekuláris erőtérrel (MMFF) optimált, szabványos 3D koordinátákat (PDB fájlokat), majd ezeket interaktív, böngészőben forgatható 3D modellként jeleníti meg.

2. **QSAR Tulajdonság-becslés (`src/qsar_model.py`):**  
   A molekulák szerkezetéből számszerűsíthető Morgan-ujjlenyomatokat (szerkezeti vektorokat) képez, és egy betanított Random Forest (Véletlen Erdő) regressziós gépi tanulási modell segítségével becsli meg a fizikai-kémiai tulajdonságokat (pl. vízoldékonyság – logS).

3. **Virtuális Szűrés (`src/virtual_screening.py`):**  
   Egy teljes molekula-könyvtárat automatikusan kiértékel a betanított QSAR modell alapján, és a megadott küszöbértékek szerint osztályozza a vegyületeket (ELFOGADVA / ELUTASÍTVA).

---

## Projekt Felépítése

```text
kemo_projekt/
│
├── data/                    # Generált 3D PDB fájlok tárolója
├── src/                     # A modulok forráskódjai
│   ├── visualization.py     # 1. Modul: 3D megjelenítés
│   ├── qsar_model.py        # 2. Modul: QSAR gépi tanulási modell
│   └── virtual_screening.py # 3. Modul: Virtuális szűrés
│
├── .gitignore               # Átmeneti fájlok kiszűrése Githez
├── requirements.txt         # Projekt függőségei
└── README.md                # Részletes dokumentációGitHub Markdown Preview