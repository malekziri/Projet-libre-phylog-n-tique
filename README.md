### 01_prepare_snv.py — Préparation des SNV

Le script lit la feuille `subs` du fichier Excel et conserve les prélèvements
sélectionnés des patients A10, A17, A21, A24, A31 et A34.

Chaque SNV est identifié par son chromosome, sa position et les allèles
de référence et alternatif. La VAF du WGS est utilisée en priorité,
puis celle du panel de validation si la première est manquante. Une VAF strictement
positive est codée 1 ; une VAF nulle ou manquante est codée 0.
Tous les statuts de validation sont conservés.

Une ligne « Normal », entièrement codée 0, est ajoutée pour l’enracinement.

Les fichiers sont enregistrés dans `Data/processed/SNV/` :
- un fichier NEXUS par patient
- `SNV_summary.csv` : nombre de prélèvements, de SNV, de sites variables
  et de sites informatifs pour la parcimonie (calculés sans la ligne Normal)
- `SNV_presence.csv` : présence ou absence de chaque SNV dans chaque
  prélèvement, avec son annotation et son numéro de colonne dans le NEXUS.

Un site variable comporte les deux états. Un site informatif pour la
parcimonie comporte au moins deux prélèvements dans chacun des deux états.

Exécution depuis la racine du projet :
python3 scripts/01_prepare_snv.py



### 02_prepare_indel.py — Préparation des indels

Le script lit la feuille `indels` du fichier Excel et conserve les prélèvements
sélectionnés des patients A10, A17, A21, A24, A31 et A34.

Chaque indel est identifié par son chromosome, sa position et les allèles
de référence et alternatif. La VAF du WGS est utilisée en priorité,
puis celle du panel de validation si la première est manquante.  
Une VAF strictement positive est codée 1 ; une VAF nulle
ou manquante est codée 0. Tous les statuts de validation sont conservés.

Une ligne « Normal », entièrement codée 0, est ajoutée pour l’enracinement.

Les fichiers sont enregistrés dans `Data/processed/INDEL/` :
- un fichier NEXUS par patient 
- `INDEL_summary.csv` : nombre de prélèvements, d’indels, de sites variables
  et de sites informatifs pour la parcimonie, calculés sans la ligne Normal 
- `INDEL_presence.csv` : présence ou absence de chaque indel dans chaque
  prélèvement, avec son annotation et son numéro de colonne dans le NEXUS.

Exécution depuis la racine du projet :

```bash
python3 scripts/02_prepare_indel.py
```



### 03_prepare_rearrangement.py — Préparation des réarrangements

Le script lit la feuille `rearrangements` du fichier Excel pour les patients
A10, A17, A21, A24, A31 et A34.

Chaque réarrangement est identifié par les chromosomes, positions et
orientations de ses deux points de cassure. Les colonnes `samples`
et `tumour-counts` associent les prélèvements à leurs nombres de lectures
de soutien. Un réarrangement est codé 1 dans un prélèvement sélectionné
si ce nombre est strictement positif, et 0 sinon.

Une ligne « Normal », entièrement codée 0, est ajoutée pour l’enracinement.

Les fichiers sont enregistrés dans `Data/processed/REARRANGEMENT/` :
- un fichier NEXUS par patient 
- `REARRANGEMENT_summary.csv` : nombre de prélèvements, de réarrangements,
  de sites variables et de sites informatifs pour la parcimonie,
  calculés sans la ligne Normal 
- `REARRANGEMENT_presence.csv` : présence ou absence de chaque réarrangement
  dans chaque prélèvement, avec son numéro de colonne dans le NEXUS.

  La colonne `annotation` reprend le champ `paired-flag` du fichier source.

Exécution depuis la racine du projet :

```bash
python3 scripts/03_prepare_rearrangement.py
```



### 04_prepare_combined.py — Combinaison des trois types d’altérations

Le script utilise les fichiers de présence produits par les trois scripts
précédents pour les patients A10, A17, A21, A24, A31 et A34.

Les matrices sont réunies dans l’ordre SNV, indels, puis réarrangements,
en conservant le même ordre des prélèvements et la ligne « Normal ».
Le type d’altération est ajouté à chaque identifiant pour distinguer
les événements.

Les fichiers sont enregistrés dans `Data/processed/SNV_INDEL_REARRANGEMENT/` :
- un fichier NEXUS par patient pour une analyse sans partition 
- `SNV_INDEL_REARRANGEMENT_summary.csv` : nombre de prélèvements,
  nombre d’événements par type, nombre total de sites, de sites variables
  et de sites informatifs pour la parcimonie
- `SNV_INDEL_REARRANGEMENT_presence.csv` : présence ou absence de chaque
  événement, avec son type, son annotation et son numéro de colonne
  dans le NEXUS combiné.

Les statistiques des sites sont calculées sur les prélèvements tumoraux,
sans la ligne Normal.

Exécution depuis la racine du projet, après les trois scripts de préparation :

```bash
python3 scripts/04_prepare_combined.py
```



### 05_run_iqtree.sh — Reconstruction des arbres avec IQ-TREE

Le script lance IQ-TREE 2.4.0 pour les six patients retenus et les quatre
jeux de données : SNV, indels, réarrangements et matrice combinée,
soit 24 analyses.

Les matrices sont déclarées comme binaires (`-st BIN`).
L’option `-m MFP` compare les modèles et sélectionne automatiquement
celui ayant le BIC le plus faible. IQ-TREE reconstruit ensuite un arbre
de maximum de vraisemblance avec 1 000 réplicats d’ultrafast bootstrap
(`-B 1000`). Les arbres sont enracinés sur « Normal » (`-o Normal`).
Les prélèvements ayant des profils identiques sont conservés pendant
l’analyse (`-keep-ident`).

Les résultats sont enregistrés dans `Results/IQTREE/`, par type de données :
- `.iqtree` : rapport, comparaison des modèles et modèle sélectionné 
- `.treefile` : arbre de maximum de vraisemblance avec supports UFBoot 
- `.contree` : arbre de consensus bootstrap 
- `.log` : journal de l’analyse.

Exécution depuis la racine du projet :

```bash
bash scripts/05_run_iqtree.sh
```



### 06_run_mrbayes_combined.sh — Analyse bayésienne des matrices combinées

Le script analyse les matrices combinant SNV, indels et réarrangements
des patients A10, A17, A21, A24, A31 et A34 avec MrBayes.

Les fichiers NEXUS sont adaptés au modèle binaire de MrBayes
(`DATATYPE=RESTRICTION`). Le modèle utilise une vitesse d’évolution
identique entre sites (`RATES=EQUAL`) et des fréquences d’états estimées
avec une loi a priori de Dirichlet (1,1).

Le réglage `CODING` est déterminé sur la matrice complète, incluant Normal :
- `VARIABLE` si elle ne contient aucun site constant 
- `ALL` si des sites constants sont présents

Chaque patient fait l’objet de deux analyses indépendantes de
10 000 000 de générations, comprenant quatres arbres. 
Les arbres et paramètres sont enregistrés toutes les 1 000 générations. 
Les premiers 25 % des enregistrements sont écartés

Les arbres sont enracinés sur « Normal ». Les résultats sont enregistrés
dans `Results/MRBAYES/SNV_INDEL_REARRANGEMENT/`, notamment :
- les fichiers NEXUS utilisés par MrBayes 
- les arbres et paramètres échantillonnés pour chaque analyse 
- les arbres de consensus avec probabilités postérieures 
- les tableaux de synthèse et diagnostics permettant d’évaluer la convergence.

Exécution depuis la racine du projet :

```bash
bash scripts/06_run_mrbayes_combined.sh
```



### 07_run_scf.sh — Calcul des facteurs de concordance des sites

Le script calcule les facteurs de concordance des sites pour les six
patients et les quatre jeux de données analysés avec IQ-TREE.

L’arbre précédemment estimé est conservé sans nouvelle recherche
(`-te`). Le modèle sélectionné selon le BIC est récupéré dans le rapport
IQ-TREE et réutilisé pour le calcul.

L’option `--scfl 100` calcule le sCF fondé sur la vraisemblance,
à partir de 100 tirages de quatre prélèvements autour de chaque
branche interne. Le sCF mesure l’accord des sites avec les regroupements
de l’arbre et complète les supports ultrafast bootsrap.

Les résultats sont enregistrés dans `Results/SCF/`, par type de données :
- `.cf.tree` : arbre portant les valeurs de sCF ajoutées aux supports existants 
- `.cf.branch` : arbre avec les identifiants des branches 
- `.cf.stat` : statistiques de concordance pour chaque branche.

Exécution depuis la racine du projet, après `05_run_iqtree.sh` :

```bash
bash scripts/07_run_scf.sh
```



