import numpy as np
import pandas as pd
from rdkit import Chem
from PIL import Image
import joblib
import pickle
from rdkit.Chem import Descriptors
import streamlit as st


######################
# Custom function
######################
## Calculate molecular descriptors

st.title("ADMETesto")

def AromaticProportion(mol):
  aromatic_atoms = [mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(mol.GetNumAtoms())]
  aa_count = [] 
  for i in aromatic_atoms:
    if i==True:
      aa_count.append(1)
  AromaticAtom = sum(aa_count)
  HeavyAtom = Descriptors.HeavyAtomCount(mol)
  AR = AromaticAtom/HeavyAtom
  return AR

def Calculate (smiles, verbose = False):
    try :
        mols = []
        for sm in smiles:
            mol = Chem.MolFromSmiles(sm)
            mols.append(mol)
        baseData= np.arange(1,1)
        i=0

        for mol in mols:

            desc_MolLogP = Descriptors.MolLogP(mol)
            desc_TPSA = Descriptors.TPSA(mol)
            desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
            desc_AromaticProportion = AromaticProportion(mol)
            row = np.array([desc_MolLogP,
                        desc_TPSA,desc_NumRotatableBonds,
                        desc_AromaticProportion])

            if(i==0):
                baseData=row
            else:
                baseData=np.vstack([baseData, row])
            i=i+1

        columnNames=["MolLogP", "TPSA","NumRotatableBonds","AromaticProportion"]
        descriptors = pd.DataFrame(data=baseData,columns=columnNames)

        return descriptors
    except:
         st.error('Please enter a valid structure')
######################
# Page Title
######################

image = Image.open('./images/AMES_logo.jpg')

st.image(image, use_column_width=True)

st.markdown("<h2 style='text-align: justify; color: black;'> ADMETesto: Accurately Predict Molecular Pharmacokenitc Properties of Small Organic Molecules !</h2>", unsafe_allow_html=True)


######################
# Input molecules (Side Panel)
######################

st.sidebar.header('User Input Features')

## Read SMILES input
SMILES_input = "CCCCO\nc1ccccc1\nCN"

SMILES = st.sidebar.text_area("SMILES input", SMILES_input)
SMILES = "C\n" + SMILES #Adds C as a dummy, first item
SMILES = SMILES.split('\n')
#st.header('Model information')
#SMILES[1:] # Skips the dummy first item

## Calculate molecular descriptors
#st.header('Computed molecular descriptors')

X = Calculate(SMILES)
#st.write(X)
#X[1:] # Skips the dummy first item

######################
# Pre-built model
######################

# Reads in saved model
ames_model = pickle.load(open('Ames_calssifcator.pkl', 'rb'))
sol_model = pickle.load(open('solubility_model.pkl', 'rb'))

# Apply model to make predictions
try:
    prediction_sol = sol_model.predict(X)
    prediction_ames = ames_model.predict(X.loc[:, ["MolLogP", "TPSA"]])
    predictions = list(zip(prediction_sol ,prediction_ames))
    prediction = pd.DataFrame(predictions, columns = ["ESOL","Amest test"], index = SMILES)
except:
    pass
#prediction_proba = load_model.predict_proba(X)

st.header('Predicted ADME values')
try:
    prediction[1:] # Skips the dummy first item 
except:
    pass