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

st.title("AmesTesto")

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

            row = np.array([desc_MolLogP,
                        desc_TPSA])

            if(i==0):
                baseData=row
            else:
                baseData=np.vstack([baseData, row])
            i=i+1

        columnNames=["MolLogP", "TPSA"]
        descriptors = pd.DataFrame(data=baseData,columns=columnNames)

        return descriptors
    except:
         st.error('Please enter a valid structure')
######################
# Page Title
######################

image = Image.open('./images/AMES_logo.jpg')

st.image(image, use_column_width=True)

st.markdown("<h2 style='text-align: justify; color: black;'> AmesTesto: Accurately Predict Molecular Mutagenicity with Adnane's Webapp !</h2>", unsafe_allow_html=True)


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
#X[1:] # Skips the dummy first item

######################
# Pre-built model
######################

# Reads in saved model
load_model = pickle.load(open('Ames_calssifcator.pkl', 'rb'))

# Apply model to make predictions
try:
    prediction = load_model.predict(X)
    prediction = pd.DataFrame(prediction, columns = ["Amest test"], index = SMILES)
except:
    pass
#prediction_proba = load_model.predict_proba(X)

st.header('Predicted Mutagenicity')
try:
    prediction[1:] # Skips the dummy first item 
except:
    pass