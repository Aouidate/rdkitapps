import numpy as np
import pandas as pd
from rdkit import Chem
from PIL import Image
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
import os
from rdkit.Chem import rdBase, RDConfig
from rdkit import Chem
from rdkit.Chem import QED
import streamlit as st


######################
# Custom function
######################
## Calculate molecular descriptors

st.title("Drug Likness Filters")

def Calculate (smiles, verbose = False):
    mols = []
    for sm in smiles:
        mol = Chem.MolFromSmiles(sm)
        if mol is not None:
            mols.append(mol)
        else :
            mols.append(None)
            print(f"Invalid structure: {sm}")
    baseData= np.arange(1,1)
    i=0
    for mol in mols:
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_MW = Descriptors.MolWt(mol)
        desc_NHA = Descriptors.NumHAcceptors(mol)
        desc_NHD =Descriptors.NumHDonors(mol)
        desc_TPSA = Descriptors.TPSA(mol)
        desc_RB = Descriptors.NumRotatableBonds(mol)
        desc_QED = QED.qed(mol) 
        row = np.array([desc_MolLogP, desc_MW, desc_NHA, desc_NHD,
                       desc_TPSA, desc_RB, desc_QED])
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1
    columnNames=["LogP", "MW", "nHA", "nHD", "TPSA", "RB", "QED"]

    descriptors = pd.DataFrame(data=baseData,columns=columnNames, index = smiles ) 

    return descriptors
#Define drug like rules

def Veber(row):

    conditions = [row['RB'] <= 10,
                    row['TPSA'] <= 140]
    true_values = sum(conditions)
    if true_values >= 2:
        return 'Yes'
    else:
        return 'No'

def Lipinski(row):

    conditions = [row['LogP'] <= 5,
                    row['MW'] <= 500,
                    row['nHA'] <= 10,
                    row['nHD'] <= 5,]
    true_values = sum(conditions)
    if true_values >= 4:
        return 'Yes'
    else:
        return 'No'

def Rule_of_3(row):

    conditions = [row['LogP'] <= 3,
                    row['MW'] <= 300,
                    row['nHA'] <= 3,
                    row['nHD'] <= 3,
                    row['RB'] <= 3]
    true_values = sum(conditions)
    if true_values >= 5:
        return 'Yes'
    else:
        return 'No'

######################
# Page Title
######################

image = Image.open('atome.png',  width=400)

st.image(image)

st.write("""# Molecular descriptors Calculation!""")
st.markdown("""<div style='text-align: justify; color: black;'> Looking for an innovative solution to predict drug likeness for your molecules? 
Look no further! Our web app takes as input a SMILES string and provides accurate predictions 
of molecular descriptors and drug likeness properties. This simple tool offers a support for 
various chemical features and is built to offer fast and easy predictions for both novice and expert users. 
Say goodbye to complex molecular models and calculations, and take advantage of our solution 
– designed to improve your drug discovery processes! </div>""", unsafe_allow_html=True)


######################
# Input molecules (Side Panel)
######################

st.sidebar.header('User Input Molecules')

## Read SMILES input
SMILES_input = "CCCCO\nc1ccccc1\nCN"

SMILES = st.sidebar.text_area("SMILES input", SMILES_input)
SMILES = "C\n" + SMILES #Adds C as a dummy, first item
SMILES = SMILES.split('\n')


## Calculate molecular descriptors
st.header('Computed molecular descriptors')
try :
    df = Calculate(SMILES)
    df["Lipinski"] = df.apply(Lipinski, axis=1)
    df["Veber"] = df.apply(Veber, axis=1)
    df["Rule of 3"] = df.apply(Rule_of_3, axis=1)
except:
    print("Invalid structure")

#df["Molecule"]= df.index.apply(Chem.MolFromSmiles) #Get molecule from smiles

st.header('Drug-likeness filters')
try:
    df[1:] # Skips the dummy first item 
except:
    st.error('Please enter a valid structure ')

