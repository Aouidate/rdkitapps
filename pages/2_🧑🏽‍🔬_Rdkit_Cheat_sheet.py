import streamlit as st
#Devid the screen into two columns
st.markdown("<h1 style='text-align: center; color: grey;'>Rdkit Cheat Sheet</h1>", unsafe_allow_html=True)
#st.title("Rdkit Cheat Sheet")

col1, col2, = st.columns(2)

with col1:
    st.header('🎈 Installation')

    code = """
    
    #Using pip
    pip install rdkit

    #Using Anaconda
    conda install -c rdkit rdkit

    #Check the rdkit version
    dkit.__version__

    #Import rdkit
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import IPythonConsole

    """
    #Inpired from byiwatobipen
    st.code(code, language="python")

with col2:
    st.header('🎈 Calculate QED')

    code = """
    
        import os
    from rdkit.Chem import rdBase, RDConfig
    from rdkit import Chem
    from rdkit.Chem import PandasTools
    from rdkit.Chem.Draw import IPythonConsole
    from rdkit.Chem.Descriptors import qed
    print( rdBase.rdkitVersion )
    
    sdfpath = os.path.join( RDConfig.RDDocsDir, "path/to/sdffile.sdf" )
    mols = [ m for m in Chem.SDMolSupplier( sdfpath ) if m != None ]
    df = PandasTools.LoadSDF( sdfpath )
    print( len( mols ))
    
    df.head( 2 )
    
    df[ "QED" ] =  df.ROMol.apply( qed )
    df.head(2 )
    
    from rdkit.Chem import QED
    for mol in mols:
        print( QED.properties( mol ) )"""
    #Inpired from byiwatobipen
    st.code(code, language="python")

with col1:
    st.header('🎈 Get a RDKit molecule from SMILES')

    code = """  smiles = 'COC(=O)c1c[nH]c2cc(OC(C)C)c(OC(C)C)cc2c1=O'
    mol = Chem.MolFromSmiles(smiles)
    print(mol)"""
    st.code(code, language="python")

    st.header('🎈 Reading sets of molecules')

    code = """  
    file_name = 'somesmilesfile.smi'
    with open(file_name, "r") as ins:
    smiles = []
    for line in ins:
        smiles.append(line.split('\n')[0])
    print('# of SMILES:', len(smiles))"""
    st.code(code, language="python")

with col2:
    st.header('🎈 Molecular fingerprints')

    code = """ 
    #Import libraries
    from rdkit import DataStructs
    from rdkit.Chem.Fingerprints import FingerprintMols

    set mols
    m1 = Chem.MolFromSmiles("CCCCO")
    m2 = Chem.MolFromSmiles("c1ccccc1CO")

    #Calculate fps
    fp1 = FingerprintMols.FingerprintMol(m1)
    fp2 = FingerprintMols.FingerprintMol(m2)

    #Get similarity
    DataStructs.FingerprintSimilarity(fp1, fp2)"""

    st.code(code, language="python")

with col1:
    st.header("🎈 Molecular Descriptors 'MDs'")

    code = """ 
        #Import libraries
    from rdkit import Chem

    #set mols
    m1 = Chem.MolFromSmiles("CCCCO")

    #Calculate MDs
    print ('H-bond donors:', Descriptors.NumHDonors(m1))
    print ('Molecular weight:', Descriptors.MolWt(m1))
    print ('LogP:', Descriptors.MolLogP(m1))
    """

    st.code(code, language="python")
with col2:
    st.header("🎈 Draw molecules ")

    code = """ 
        #Import libraries
    from rdkit import Chem
    from rdkit import Draw

    #read mols 

    sdfpath = os.path.join( RDConfig.RDDocsDir, "path/to/sdffile.sdf" )
    mols = [ m for m in Chem.SDMolSupplier( sdfpath ) if m != None ]

    #Draw mols
    Draw.MolsToGridImage(mols, molsPerRow=8)
    """

    st.code(code, language="python")
