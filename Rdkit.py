import streamlit as st
#Devid the screen into two columns

col1, col2 = st.columns(2)

with col1:
    st.title('🎈 Installation')

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

with col1:
    st.title('🎈 Calculate QED')

    code = """
    
        import os
    from rdkit.Chem import rdBase, RDConfig
    from rdkit import Chem
    from rdkit.Chem import PandasTools
    from rdkit.Chem.Draw import IPythonConsole
    from rdkit.Chem.Descriptors import qed
    print( rdBase.rdkitVersion )
    
    sdfpath = os.path.join( RDConfig.RDDocsDir, "Book/data/cdk2.sdf" )
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

with col2:
    st.title('🎈 Get a RDKit molecule from SMILES')

    code = """  smiles = 'COC(=O)c1c[nH]c2cc(OC(C)C)c(OC(C)C)cc2c1=O'
    mol = Chem.MolFromSmiles(smiles)
    print(mol)"""
    st.code(code, language="python")
