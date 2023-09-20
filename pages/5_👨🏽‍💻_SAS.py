import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from rdkit import Chem, DataStructs
#from rdkit.Chem import PandasTools
from rdkit.Chem import rdFMCS, AllChem
import streamlit as st
import base64
import io
from io import BytesIO

st.markdown("<h2 style='text-align: center; color: grey;'>Structure-Activity Similarity (SAS) map</h2>", unsafe_allow_html=True)

st.markdown("<div style='text-align: justify; color: black; font-size: large'>SALI (Structure Activity Landscape Index) is a technique that was published by Rajarshi Guha and Jonn Van Drie In 2008 that allows for the easy detection of Activity cliffs (ACs) compound pairs. These pairs are unique in that a small change to their chemical structure can result in a significant difference in their physical properties or biological activity. Frequently, these modifications can greatly assist us in identifying the key functional components of the molecule that significantly contribute to its activity. This web application allows users to scan Structure-Activity Relationships using Structure-Activity Similarity. It generates a Structure-Activity Similarity (SAS) map, enabling the identification of molecules with the highest SALI values. </div>", unsafe_allow_html=True)

st.set_option('deprecation.showfileUploaderEncoding', False)
#st.markdown("<h2 style='text-align: center; color: grey;'>Structure-Activity Similarity (SAS) maps</h2>", unsafe_allow_html=True)

st.markdown("<p style='text-align: justify;'>   Note: The current version uses the Extended-connectivity fingerprints (ECFPs) with Radius 4 </p>", unsafe_allow_html=True)
st.sidebar.header('*User Input Parameters*')

def user_input_features_1():
    a = st.sidebar.slider('Activity Difference', 2.0,3.0,2.5)
    data = {'Activity Difference': a}
    features_1 = float(a)
    return features_1
a= user_input_features_1()

def user_input_features_2():
    b = st.sidebar.slider('Similarity threshold', 0.5,1.0,0.9)
    data = {'Similarity threshold': b}
    features_2 = float(b)
    return features_2
b = user_input_features_2()

def calc(df, a=a, b=b, to_plot = True):
    try:
        df['mol'] = [Chem.MolFromSmiles(x) for x in df.smiles]
        df['fp'] = [Chem.RDKFingerprint(x) for x in df.mol]
        pIC50_list = df.pIC50.values
        smiles_list = df.smiles.values
        fp_list = df.fp.values
        sal_list = []
        for i,fp in enumerate(df.fp):
            Similarity = DataStructs.BulkTanimotoSimilarity(fp,fp_list)
            for j in range(0,i):
                ii,jj = i,j
                if pIC50_list[i] >= pIC50_list[j]:
                    ii,jj = i,j
                else:
                    jj,ii = i,j
                delta = abs(pIC50_list[i]-pIC50_list[j])
                sim = Similarity[j]
                sal_list.append([ii,jj,sim,delta/(1-sim+0.001),delta, smiles_list[ii],smiles_list[jj],pIC50_list[ii],pIC50_list[jj]])
        sal_df = pd.DataFrame(sal_list,columns=["i","j","Similarity","SAL","delta","SMILES_i","SMILES_j","pIC50_i","pIC50_j"])
        sal_df.rename(columns={"SAL":"SALI","delta":"Delta"}, inplace=True)
        sal_df.sort_values("SALI",ascending=False,inplace=True)
        sal_df = sal_df[sal_df["Similarity"] !=1]
        clifs = sal_df.loc[(sal_df["Similarity"]>=b) & (sal_df["Delta"]>=a)]

        if to_plot:
                fig, ax=plt.subplots(figsize=(10, 8))
                im =ax.scatter(data=sal_df, x="Similarity", y="Delta", s = "SALI", c ="SALI", cmap= "jet")
                fig.colorbar(im,ax=ax, label = "SALI value")
                plt.axhline(y=a) 
                plt.axvline(x=b) 
                plt.xlabel('Structural Similarity')
                plt.ylabel('Activity Difference')
                plt.show()
                st.pyplot(plt)
                st.markdown(imagedownload (plt,'sali.png'), unsafe_allow_html=True)
                st.write("Here, we show up to 6 top compound pairs that may present an activity cliffs:" )
                st.dataframe(clifs.head(6))
        return clifs.head(6)
    except (ZeroDivisionError, UnboundLocalError, ValueError, AttributeError) :
        print('Please upload your dataset in CSV format')
    pass

def convert_df(clifs):
   return clifs.to_csv(index=False).encode('utf-8')
try:
    csv = convert_df(clifs)
    st.download_button(
    "Press to Download",
    csv,
    "file.csv",
    "text/csv",
    key='download-csv')
except:
    pass
def imagedownload (plt, filename):
    s = io.BytesIO()
    plt.savefig(s, format = 'png', bbox_inches = 'tight')
    plt.close()
    b64 = base64.b64encode(s.getvalue()).decode()
    href = f'<a href="data:image/png;base64,{b64}" download={filename}>Download {filename} File</a>'
    return href

uploaded_file = st.sidebar.file_uploader("""Upload your dataframe""")
if uploaded_file is not None:
  uploaded_file = pd.read_csv(uploaded_file, index_col=1, encoding = "ISO-8859-1")
  #st.write(uploaded_file.head())

calc(df=uploaded_file, to_plot=True)
# if st.button('Calculate AD '):
#williams_plot(X_train = uploaded_file1 , X_test = uploaded_file2, Y_true_train = uploaded_file3, Y_true_test = uploaded_file4, y_pred_train = uploaded_file5,y_pred_test = uploaded_file6, toPrint =True,toPlot=True)
   #st.write(results)
st.markdown("""*[The application was inspired by the works of Pat Walters.](http://practicalcheminformatics.blogspot.com/)*""")
st.markdown("""*[For more information about the method click here](https://pubs.acs.org/doi/10.1021/ci7004093)*""")