import streamlit as st
from PIL import Image

st.markdown("<h2 style='text-align: center; color: grey;'>Exploring Various Chemoinformatic Applications and Codes for Drug Discovery</h2>", unsafe_allow_html=True)

st.markdown("""<div style='text-align: justify; color: black;'>Welcome to our Streamlit platform, where you can explore the fascinating world of 
chemoinformatics with ease. Our multi-app platform provides you with an all-in-one solution for all your chemoinformatics needs, including QSAR 
modeling of Ames test, a comprehensive Rdkit cheat sheet, molecular visualization, and much more - all in a single, seamless interface.


As an avid chemoinformatics enthusiast, I'm thrilled to be sharing with you exciting blogs, codes, and articles through this small blog. 
With my help, you'll be able to discover the incredible world of chemoinformatics, starting with our Streamlit multi-web app. It features 
an array of apps and code snippets that will help you dive deep into the world of chemoinformatics.
 </div>""", unsafe_allow_html=True)
#st.markdown("Explore the fascinating world of chemoinformatics 🖥️ 🧪 through our multi-app Streamlit platform, where you can dive into diverse things such as QSAR modeling of Ames test, Rdkit cheat sheet, molecular visualization, and more, all in one seamless interface. I am Adnane Aouidate and I will be sharing with you some exciting blogs, codes, articles and also help you to discover the chemoinformatics' world: Here is the screenshot of a streamlit web app that I started to develop which will find many apps and code snippets I am an organic chemist by training, I statred my journey in molecular modelling and CADD from 2014 during my thesis since 2019 I'm working as postdoctoral researcher (chemoinformatician).")

st.markdown("""
    I share some code snippets about chemoinformatics, Drug Design, Data Science, AI  and related topics. 
    The posts and codes will mostly be shared on github and linkedin platforms.
    
    You can find my scientific articles on my [google scholar](https://scholar.google.com/citations?user=Yngy4o4AAAAJ&hl=en)
    and if you would like to check when I publish new codes, you can 
    visit my github repos : 
    [Github](https://github.com/Aouidate/).

    Below is a github repo you might find interesting...
""")

with st.container():
    image_col, text_col = st.columns((1,2))
    with image_col:
        image = Image.open('ADN_Chemoiinformatics.png')
        st.image(image, use_column_width=True)

    with text_col:
        st.subheader("Chemoinformatics Playground: Exploring the World of Molecular Data Analysis and Modeling")
        st.write("""Compilation of chemoinformatics and machine learning techniques for drug discovery. 
        In the present repository, you can find a set of tutorials that are available for every one that is passionate about data science, 
        chemoinformatics or drug discovery. The main objective is to create a compilation of different 
        tutorials with explainations for the different required techniques to ghuide newbies in this field of research.
            """)
        st.markdown("[Read more...](https://github.com/Aouidate/Chemoinformatics-tutos/tree/master)")
