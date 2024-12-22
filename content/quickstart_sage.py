from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components

from src.common.common import page_setup, v_space

page_setup(page="main")

st.markdown("# 👋 Quick Start")
st.markdown("## StreamSage: A Streamlit-based GUI for Sage using the **[OpenMS WebApp Template](https://github.com/OpenMS/streamlit-template/tree/main)**")
c1, c2 = st.columns(2)
c1.markdown(
    """
## ⭐ Features
       
- Upload (200Mb limit) raw data files or use local raw data files
- Configure Sage parameters
- Execute Sage in a background process
- Visualize results
"""
)
v_space(1, c2)
c2.image("assets/pyopenms_transparent_background.png", width=300)
if Path("OpenMS-App.zip").exists():
    st.subheader(
        """
Download the latest version for Windows here by clicking the button below.
"""
    )
    with open("OpenMS-App.zip", "rb") as file:
        st.download_button(
            label="Download for Windows",
            data=file,
            file_name="OpenMS-App.zip",
            mime="archive/zip",
            type="primary",
        )
    st.markdown(
        """
Extract the zip file and run the installer (.msi) file to install the app. The app can then be launched using the corresponding desktop icon.
"""
    )
    
st.markdown("## 📖 SageDocumentation")
st.markdown(
    f"""
Sage has great documentation. You can find it [here](https://sage-docs.vercel.app/docs) or see below.
"""
)
iframe_src = "https://sage-docs.vercel.app/docs"
components.iframe(iframe_src, height=1500, scrolling=True)
