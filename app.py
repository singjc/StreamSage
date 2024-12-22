import streamlit as st
from pathlib import Path

# For some reason the windows version only works if this is imported here
import pyopenms


if __name__ == "__main__":

    pages = {
        "Proteomics Database Search": [
            st.Page(Path("content", "quickstart_sage.py"), title="Quickstart", icon="👋"),
            st.Page(Path("content", "sageworkflow.py"), title="Sage", icon="🚀"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()
