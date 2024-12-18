import streamlit as st
from pathlib import Path

# For some reason the windows version only works if this is imported here
import pyopenms


if __name__ == "__main__":

    pages = {
        "StreamSage Web App": [
            st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="👋"),
            st.Page(
                Path("content", "documentation.py"), title="Documentation", icon="📖"
            ),
        ],
        "Proteomics Database Search": [
            st.Page(Path("content", "sageworkflow.py"), title="Sage", icon="🚀"),
        ],
    }

    pg = st.navigation(pages)
    pg.run()
