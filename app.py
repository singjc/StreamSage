import streamlit as st
from pathlib import Path

pages = {
    "StreamSage Web App" : [
        st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="👋"),
        st.Page(Path("content", "documentation.py"), title="Documentation", icon="📖"),
    ],
    "Proteomics Database Search": [
        st.Page(Path("content", "sageworkflow.py"), title="Sage", icon="🚀"),
    ]
}

pg = st.navigation(pages)
pg.run()