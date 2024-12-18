
from typing import Any
import streamlit as st
import time
from streamlit.components.v1 import html
import pyopenms as poms

def load_fasta():
    """
    Load the FASTA database file into the session state.
    """
    
    entries = []
    f = poms.FASTAFile()
    f.load(st.session_state["sage_config"]['database']['fasta'], entries)
    if "fasta_database" not in st.session_state:
        st.session_state["fasta_database"] = entries