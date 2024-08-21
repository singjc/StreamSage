import streamlit as st

from src.common import page_setup, save_params, show_table
from src import sageworkflow

# Page name "workflow" will show mzML file selector in sidebar
params = page_setup()

cols = st.columns(3)
cols[1].image("https://github.com/lazear/sage/blob/master/figures/logo.png?raw=true")
st.markdown("Sage, a proteomics database search engine, is an open-source, high-performance proteomics search engine designed to accelerate workflows by integrating features such as retention time prediction, quantification, peptide-spectrum match rescoring, and false discovery rate control, all while providing rapid and scalable searching capabilities for mass spectrometry data. For more information, visit the [Sage GitHub repository](https://github.com/lazear/sage)")

wf = sageworkflow.SageWorkflow()

t = st.tabs(["📁 **File Upload**", "⚙️ **Configure**", "🚀 **Run**", "📊 **Results**"])
with t[0]:
    wf.show_file_upload_section()
    
with t[1]:
    wf.show_parameter_section()
    
with t[2]:
    wf.show_execution_section()
    
with t[3]:
    wf.show_results_section()