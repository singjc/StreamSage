import time

import numpy as np
import pandas as pd
import streamlit as st

from pathlib import Path
import json

import plotly.express as px
from .view import load_ms_file, msspectrum_get_df, view_identifications, view_quantification, get_theo_spectrum, SpectrumAlignment, highlight_peptides
from .common import show_fig
from .workflow.WorkflowManager import WorkflowManager


@st.cache_data
def generate_random_table(x, y):
    """Example for a cached table"""
    df = pd.DataFrame(np.random.randn(x, y))
    time.sleep(2)
    return df

class SageWorkflow(WorkflowManager):
    # Setup pages for upload, parameter, execution and results.
    # For layout use any streamlit components such as tabs (as shown in example), columns, or even expanders.
    def __init__(self) -> None:
        # Initialize the parent class with the workflow name.
        super().__init__("Sage Workflow", st.session_state["workspace"])
        
        if "workflow_dir" not in st.session_state:
            st.session_state["workflow_dir"] = self.workflow_dir
        
    def upload(self) -> None:
        t = st.tabs(["MS data", "FASTA database", "Sage Config"])
        with t[0]:
            # Use the upload method from StreamlitUI to handle mzML file uploads.
            self.ui.upload_widget(
                key="mzML-files",
                name="MS data",
                file_types=["mzML", "mzXML", "d"],
                fallback=[str(f) for f in Path("example-data", "mzML").glob("*.mzML")],
            )
            
        with t[1]:
            # Use the upload method from StreamlitUI to handle FASTA database uploads.
            self.ui.upload_widget(
                key="fasta_database",
                name="FASTA database",
                file_types=["fasta"],
                fallback=[str(f) for f in Path("example-data", "fasta").glob("*.fasta")],
            )
            
        with t[2]:
            # Use the upload method from StreamlitUI to handle Sage config uploads.
            self.ui.upload_widget(
                key="sage-config",
                name="Sage Config",
                file_types=["json"],
                fallback=[str(f) for f in Path("example-data", "sage_config_template.json").glob("*.json")],
            )
            
    @st.experimental_fragment
    def configure(self) -> None:
        # Allow users to select mzML files for the analysis.
        self.ui.select_input_file("mzML-files", multiple=True)

        # Create tabs for different analysis steps.
        t = st.tabs(
            ["**Sage**", "**SageAdapter**"]
        )
        with t[0]:
            
            self.ui.input_exec('sage', Path('./assets/sage_config_template.json').resolve(), 4)
        
        # with t[1]:
        #     try:
        #         self.ui.input_TOPP(
        #             "SageAdapter",
        #             custom_defaults={},
        #         )
        #     except FileNotFoundError as e:
        #         st.error(f"An error occurred while trying to configure SageAdapter: {e}")
            
    @st.experimental_fragment
    def execution(self) -> None:
        # Any parameter checks, here simply checking if mzML files are selected
        if not self.params["mzML-files"]:
            self.logger.log("ERROR: No mzML files selected.")
            return

        # Get mzML files with FileManager
        in_mzML = self.file_manager.get_files(self.params["mzML-files"])

        # Log any messages.
        self.logger.log(f"Number of input mzML files: {len(in_mzML)}")

        # Run FeatureFinderMetabo tool with input and output files.
        self.logger.log("Detecting features...")
        
        st.session_state["sage_config"]["mzml_paths"] = in_mzML
        st.session_state["sage_config"]["output_directory"] = str(Path(self.file_manager.workflow_dir, "results"))
        
        # Write out config to json file
        with open(Path(self.file_manager.workflow_dir, "sage_config.json"), "w") as json_file:
            json.dump(st.session_state["sage_config"], json_file, indent=2)
        
        self.executor.run_exec(st.session_state['sage-exec-path'], str(Path(self.file_manager.workflow_dir, "sage_config.json")), st.session_state['batch-size'])
        
    @st.experimental_fragment
    def results(self) -> None:
        @st.experimental_fragment
        def show_consensus_features():
            df = pd.read_csv(file, sep="\t", index_col=0)
            df['label_mapped'] = np.where(df['label'] == 1, 'Target', 'Decoy')
            st.metric("number of PSMs", df.shape[0])
            
            c1, c2 = st.columns(2)
            c1.dataframe(df)
            
            color_blind_friendly_colors = ['#0072B2', '#D55E00']  

            fig = px.histogram(df, x='sage_discriminant_score', color='label_mapped', 
                            histnorm='density', 
                            title='Density Histogram of Discriminant Scores',
                            labels={'sage_discriminant_score': 'Discriminant Score', 'label_mapped': 'Label'},
                            barmode='overlay', 
                            color_discrete_sequence=color_blind_friendly_colors)

            fig.update_traces(opacity=0.75)  
            with c2:
                show_fig(fig, "Hyperscore Density Plot")
            
            st.title("PSMs with q-values <= 0.01")
            
            # Step 1: Filter the DataFrame
            filtered_df = df[(df["peptide_q"] <= 0.01) & (df["protein_q"] <= 0.01) & (df["label"] == 1)]
            st.metric("Number of PSMs with q-values <= 0.01", filtered_df.shape[0])
            
            cols = st.columns(2)
            with cols[0]:
                
                fig = view_identifications(filtered_df)
                show_fig(fig, "Identifications")
            
            with cols[1]:
                # log 2 transform ms2_intensity
                filtered_df['log2_ms2_intensity'] = np.log2(filtered_df['ms2_intensity'])
                fig = view_quantification(filtered_df)
                show_fig(fig, "Quantification")
            
            selected_mzml_file = st.selectbox("Select file", filtered_df['filename'].unique())
            
            single_file_df = filtered_df[filtered_df['filename'] == selected_mzml_file]
            # Arrange by peptide
            single_file_df = single_file_df.sort_values(by='proteins')
            
            c1, c2 = st.columns(2)
            c1.metric("Number of PSMs with q-values <= 0.01:", single_file_df.shape[0])
            rows = c1.dataframe(single_file_df, selection_mode="single-row", on_select="rerun")[
                "selection"
            ]["rows"]
            
            if rows:
                selected_row = single_file_df.iloc[rows, ]
                selected_mzml_file = str(Path(st.session_state.workspace, "sage-workflow/input-files/mzML-files", selected_row['filename'].values[0]))
                
                od_exp, meta_data = load_ms_file(selected_mzml_file)
                
                spectrum = od_exp.getSpectrumByNativeId(selected_row['scannr'].values[0])
                                
                spec_df = msspectrum_get_df(spectrum)
                spec_df['protein'] = selected_row['proteins'].values[0]
                spec_df['peptide'] = selected_row['peptide'].values[0]
                
                # get theoretical spectrum
                spec_theo = get_theo_spectrum(selected_row['peptide'].values[0])
                spec_theo_df = msspectrum_get_df(spec_theo)
                
                spec_alignment = SpectrumAlignment(spectrum, spec_theo)
                match_peaks_observed, match_peaks_theoretical = list(zip(*spec_alignment.alignment))
                
                obs_theo_match_df = spec_alignment.inspect()
                obs_theo_match_df['observed m/z'] = obs_theo_match_df['observed m/z'].astype(float)
                
                # Merge the DataFrames on the observed m/z and mz columns
                merged_df = spec_df.merge(obs_theo_match_df[['observed m/z', 'ion']], 
                                            left_on='mz', 
                                            right_on='observed m/z', 
                                            how='left')

                # Fill ion_annotation with the ion values from df1 where there is a match
                spec_df['ion_annotation'] = merged_df['ion']

                spec_df['peak_color'] = np.where(spec_df.index.isin(match_peaks_theoretical), 'black', 'grey')
                
                # spec_theo_df['peak_color'] = np.where(spec_theo_df.index.isin(match_peaks_observed), 'black', 'grey')
                
                # st.write()
                fig = spec_df.plot(x='mz', 
                                   y='intensity', 
                                   title=f"{selected_row['proteins'].values[0]} | {selected_row['peptide'].values[0]} <br><sup>{selected_row['scannr'].values[0]} @ RT = {selected_row['rt'].values[0]}</sup>", 
                                   kind='spectrum', 
                                   ion_annotation="ion_annotation",
                                   peak_color="peak_color",
                                   reference_spectrum=spec_theo_df, mirror_spectrum=True,
                                   backend='ms_plotly', 
                                   grid=False, show_plot=False)
                
                
                
                with c2:
                    st.metric("Number of matched peaks: ", str(len(spec_alignment.alignment)))
                    st.write(obs_theo_match_df)
                
                
                
                show_fig(fig.fig, f"mirror_spectrum_{selected_row['proteins'].values[0]}_{selected_row['peptide'].values[0]}_-_{selected_row['scannr'].values[0]}_at_RT = {selected_row['rt'].values[0]}")              
                
                # Get protein sequence from entries based on protein ID for current result (selected_row['proteins'].values[0])
                protein_sequence = [entry.sequence for entry in st.session_state["fasta_database"] if entry.identifier == selected_row['proteins'].values[0]]
                
                # st.write(protein_sequence)
                # st.write(selected_row['peptide'].values[0])
                
                # get peptide sequence from filtered_df for current protein
                filtered_df_peptides = filtered_df[filtered_df['proteins'] == selected_row['proteins'].values[0]][['peptide']].values
                filtered_df_peptides = list(np.unique(filtered_df_peptides))
                
                # st.write(filtered_df_peptides)
                
                highlighted_protein = highlight_peptides(protein_sequence[0], filtered_df_peptides, selected_row['peptide'].values[0])

                st.markdown(f"<p style='font-family:monospace;'>{highlighted_protein}</p>", unsafe_allow_html=True)
                
            else:
                st.info(
                    "💡 Select one ore more rows in the table to show the spectrum plot."
                )
            

        file = Path(
            self.workflow_dir, "results", "results.sage.tsv"
        )
        if file.exists():
            show_consensus_features()
        else:
            st.warning("No consensus feature file found. Please run workflow first.")