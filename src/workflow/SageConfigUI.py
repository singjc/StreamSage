import streamlit as st
import json

class SageConfigUI:
    def __init__(self):
        # st.session_state["sage_config"] = config
        pass

    def display(self):
        st.title("Configuration Settings")

        # Database Configurations
        st.header("Database Configuration")
        st.session_state["sage_config"]['database']['bucket_size'] = st.number_input(
            "Bucket Size", value=st.session_state["sage_config"]['database']['bucket_size'], step=1
        )
        with st.expander("Enzyme: Expand for more options"):
            cols = st.columns(3)
            enzyme = st.session_state["sage_config"]['database']['enzyme']
            enzyme['missed_cleavages'] = cols[0].number_input(
                "Missed Cleavages", value=enzyme['missed_cleavages'], step=1
            )
            enzyme['min_len'] = cols[1].number_input(
                "Minimum Length", value=enzyme['min_len'], step=1
            )
            enzyme['max_len'] = cols[2].number_input(
                "Maximum Length", value=enzyme['max_len'], step=1
            )
            cols = st.columns(4)
            enzyme['cleave_at'] = cols[0].text_input(
                "Cleave At", value=enzyme['cleave_at']
            )
            enzyme['restrict'] = cols[1].text_input(
                "Restrict", value=enzyme['restrict']
            )
            enzyme['c_terminal'] = cols[2].checkbox(
                "C Terminal", value=enzyme['c_terminal']
            )
            enzyme['semi_enzymatic'] = cols[3].checkbox(
                "Semi Enzymatic", value=enzyme['semi_enzymatic']
            )

        cols = st.columns(3)
        st.session_state["sage_config"]['database']['fragment_min_mz'] = cols[0].number_input(
            "Fragment Min m/z", value=st.session_state["sage_config"]['database']['fragment_min_mz']
        )
        st.session_state["sage_config"]['database']['fragment_max_mz'] = cols[1].number_input(
            "Fragment Max m/z", value=st.session_state["sage_config"]['database']['fragment_max_mz']
        )
        st.session_state["sage_config"]['database']['peptide_min_mass'] = cols[2].number_input(
            "Peptide Min Mass", value=st.session_state["sage_config"]['database']['peptide_min_mass']
        )
        cols = st.columns(3)
        st.session_state["sage_config"]['database']['peptide_max_mass'] = cols[0].number_input(
            "Peptide Max Mass", value=st.session_state["sage_config"]['database']['peptide_max_mass']
        )
        st.session_state["sage_config"]['database']['ion_kinds'] = cols[1].multiselect(
            "Ion Kinds", options=["b", "y"], default=st.session_state["sage_config"]['database']['ion_kinds']
        )
        st.session_state["sage_config"]['database']['min_ion_index'] = cols[2].number_input(
            "Min Ion Index", value=st.session_state["sage_config"]['database']['min_ion_index'], step=1
        )

        cols = st.columns(2)
         # Static Mods
        cols[0].subheader("Static Modifications")
        static_mods_input = cols[0].text_area(
            "Static Mods (format: key\:value, one per line)",
            value="\n".join([f"{key}:{value}" for key, value in st.session_state["sage_config"]['database']['static_mods'].items()])
        )
        static_mods = {}
        for line in static_mods_input.splitlines():
            if ':' in line:
                key, value = line.split(':')
                static_mods[key.strip()] = float(value.strip())
        st.session_state["sage_config"]['database']['static_mods'] = static_mods

        # Variable Mods
        cols[1].subheader("Variable Modifications")
        variable_mods_input = cols[1].text_area(
            "Variable Mods (format: key\:value1,value2,..., one per line)",
            value="\n".join([f"{key}:{','.join(map(str, values))}" for key, values in st.session_state["sage_config"]['database']['variable_mods'].items()])
        )
        variable_mods = {}
        for line in variable_mods_input.splitlines():
            if ':' in line:
                key, values = line.split(':')
                variable_mods[key.strip()] = [float(v.strip()) for v in values.split(',')]
        st.session_state["sage_config"]['database']['variable_mods'] = variable_mods

        st.session_state["sage_config"]['database']['max_variable_mods'] = st.number_input(
            "Max Variable Mods", value=st.session_state["sage_config"]['database']['max_variable_mods'], step=1
        )
        
        cols = st.columns(2)
        st.session_state["sage_config"]['database']['decoy_tag'] = cols[0].text_input(
            "Decoy Tag", value=st.session_state["sage_config"]['database']['decoy_tag']
        )
        st.session_state["sage_config"]['database']['generate_decoys'] = cols[1].checkbox(
            "Generate Decoys", value=st.session_state["sage_config"]['database']['generate_decoys']
        )
        st.session_state["sage_config"]['database']['fasta'] = st.text_input(
            "FASTA Path", value=st.session_state["sage_config"]['database']['fasta']
        )

        # Quant Configurations
        st.header("Quant Configuration")
        cols = st.columns(2)    
        tmt_on = cols[0].checkbox("TMT", key="tmt_on")
        lfq_on = cols[1].checkbox("LFQ", value=True, key="lfq_on")
        
        if tmt_on:
            with st.expander("TMT: Expand for more options"):
                cols = st.columns(3)
                st.session_state["sage_config"]['quant']['tmt'] = cols[0].selectbox(
                    "TMT", options=["Tmt6", "Tmt10", "Tmt11", "Tmt16", "Tmt18"], index=["Tmt6", "Tmt10", "Tmt11", "Tmt16", "Tmt18"].index(st.session_state["sage_config"]['quant']['tmt'])
                )
                tmt_settings = st.session_state["sage_config"]['quant']['tmt_settings']
                tmt_settings['level'] = cols[1].number_input(
                    "TMT Level", value=tmt_settings['level'], step=1
                )
                tmt_settings['sn'] = cols[2].checkbox(
                    "TMT SN", value=tmt_settings['sn']
                )
        else:
            st.session_state["sage_config"]['quant']['tmt'] = None
        
        if lfq_on:
            with st.expander("LFQ: Expand for more options"):
                cols = st.columns(3)
                st.session_state["sage_config"]['quant']['lfq'] = cols[0].checkbox(
                    "LFQ", value=st.session_state["sage_config"]['quant']['lfq']
                )
                lfq_settings = st.session_state["sage_config"]['quant']['lfq_settings']
                lfq_settings['peak_scoring'] = cols[1].text_input(
                    "LFQ Peak Scoring", value=lfq_settings['peak_scoring']
                )
                lfq_settings['integration'] = cols[2].selectbox(
                    "LFQ Integration", options=["Sum", "Apex"], index=["Sum", "Apex"].index(lfq_settings['integration'])
                )
                cols = st.columns(3)
                lfq_settings['spectral_angle'] = cols[0].number_input(
                    "LFQ Spectral Angle", value=lfq_settings['spectral_angle']
                )
                lfq_settings['ppm_tolerance'] = cols[1].number_input(
                    "LFQ PPM Tolerance", value=lfq_settings['ppm_tolerance']
                )
                lfq_settings['combine_charge_states'] = cols[2].checkbox(
                    "LFQ Combine Charge States", value=lfq_settings['combine_charge_states']
                )
        else:
            st.session_state["sage_config"]['quant']['lfq'] = None

        # Tolerances
        st.header("Tolerances")
        cols = st.columns(2)
        st.session_state["sage_config"]['precursor_tol']['ppm'] = cols[0].text_input(
            "Precursor Tolerance (ppm)", value=", ".join(map(str, st.session_state["sage_config"]['precursor_tol']['ppm']))
        )
        st.session_state["sage_config"]['precursor_tol']['ppm'] = [int(x.strip()) for x in st.session_state["sage_config"]['precursor_tol']['ppm'].split(',')]

        st.session_state["sage_config"]['fragment_tol']['ppm'] = cols[1].text_input(
            "Fragment Tolerance (ppm)", value=", ".join(map(str, st.session_state["sage_config"]['fragment_tol']['ppm']))
        )
        st.session_state["sage_config"]['fragment_tol']['ppm'] = [int(x.strip()) for x in st.session_state["sage_config"]['fragment_tol']['ppm'].split(',')]

        # Other Configurations
        st.header("Other Configurations")
        cols = st.columns(4)
        st.session_state["sage_config"]['precursor_charge'] = cols[0].text_input(
            "Precursor Charge", value=", ".join(map(str, st.session_state["sage_config"]['precursor_charge']))
        )
        st.session_state["sage_config"]['precursor_charge'] = [int(x.strip()) for x in st.session_state["sage_config"]['precursor_charge'].split(',')]

        st.session_state["sage_config"]['isotope_errors'] = cols[1].text_input(
            "Isotope Errors", value=", ".join(map(str, st.session_state["sage_config"]['isotope_errors']))
        )
        st.session_state["sage_config"]['isotope_errors'] = [int(x.strip()) for x in st.session_state["sage_config"]['isotope_errors'].split(',')]

        st.session_state["sage_config"]['deisotope'] = cols[2].checkbox(
            "Deisotope", value=st.session_state["sage_config"]['deisotope']
        )
        st.session_state["sage_config"]['chimera'] = cols[3].checkbox(
            "Chimera", value=st.session_state["sage_config"]['chimera']
        )
        cols = st.columns(4)
        st.session_state["sage_config"]['wide_window'] = cols[0].checkbox(
            "Wide Window", value=st.session_state["sage_config"]['wide_window']
        )
        st.session_state["sage_config"]['predict_rt'] = cols[1].checkbox(
            "Predict RT", value=st.session_state["sage_config"]['predict_rt']
        )
        st.session_state["sage_config"]['min_peaks'] = cols[2].number_input(
            "Min Peaks", value=st.session_state["sage_config"]['min_peaks'], step=1
        )
        st.session_state["sage_config"]['max_peaks'] = cols[3].number_input(
            "Max Peaks", value=st.session_state["sage_config"]['max_peaks'], step=1
        )
        cols = st.columns(4)
        st.session_state["sage_config"]['min_matched_peaks'] = cols[0].number_input(
            "Min Matched Peaks", value=st.session_state["sage_config"]['min_matched_peaks'], step=1
        )
        st.session_state["sage_config"]['max_fragment_charge'] = cols[1].number_input(
            "Max Fragment Charge", value=st.session_state["sage_config"]['max_fragment_charge'], step=1
        )
        st.session_state["sage_config"]['report_psms'] = cols[2].number_input(
            "Report PSMs", value=st.session_state["sage_config"]['report_psms'], step=1
        )
        st.session_state["sage_config"]['output_directory'] = cols[3].text_input(
            "Output Directory", value=st.session_state["sage_config"]['output_directory']
        )
        st.session_state["sage_config"]['mzml_paths'] = st.text_area(
            "mzML Paths", value="\n".join(st.session_state["sage_config"]['mzml_paths'])
        )
        st.session_state["sage_config"]['mzml_paths'] = st.session_state["sage_config"]['mzml_paths'].splitlines()

        # # Save Configuration
        # if st.button("Save Configuration"):
        #     with open("updated_config.json", "w") as json_file:
        #         json.dump(st.session_state["sage_config"], json_file, indent=2)
        #     st.success("Configuration saved!")