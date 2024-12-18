import re
import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import pyopenms as poms
# from .plotting.MSExperimentPlotter import plotMSExperiment
from src.common.common import show_fig

from typing import Union

def remove_modifications(peptide_sequences):
    modified_sequences = []
    for sequence in peptide_sequences:
        modified_sequence = re.sub(r'\[[^\]]*\]', '', sequence)
        if modified_sequence not in modified_sequences:
            modified_sequences.append(modified_sequence)
    return modified_sequences

def highlight_peptides(protein_seq, identified_peptides, target_peptide):
    highlighted_seq = protein_seq
    
    target_peptide = re.sub(r'\[[^\]]*\]', '', target_peptide)
    
    identified_peptides = remove_modifications(identified_peptides)
    
    # pop target peptide from identified peptides
    if target_peptide in identified_peptides:
        identified_peptides.remove(target_peptide)
    
    # Highlight the target peptide in yellow
    start = 0
    while True:
        start = protein_seq.find(target_peptide, start)
        if start == -1:
            break
        end = start + len(target_peptide)
        highlighted_seq = highlighted_seq[:start] + f"<span style='background-color:yellow'>{protein_seq[start:end]}</span>" + highlighted_seq[end:]
        start = end

    # # Highlight identified peptides in grey
    # if len(identified_peptides) > 0:
    #     for peptide in identified_peptides:
    #         start = 0
    #         while True:
    #             start = highlighted_seq.find(peptide, start)
    #             if start == -1:
    #                 break
    #             end = start + len(peptide)
    #             highlighted_seq = (
    #                 highlighted_seq[:start] +
    #                 f"<span style='background-color:lightgrey'>{highlighted_seq[start:end]}</span>" +
    #                 highlighted_seq[end:]
    #             )
    #             start = end

    return highlighted_seq

@st.cache_resource
def load_ms_file(file: Union[str, Path]):                
    od_exp = poms.OnDiscMSExperiment()
    od_exp.openFile(file)
    meta_data = od_exp.getMetaData()
    return od_exp, meta_data

def _add_meta_values(df: pd.DataFrame, object: any) -> pd.DataFrame:
    """
    Adds metavalues from given object to given DataFrame.
    
    Args:
        df (pd.DataFrame): DataFrame to which metavalues will be added.
        object (any): Object from which metavalues will be extracted.
    
    Returns:
        pd.DataFrame: DataFrame with added meta values.
    """
    mvs = []
    object.getKeys(mvs)
    for k in mvs:
        v = object.getMetaValue(k)
        dtype = 'U100'
        try:
            v = int(v)
            dtype = int
        except ValueError:
            try:
                v = float(v)
                dtype = 'double'
            except ValueError:
                dtype = f'U{len(v)}'
        
        df[k.decode()] = np.full(df.shape[0], v, dtype=np.dtype(dtype))

    return df

def msspectrum_get_df(spec, export_meta_values: bool = True) -> pd.DataFrame:
        """
        Returns a DataFrame representation of the MSSpectrum.

        Args:
            export_meta_values (bool): Whether to export meta values.

        Returns:
            pd.DataFrame: DataFrame representation of the MSSpectrum.
        """
        mzs, intensities = spec.get_peaks()

        df = pd.DataFrame({'mz': mzs, 'intensity': intensities})

        cnt = df.shape[0]
        
        # ion mobility
        df['ion_mobility'] = np.array([i for i in spec.getFloatDataArrays()[0]]) if spec.containsIMData() else np.nan
        df['ion_mobility_unit'] = np.full(cnt, spec.getDriftTimeUnitAsString(), dtype=np.dtype(f'U{len(spec.getDriftTimeUnitAsString())}'))

        df['ms_level'] = np.full(cnt, spec.getMSLevel(), dtype=np.dtype('uint16'))

        precs = spec.getPrecursors()
        df['precursor_mz'] = np.full(cnt, (precs[0].getMZ() if precs else 0.0), dtype=np.dtype('double'))
        df['precursor_charge'] = np.full(cnt, (precs[0].getCharge() if precs else 0), dtype=np.dtype('uint16'))
        
        df['native_id'] = np.full(cnt, spec.getNativeID(), dtype=np.dtype('U100'))

        # peptide sequence
        peps = spec.getPeptideIdentifications()  # type: list[PeptideIdentification]
        seq = ''
        if peps:
            hits = peps[0].getHits()
            if hits:
                seq = hits[0].getSequence().toString()
        df['sequence'] = np.full(cnt, seq, dtype=np.dtype(f'U{len(seq)}'))

        # ion annotations in string data array with names IonName or IonNames
        ion_annotations = np.full(cnt, '', dtype=np.dtype('U1'))
        for sda in spec.getStringDataArrays():
            if sda.getName() == 'IonNames':
                decoded = [ion.decode() for ion in sda]
                if len(decoded) == df.shape[0]:
                    ion_annotations = np.array(decoded, dtype=np.dtype(f'U{len(max(decoded))}'))
                    break
        df['ion_annotation'] = ion_annotations

        if export_meta_values:
            df = _add_meta_values(df, spec)

        return df

def get_theo_spectrum(peptide):
    tsg = poms.TheoreticalSpectrumGenerator()

    theo_spectrum = poms.MSSpectrum()

    p = tsg.getParameters()

    p.setValue("add_y_ions", "true")

    p.setValue("add_b_ions", "true")

    p.setValue("add_metainfo", "true")

    tsg.setParameters(p)

    peptide = poms.AASequence.fromString(peptide)

    tsg.getSpectrum(theo_spectrum, peptide, 1, 2)
    
    return theo_spectrum

class SpectrumAlignment:
    def __init__(self, observed_spectrum, theo_spectrum):
        self.observed_spectrum = observed_spectrum
        self.theo_spectrum = theo_spectrum
        self.alignment = self.spectrum_alignment(observed_spectrum, theo_spectrum)

    @staticmethod
    def spectrum_alignment(observed_spectrum, theo_spectrum):
        alignment = []

        spa = poms.SpectrumAlignment()

        p = spa.getParameters()

        # use 0.5 Da tolerance (Note: for high-resolution data we could also use ppm by setting the is_relative_tolerance value to true)

        p.setValue("tolerance", 0.5)

        p.setValue("is_relative_tolerance", "false")

        spa.setParameters(p)

        # align both spectra
        spa.getSpectrumAlignment(alignment, theo_spectrum, observed_spectrum)
        
        return alignment
    
    def inspect(self):
        t = []

        for theo_idx, obs_idx in self.alignment:
            ion_name = self.theo_spectrum.getStringDataArrays()[0][theo_idx].decode()
            ion_charge = self.theo_spectrum.getIntegerDataArrays()[0][theo_idx]

            t.append(
                [
                    ion_name,
                    str(ion_charge),
                    str(self.theo_spectrum[theo_idx].getMZ()),
                    str(self.observed_spectrum[obs_idx].getMZ()),
                ]
            )

        df = pd.DataFrame(t, columns=["ion", "charge", "theo. m/z", "observed m/z"])
        
        return df

        
def view_identifications(df, level="psm"):
    """
    Display identifications in a DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing identifications.
        level (str): Level of identifications (psm or peptide).

    Returns:
        None
    """
    # Group by filename and count number of rows
    counts = df.groupby("filename").size().reset_index(name="count")
    
    # Plot number of identifications per file
    fig = px.bar(counts, x="filename", y="count", title=f"Number of {level}s per file")
    return fig

def view_quantification(df, level="psm", value="log2_ms2_intensity"):
    """"
    Display boxplot of quantification values per file
    """
    
    # plot boxplot
    fig = px.box(df, x="filename", y=value, title=f"Quantification values per {level}")
    
    return fig