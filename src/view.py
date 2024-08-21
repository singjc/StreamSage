import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import pyopenms as poms
from .plotting.MSExperimentPlotter import plotMSExperiment
from .common import show_fig

from typing import Union


def get_df(file: Union[str, Path]) -> pd.DataFrame:
    """
    Load a Mass Spectrometry (MS) experiment from a given mzML file and return
    a pandas dataframe representation of the experiment.

    Args:
        file (Union[str, Path]): The path to the mzML file to load.

    Returns:
        pd.DataFrame: A pandas DataFrame with the following columns: "mslevel",
        "precursormz", "mzarray", and "intarray". The "mzarray" and "intarray"
        columns contain NumPy arrays with the m/z and intensity values for each
        spectrum in the mzML file, respectively.
    """
    exp = poms.MSExperiment()
    poms.MzMLFile().load(str(file), exp)
    df_spectra = exp.get_df()
    df_spectra["MS level"] = [spec.getMSLevel() for spec in exp]
    precs = []
    for spec in exp:
        p = spec.getPrecursors()
        if p:
            precs.append(p[0].getMZ())
        else:
            precs.append(np.nan)
    df_spectra["precursor m/z"] = precs
    df_spectra["max intensity m/z"] = df_spectra.apply(
        lambda x: x["mzarray"][x["intarray"].argmax()], axis=1
    )
    if not df_spectra.empty:
        st.session_state["view_spectra"] = df_spectra
    else:
        st.session_state["view_spectra"] = pd.DataFrame()
    exp_ms2 = poms.MSExperiment()
    exp_ms1 = poms.MSExperiment()
    for spec in exp:
        if spec.getMSLevel() == 1:
            exp_ms1.addSpectrum(spec)
        elif spec.getMSLevel() == 2:
            exp_ms2.addSpectrum(spec)
    if not exp_ms1.empty():
        st.session_state["view_ms1"] = exp_ms1.get_df(long=True)
    else:
        st.session_state["view_ms1"] = pd.DataFrame()
    if not exp_ms2.empty():
        st.session_state["view_ms2"] = exp_ms2.get_df(long=True)
    else:
        st.session_state["view_ms2"] = pd.DataFrame()

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

def plot_bpc_tic() -> go.Figure:
    """Plot the base peak and total ion chromatogram (TIC).

    Returns:
        A plotly Figure object containing the BPC and TIC plot.
    """
    fig = go.Figure()
    if st.session_state.view_tic:
        df = st.session_state.view_ms1.groupby("RT").sum().reset_index()
        fig.add_scatter(
            x=df["RT"],
            y=df["inty"],
            mode="lines",
            line=dict(color="#f24c5c", width=3),  # OpenMS red
            name="TIC",
            showlegend=True,
        )
    if st.session_state.view_bpc:
        df = st.session_state.view_ms1.groupby("RT").max().reset_index()
        fig.add_scatter(
            x=df["RT"],
            y=df["inty"],
            mode="lines",
            line=dict(color="#2d3a9d", width=3),  # OpenMS blue
            name="BPC",
            showlegend=True,
        )
    if st.session_state.view_eic:
        df = st.session_state.view_ms1
        target_value = st.session_state.view_eic_mz.strip().replace(",", ".")
        try:
            target_value = float(target_value)
            ppm_tolerance = st.session_state.view_eic_ppm
            tolerance = (target_value * ppm_tolerance) / 1e6

            # Filter the DataFrame
            df_eic = df[(df['mz'] >= target_value - tolerance) & (df['mz'] <= target_value + tolerance)]
            if not df_eic.empty:
                fig.add_scatter(
                    x=df_eic["RT"],
                    y=df_eic["inty"],
                    mode="lines",
                    line=dict(color="#f6bf26", width=3),
                    name="XIC",
                    showlegend=True,
                )
        except ValueError:
            st.error("Invalid m/z value for XIC provided. Please enter a valid number.")

    fig.update_layout(
        title=f"{st.session_state.view_selected_file}",
        xaxis_title="retention time (s)",
        yaxis_title="intensity",
        plot_bgcolor="rgb(255,255,255)",
        height=500,
    )
    fig.layout.template = "plotly_white"
    return fig


@st.cache_resource
def plot_ms_spectrum(spec, title, color):
    """
    Takes a pandas Series (spec) and generates a needle plot with m/z and intensity dimension.

    Args:
        spec: Pandas Series representing the mass spectrum with "mzarray" and "intarray" columns.
        title: Title of the plot.
        color: Color of the line in the plot.

    Returns:
        A Plotly Figure object representing the needle plot of the mass spectrum.
    """

    # Every Peak is represented by three dots in the line plot: (x, 0), (x, y), (x, 0)
    def create_spectra(x, y, zero=0):
        x = np.repeat(x, 3)
        y = np.repeat(y, 3)
        y[::3] = y[2::3] = zero
        return pd.DataFrame({"mz": x, "intensity": y})

    df = create_spectra(spec["mzarray"], spec["intarray"])
    fig = px.line(df, x="mz", y="intensity")
    fig.update_traces(line_color=color)
    fig.add_hline(0, line=dict(color="#DDDDDD"), line_width=3)
    fig.update_layout(
        showlegend=False,
        title_text=title,
        xaxis_title="m/z",
        yaxis_title="intensity",
        plot_bgcolor="rgb(255,255,255)",
        dragmode="select",
    )
    # add annotations
    top_indices = np.argsort(spec["intarray"])[-5:][::-1]
    for index in top_indices:
        mz = spec["mzarray"][index]
        i = spec["intarray"][index]
        fig.add_annotation(
            dict(
                x=mz,
                y=i,
                text=str(round(mz, 5)),
                showarrow=False,
                xanchor="left",
                font=dict(
                    family="Open Sans Mono, monospace",
                    size=12,
                    color=color,
                ),
            )
        )
    fig.layout.template = "plotly_white"
    # adjust x-axis limits to not cut peaks and annotations
    x_values = [trace.x for trace in fig.data]
    xmin = min([min(values) for values in x_values])
    xmax = max([max(values) for values in x_values])
    padding = 0.15 * (xmax - xmin)
    fig.update_layout(
        xaxis_range=[
            xmin - padding,
            xmax + padding,
        ]
    )
    return fig


@st.experimental_fragment
def view_peak_map():
    df = st.session_state.view_ms1
    if "view_peak_map_selection" in st.session_state:
        box = st.session_state.view_peak_map_selection.selection.box
        if box:
            df = st.session_state.view_ms1.copy()
            df = df[df["RT"] > box[0]["x"][0]]
            df = df[df["mz"] > box[0]["y"][1]]
            df = df[df["mz"] < box[0]["y"][0]]
            df = df[df["RT"] < box[0]["x"][1]]
    peak_map = plotMSExperiment(
        df, plot3D=False, title=st.session_state.view_selected_file
    )
    c1, c2 = st.columns(2)
    with c1:
        st.info(
            "💡 Zoom in via rectangular selection for more details and 3D plot. Double click plot to zoom back out."
        )
        show_fig(
            peak_map,
            f"peak_map_{st.session_state.view_selected_file}",
            selection_session_state_key="view_peak_map_selection",
        )
    with c2:
        if df.shape[0] < 2500:
            peak_map_3D = plotMSExperiment(df, plot3D=True, title="")
            st.pyplot(peak_map_3D, use_container_width=True)


@st.experimental_fragment
def view_spectrum():
    cols = st.columns([0.34, 0.66])
    with cols[0]:
        df = st.session_state.view_spectra.copy()
        df["spectrum ID"] = df.index + 1
        event = st.dataframe(
            df,
            column_order=[
                "spectrum ID",
                "RT",
                "MS level",
                "max intensity m/z",
                "precursor m/z",
            ],
            selection_mode="single-row",
            on_select="rerun",
            use_container_width=True,
            hide_index=True,
        )
        rows = event.selection.rows
    with cols[1]:
        if rows:
            df = st.session_state.view_spectra.iloc[rows[0]]
            if "view_spectrum_selection" in st.session_state:
                box = st.session_state.view_spectrum_selection.selection.box
                if box:
                    mz_min, mz_max = sorted(box[0]["x"])
                    mask = (df["mzarray"] > mz_min) & (df["mzarray"] < mz_max)
                    df["intarray"] = df["intarray"][mask]
                    df["mzarray"] = df["mzarray"][mask]

            if df["mzarray"].size > 0:
                title = f"{st.session_state.view_selected_file}  spec={rows[0]+1}  mslevel={df['MS level']}"
                if df["precursor m/z"] > 0:
                    title += f" precursor m/z: {round(df['precursor m/z'], 4)}"
                fig = plot_ms_spectrum(df, title, "#2d3a9d")
                show_fig(fig, title.replace(" ", "_"), True, "view_spectrum_selection")
            else:
                st.session_state.pop("view_spectrum_selection")
                st.rerun()
        else:
            st.info("💡 Select rows in the spectrum table to display plot.")


@st.experimental_fragment()
def view_bpc_tic():
    cols = st.columns(5)
    cols[0].checkbox(
        "Total Ion Chromatogram (TIC)", True, key="view_tic", help="Plot TIC."
    )
    cols[1].checkbox(
        "Base Peak Chromatogram (BPC)", True, key="view_bpc", help="Plot BPC."
    )
    cols[2].checkbox(
        "Extracted Ion Chromatogram (EIC/XIC)", True, key="view_eic", help="Plot extracted ion chromatogram with specified m/z."
    )
    cols[3].text_input(
        "XIC m/z",
        "235.1189",
        help="m/z for XIC calculation.",
        key="view_eic_mz",
    )
    cols[4].number_input(
        "XIC ppm tolerance",
        0.1, 50.0, 10.0, 1.0,
        help="Tolerance for XIC calculation (ppm).",
        key="view_eic_ppm"
    )
    fig = plot_bpc_tic()
    show_fig(fig, f"BPC-TIC-{st.session_state.view_selected_file}")

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