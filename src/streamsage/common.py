
import os
import sys
import re
from pathlib import Path
import requests
import zipfile
import tarfile
import streamlit as st
import pyopenms as poms

from src.common.common import OS_PLATFORM

def load_fasta():
    """
    Load the FASTA database file into the session state.
    """
    
    entries = []
    f = poms.FASTAFile()
    f.load(str(st.session_state["sage_config"]['database']['fasta']), entries)
    if "fasta_database" not in st.session_state:
        st.session_state["fasta_database"] = entries
        

def download_and_unpack_sage_exec(target_dir="./bin/", version="v0.14.7"):
    """
    Download the Sage executable, unpack it, and store it in a specific directory.
    Returns the full path of the unpacked Sage executable.
    """
    
    # Define the base URL for the GitHub release
    base_url = "https://github.com/lazear/sage/releases/download"

    # Determine the URL based on the platform
    if sys.platform == "win32":
        # NOTE: The v0.14.7 release for windows does not work.
        # url = f"{base_url}/{version}/sage-{version}-x86_64-pc-windows-msvc.zip"
        url = "https://github.com/singjc/StreamSage/releases/download/v0.0.1-alpha/sage-v0.14.7-x86_64-pc-windows-msvc.zip"
        file_extension = ".zip"
        executable_name = "sage.exe"
    elif sys.platform == "linux" or sys.platform == "linux2":
        url = f"{base_url}/{version}/sage-{version}-x86_64-unknown-linux-gnu.tar.gz"
        file_extension = ".tar.gz"
        executable_name = "sage"
    elif sys.platform == "darwin":
        url = f"{base_url}/{version}/sage-{version}-x86_64-apple-darwin.tar.gz"
        file_extension = ".tar.gz"
        executable_name = "sage"
    else:
        raise ValueError(f"Unsupported platform: {sys.platform}. Try manually downloading the Sage executable for your platform from the Sage GitHub release page.")

    # Download the file
    print(f"Downloading Sage executable from: {url}")
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code != 200:
        raise Exception(f"Failed to download file, status code: {response.status_code}")

    # Ensure the target directory exists
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Determine the full path to save the downloaded file
    file_name = os.path.join(target_dir, f"sage{file_extension}")
    
    # Save the file to the target directory
    with open(file_name, 'wb') as file:
        file.write(response.content)
    print(f"Downloaded file saved as: {file_name}")
    
    # Unpack the downloaded file
    if file_extension == ".zip":
        print("Unpacking ZIP file...")
        with zipfile.ZipFile(file_name, 'r') as zip_ref:
            zip_ref.extractall(target_dir)
    elif file_extension == ".tar.gz":
        print("Unpacking TAR.GZ file...")
        with tarfile.open(file_name, 'r:gz') as tar_ref:
            tar_ref.extractall(target_dir)
    
    print(f"Sage executable unpacked and stored in: {target_dir}")

    # Perform a recursive search for the executable in the target directory
    def find_executable_in_directory(directory, exec_name):
        """
        Recursively search for the executable in the specified directory.
        """
        for root, dirs, files in os.walk(directory):
            if exec_name in files:
                return os.path.join(root, exec_name)
        return None

    # Search for the Sage executable
    exec_path = find_executable_in_directory(target_dir, executable_name)

    if not exec_path:
        raise Exception(f"Executable not found in the unpacked directory: {target_dir}")

    print(f"Full path to the Sage executable: {exec_path}")
    st.success(f"Sage executable downloaded to: {exec_path}")
    
    # Return the full executable path
    return exec_path
        
        