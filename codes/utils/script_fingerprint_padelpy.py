"""
PaDEL Descriptor Calculation 
"""

import os
import glob
import pandas as pd
from padelpy import padeldescriptor
import logging

# Configuration
FP_ALL = [
    'AtomPairs2DCount', 
    'AtomPairs2D', 
    'EState', 
    'CDKextended', 
    'CDK', 
    'CDKgraphonly',
    'KlekotaRothCount', 
    'KlekotaRoth', 
    'MACCS', 
    'PubChem', 
    'SubstructureCount', 
    'Substructure'
]
FP_LIST = [
    'MACCS', 'PubChem', 'Substructure'
]  # Selected fingerprint types

TMP_SMI_DIR = "path/to/tmp/tmp_{}.smi"
FP_OUTPUT_DIR = "path/to/fp/{}"
XML_FILES_DIR = "path/to/fp/fingerprints_xml/*.xml"
FP_META_PATH = "path/to/fp/fp_meta.csv"

# Logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def setup_directories(cid):
    """Create necessary directories for storing temporary and output files."""
    output_dir = FP_OUTPUT_DIR.format(cid)
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def write_smiles_to_file(cid, smile):
    """Write the SMILES string and CID to a temporary file."""
    tmp_smi_path = TMP_SMI_DIR.format(cid)
    with open(tmp_smi_path, "w") as f:
        f.write(f"{smile}\t{cid}")
    return tmp_smi_path

def calculate_fingerprints(cid, tmp_smi_path, output_dir, fp_list, fp_all, xml_files_dir, fp_meta_path):
    """Calculate fingerprints using PaDEL-Descriptor for the given SMILES."""
    xml_files = glob.glob(xml_files_dir)
    xml_files.sort()
    fp_xml_dict = dict(zip(fp_all, xml_files))
    fp_meta = pd.read_csv(fp_meta_path)

    for fp in fp_list:
        output_file = os.path.join(output_dir, f"{fp}.csv")
        try:
            padeldescriptor(
                mol_dir=tmp_smi_path,
                d_file=output_file,
                descriptortypes=fp_xml_dict[fp],
                detectaromaticity=True,
                standardizenitro=True,
                standardizetautomers=True,
                removesalt=True,
                log=False,
                fingerprints=True
            )
            logger.info(f"Successfully calculated {fp} fingerprints for CID {cid}")
        except Exception as e:
            logger.error(f"Failed to calculate {fp} fingerprints for CID {cid}: {e}")
            fp_meta_sub = fp_meta[fp_meta.Type.isin([fp])]
            output_df = pd.DataFrame([cid] + [None] * fp_meta_sub.shape[0]).transpose()
            output_df.columns = ["Name"] + list(fp_meta_sub.Name)
            output_df.to_csv(output_file, index=False)
            logger.info(f"Created placeholder file for {fp} fingerprints")

def merge_fingerprints(cid, output_dir, fp_list):
    """Merge all calculated fingerprints into a single DataFrame."""
    fp_feat = pd.DataFrame({"Name": [int(cid)]})
    for fp in fp_list:
        fp_file = os.path.join(output_dir, f"{fp}.csv")
        if os.path.exists(fp_file):
            feat_tmp = pd.read_csv(fp_file)
            fp_feat = pd.merge(fp_feat, feat_tmp, on="Name", how="left")
        else:
            logger.warning(f"Fingerprint file {fp_file} not found")
    return fp_feat

def get_padel_fp(cid, smile):
    """Main function to calculate and merge fingerprints for a given CID and SMILES."""
    logger.info(f"Processing compound CID: {cid}, SMILES: {smile}")

    # Setup directories and write SMILES to file
    output_dir = setup_directories(cid)
    tmp_smi_path = write_smiles_to_file(cid, smile)

    # Calculate fingerprints
    calculate_fingerprints(cid, tmp_smi_path, output_dir, FP_LIST, FP_ALL, XML_FILES_DIR, FP_META_PATH)

    # Merge fingerprints
    fp_feat = merge_fingerprints(cid, output_dir, FP_LIST)

    # Clean up temporary file
    os.remove(tmp_smi_path)
    logger.info(f"Temporary file {tmp_smi_path} removed")

    return fp_feat

# Example usage
if __name__ == "__main__":
    cid = "86287635"
    smile = "C[C@@H]1CC2=C([C@H](N1CC(C)(C)F)C3=C(C=C(C=C3F)/C=C/C(=O)O)F)NC4=CC=CC=C24"
    result = get_padel_fp(cid, smile)
    print(result)