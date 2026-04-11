# =============================================================================
# COLAB SETUP CELL
# Paste this at the top of each notebook when running in Google Colab.
# Skip this cell if running locally or on a cloud VM.
# =============================================================================

import os
import sys

IN_COLAB = "google.colab" in sys.modules

if IN_COLAB:
    # ── Mount Google Drive (recommended for persistent storage) ───────────────
    # Your processed .h5ad files will survive session disconnections
    from google.colab import drive
    drive.mount("/content/drive")

    # Set this to your preferred folder in Google Drive
    PROJECT_DIR = "/content/drive/MyDrive/ipf-scrna-seq-integration"
    DATA_DIR    = f"{PROJECT_DIR}/data"
    RESULTS_DIR = f"{PROJECT_DIR}/results"

    os.makedirs(f"{DATA_DIR}/raw", exist_ok=True)
    os.makedirs(f"{DATA_DIR}/processed", exist_ok=True)
    os.makedirs(f"{RESULTS_DIR}/figures", exist_ok=True)

    # ── Clone the repo (first time only) ─────────────────────────────────────
    if not os.path.exists(f"{PROJECT_DIR}/.git"):
        os.system(
            f"git clone https://github.com/YOUR_USERNAME/ipf-scrna-seq-integration.git "
            f"{PROJECT_DIR}"
        )
    else:
        os.system(f"cd {PROJECT_DIR} && git pull")

    # Add src/ to path so local modules are importable
    sys.path.insert(0, f"{PROJECT_DIR}/src")

    # ── Install dependencies ──────────────────────────────────────────────────
    # Colab has some packages pre-installed; this installs the rest
    # Only needs to run once per session (not once per cell)
    os.system("pip install -q scanpy anndata harmonypy scvi-tools liana-py pyscenic sccoda gseapy GEOparse")

    print(f"✓ Colab setup complete.")
    print(f"  Project dir: {PROJECT_DIR}")
    print(f"  Data dir:    {DATA_DIR}")

else:
    # ── Local / cloud VM paths ────────────────────────────────────────────────
    # Adjust PROJECT_DIR to wherever you cloned the repo
    PROJECT_DIR = os.path.dirname(os.path.abspath("__file__"))
    DATA_DIR    = f"{PROJECT_DIR}/data"
    RESULTS_DIR = f"{PROJECT_DIR}/results"

    sys.path.insert(0, f"{PROJECT_DIR}/src")

    print(f"✓ Local setup complete.")
    print(f"  Project dir: {PROJECT_DIR}")
    print(f"  Data dir:    {DATA_DIR}")


# =============================================================================
# SUBSAMPLE FLAG — set to True during development to work on a small slice
# =============================================================================
SUBSAMPLE      = False   # set True to use a fast 20k-cell subset
SUBSAMPLE_N    = 20_000  # cells to keep when subsampling

print(f"  Subsample mode: {'ON (' + str(SUBSAMPLE_N) + ' cells)' if SUBSAMPLE else 'OFF (full dataset)'}")
