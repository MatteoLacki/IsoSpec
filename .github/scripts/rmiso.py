import shutil
from pathlib import Path

path = Path(__file__).parent.parent / "IsoSpecPy"
print(f"Removing IsoSpecPy directory at: {path}")
shutil.rmtree(path, ignore_errors=True)
print("Done.")