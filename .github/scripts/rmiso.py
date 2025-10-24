import shutil
from pathlib import Path

path = Path(__file__).parent.parent.parent / "IsoSpecPy"
print(f"Removing IsoSpecPy directory at: {path}")
print("Its contents are:", list(path.iterdir()))
shutil.rmtree(path, ignore_errors=True)
print("Done.")