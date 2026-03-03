import sys
from pathlib import Path

# Make scripts/ importable without a package __init__
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
