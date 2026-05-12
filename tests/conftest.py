import sys
from pathlib import Path

# Make workflow/scripts/ importable without a package __init__
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))
