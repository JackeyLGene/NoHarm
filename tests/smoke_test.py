from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from noharm.scan import ScanParams, scan_fasta


def test_demo_scan_has_branch_gene():
    genes, isoforms, summary = scan_fasta(ROOT / "data" / "demo_isoforms.fa", ScanParams())
    assert summary["multi_isoform_genes"] == 3
    assert len(isoforms) == 6
    assert genes[0].gene == "GENE_BRANCH"
    assert genes[0].ch_range > 0.0


if __name__ == "__main__":
    test_demo_scan_has_branch_gene()
    print("smoke test passed")

