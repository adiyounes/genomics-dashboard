import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from modules.annotation.annotator import calculate_clinvar_score, calculate_pharmgkb_score, classify_significance


def test_calculate_clinvar_score_pathogenic():
    assert calculate_clinvar_score("Pathogenic", "homozygous_alt") == 1.0
    assert calculate_clinvar_score("Pathogenic", "heterozygous") == 0.7

def test_calculate_clinvar_score_vus():
    assert calculate_clinvar_score("uncertain significance", "homozygous_alt") == 0.5
    assert calculate_clinvar_score("uncertain significance", "heterozygous") == 0.35

def test_calculate_pharmgkb_score ():
    assert calculate_pharmgkb_score("1A", "homozygous_alt") == 0.9

def test_classify_significance():
    assert classify_significance("") == "unknown"
    assert classify_significance("assac pathogenic asda") == "pathogenic"
    assert classify_significance("as likely pathogenic qw") == "likely_pathogenic"
    assert classify_significance("uncertain as") == "vus"
    assert classify_significance("unc  vus ertain as") == "vus"
    assert classify_significance("uncert benign ain as") == "benign"