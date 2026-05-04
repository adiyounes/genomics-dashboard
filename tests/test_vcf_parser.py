import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from modules.ingestion.vcf_parser import detect_assembly, decode_genotype, validate_variant, parse_info_field

def test_detect_assembly_grch38():
    header_lines = [
    "##fileformat=VCFv4.2",
    "##reference=GRCh38",
    "##contig=<ID=chr1,length=248956422>"
    ]
    result = detect_assembly(header_lines)
    assert result == 'GRCh38'

def test_detect_assembly_rejected():
    header_lines = [
    "##fileformat=VCFv4.2",
    "##reference=GRCh37",
    "##contig=<ID=chr1,length=248956422>"
    ]
    result = detect_assembly(header_lines)
    assert result is None

def test_decode_genotype():
    assert decode_genotype("0/0") == "homozygous_ref"
    assert decode_genotype("0/1") == "heterozygous"
    assert decode_genotype("1/0") == "heterozygous"
    assert decode_genotype("1/1") == "homozygous_alt"
    assert decode_genotype("0|1") == "heterozygous"
    assert decode_genotype("?") == "unknown"

def test_validate_variant():
    assert validate_variant("chr17",123595,"A","T") == True
    assert validate_variant("crr17",123595,"A","T") == False
    assert validate_variant("chr17",-123595,"A","T") == False
    assert validate_variant("chr17",123595,"B","T")== False
    assert validate_variant("chr17",123595,"A","K")== False

def test_parse_info_field():
    assert parse_info_field("DP=30;AF=0.5") == {"DP": "30", "AF": "0.5"}
    assert parse_info_field("SOMATIC") == {}