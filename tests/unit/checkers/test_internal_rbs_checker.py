import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


@pytest.fixture
def checker():
    checker = InternalRBSChecker()
    checker.initiate()
    return checker


def test_internal_rbs_detected(checker):
    """
    Sequences that contain a Shine-Dalgarno-like sequence followed by ATG
    within the correct spacer distance. All should return False.
    """
    bad_seqs = [
        # Canonical AGGAGG + 6bp spacer + ATG
        "ATGCGTAAAGGAGGTTTTTTTATGCCCGTA",
        # AGGAG + 5bp spacer + ATG
        "CCCAGGAGCCCCCATGGGGTTTAAA",
        # GAGG + 7bp spacer + ATG
        "TTTGAGGTTTTTTTTATGAAACCC",
        # AAGG + 5bp spacer + ATG
        "GGGAAGGAAAAATGCGTTTACCC",
        # AGGA + 6bp spacer + ATG
        "CCCAGGATTTTTTTATGGGGTTT",
        # GGAG + 5bp spacer + ATG
        "AAAGGAGCCCCCATGTTTTGGG",
        # Canonical SD near start of sequence
        "AGGAGGTAAAAATGCCCAAATTT",
        # SD on reverse complement strand
        "CATTTTATTTTTCCCTCCTTTGCG",  # RC contains AGGAGG near ATG
    ]

    print(">> Testing sequences with internal RBS (expected False)")
    for seq in bad_seqs:
        result, site = checker.run(seq.upper())
        print(f"result: {result}, site: {site} on {seq}")
        assert result == False, f"Expected False for sequence with internal RBS: {seq}"


def test_no_internal_rbs(checker):
    """
    Sequences with no Shine-Dalgarno pattern close to an ATG.
    All should return True.
    """
    good_seqs = [
        # No SD pattern anywhere
        "ATGCGTTTACCCGGTCTGCCCAAATTT",
        # Has ATG but no upstream SD
        "CCCTTTTTTCCCGGTCTGCCCATGAAA",
        # Has SD-like but ATG is too far away (>10bp)
        "AGGAGGTTTTTTTTTTTTTTATGCCC",
        # Has SD-like but ATG is too close (<4bp spacer)
        "AGGAGGATGCCC",
        # Random sequence with no meaningful patterns
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        # All same nucleotide (no SD possible)
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        # Typical codon-optimized stretch with no SD
        "ATGCTGGTACAGCCGACCGTTGAAGGT",
    ]

    print("\n>> Testing sequences without internal RBS (expected True)")
    for seq in good_seqs:
        result, site = checker.run(seq.upper())
        print(f"result: {result}, site: {site} on {seq}")
        assert result == True, f"Expected True for clean sequence: {seq}"


def test_invalid_input_raises(checker):
    """
    Invalid inputs should raise a ValueError.
    """
    with pytest.raises(ValueError):
        checker.run("")

    with pytest.raises(ValueError):
        checker.run(None)