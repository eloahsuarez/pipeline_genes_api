from bioinformatics import find_conserved_regions, generate_primer_pairs, PrimerDesignParams


def test_conserved_region_exact():
    alignment = ">a\nAACCGGTTAA\n>b\nAACCGGTTAA\n>c\nAACCGGTTAA\n"
    regions = find_conserved_regions(alignment, 1.0, 1.0, 4)
    assert len(regions) == 1
    assert regions[0].sequence == "AACCGGTTAA"
    assert regions[0].length == 10


def test_gap_breaks_exact_region():
    alignment = ">a\nAACCGGTTAA\n>b\nAACC-GTTAA\n>c\nAACCGGTTAA\n"
    regions = find_conserved_regions(alignment, 1.0, 1.0, 4)
    assert [r.sequence for r in regions] == ["AACC", "GTTAA"]


def test_primer_generation_returns_pairs():
    seq = "ACGTCGATGCTAGCTACGATCGTACGATGCTAGCTACGATCGTACGATGCTAGCTACGATCGTACGATGCTAGCTACGATCGTACGATGCTAGCTACGATCGT"
    alignment = f">a\n{seq}\n>b\n{seq}\n>c\n{seq}\n"
    params = PrimerDesignParams(
        min_primer_len=18, max_primer_len=22, min_gc=25, max_gc=75,
        min_tm=40, max_tm=80, min_amplicon=60, max_amplicon=110,
        target_amplicon=80, top_pairs=5,
    )
    pairs = generate_primer_pairs(alignment, params)
    assert pairs
    assert pairs[0].forward.orientation == "F"
    assert pairs[0].reverse.orientation == "R"
