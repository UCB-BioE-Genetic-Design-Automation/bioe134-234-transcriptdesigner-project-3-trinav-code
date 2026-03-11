from genedesign.seq_utils.reverse_complement import reverse_complement


class InternalRBSChecker:
    """
    Checks a DNA sequence for internal ribosome binding sites (RBS).

    An internal RBS is a Shine-Dalgarno-like sequence (e.g. AGGAGG or close
    variants) followed by a start codon (ATG) within 5-10 bp downstream.
    If such a site exists inside a CDS, the ribosome can initiate translation
    at the wrong position, producing truncated or aberrant proteins.

    The checker scans both the forward and reverse complement strands.

    Returns:
        (True, None)         if no internal RBS is found
        (False, site_seq)    if an internal RBS is found, with the offending window
    """

    def __init__(self):
        self.sd_patterns = []       # Shine-Dalgarno consensus sequences to scan for
        self.min_spacer = 4         # Minimum bp between SD end and ATG start
        self.max_spacer = 10        # Maximum bp between SD end and ATG start

    def initiate(self):
        """
        Sets up the Shine-Dalgarno consensus patterns to scan for.

        These are the core SD sequences known to promote ribosome binding in E. coli.
        The canonical SD is AGGAGG, but partial matches (AGGAG, GAGG, AAGG) are
        also recognized as functional by the ribosome.
        """
        self.sd_patterns = [
            "AGGAGG",   # canonical full SD
            "AGGAG",    # strong partial
            "GAGG",     # strong partial
            "AAGG",     # moderate partial
            "AGGA",     # moderate partial
            "GGAG",     # moderate partial
        ]

    def run(self, dnaseq: str) -> tuple[bool, str | None]:
        """
        Scans the sequence (and its reverse complement) for internal RBS sites.

        An internal RBS hit is defined as any SD pattern followed by ATG
        within self.min_spacer to self.max_spacer bases downstream.

        Parameters:
            dnaseq (str): The DNA sequence (CDS) to check.

        Returns:
            tuple: (bool, str or None)
                - True, None if no internal RBS found
                - False, offending_window if an internal RBS is found
        """
        if not dnaseq or not isinstance(dnaseq, str):
            raise ValueError(f"Invalid input: dnaseq must be a non-empty string, got {type(dnaseq)}")

        seq = dnaseq.upper()
        rc = reverse_complement(seq)

        # Check both strands; use 'x' separator so no window spans the junction
        combined = seq + "x" + rc

        for sd in self.sd_patterns:
            start = 0
            while True:
                pos = combined.find(sd, start)
                if pos == -1:
                    break

                # Look for ATG in the spacer window after this SD
                sd_end = pos + len(sd)
                for spacer in range(self.min_spacer, self.max_spacer + 1):
                    atg_pos = sd_end + spacer
                    if atg_pos + 3 <= len(combined):
                        if combined[atg_pos:atg_pos + 3] == "ATG":
                            # Return the full offending window: SD + spacer + ATG
                            offending = combined[pos:atg_pos + 3]
                            if 'x' not in offending:  # Ensure it's not crossing the junction
                                return False, offending

                start = pos + 1

        return True, None


if __name__ == "__main__":
    checker = InternalRBSChecker()
    checker.initiate()

    # Should detect internal RBS: AGGAGG + 6bp spacer + ATG
    bad_seq = "ATGCGTAAAGGAGGTTTTTTTATGCCCGTA"
    result, site = checker.run(bad_seq)
    print(f"Bad seq -> result: {result}, site: {site}")  # Expected: False

    # Clean sequence with no SD-like pattern near ATG
    good_seq = "ATGCGTAAACCCTTTTTTCCCGGTCTGCCC"
    result, site = checker.run(good_seq)
    print(f"Good seq -> result: {result}, site: {site}")  # Expected: True