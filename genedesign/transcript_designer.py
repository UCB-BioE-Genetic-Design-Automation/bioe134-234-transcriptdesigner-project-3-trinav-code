import random
import csv
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker


class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and selects an RBS.

    Uses a sliding window codon selection strategy: codons are sampled one at a time
    weighted by RSCU frequency, and a local hairpin check is performed at each step.
    If the newly added codon creates a hairpin in the local window, it is resampled
    immediately before proceeding. This prevents hairpin problems from accumulating
    across the sequence. A final pass of all checkers validates the complete CDS,
    with full resampling retries if needed.
    """

    def __init__(self):
        self.codonTable = {}
        self.rbsChooser = None
        self.forbiddenChecker = None
        self.promoterChecker = None
        self.codonChecker = None
        self.internalRBSChecker = None
        self.rng = random.Random(42)

    def initiate(self) -> None:
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        self.forbiddenChecker = ForbiddenSequenceChecker()
        self.forbiddenChecker.initiate()

        self.promoterChecker = PromoterChecker()
        self.promoterChecker.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.internalRBSChecker = InternalRBSChecker()
        self.internalRBSChecker.initiate()

        codon_usage_file = 'genedesign/data/codon_usage.txt'
        raw = {}
        with open(codon_usage_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 4:
                    continue
                codon = row[0].strip()
                aa = row[1].strip()
                freq = float(row[2].strip())
                if freq > 0:
                    raw[codon] = (aa, freq)

        for codon, (aa, freq) in raw.items():
            if aa not in self.codonTable:
                self.codonTable[aa] = []
            self.codonTable[aa].append((codon, freq))

    def _sliding_window_sample(self, peptide: str) -> list[str]:
        """
        Builds a codon list one codon at a time, checking for hairpins locally
        at each step using a sliding window of the last 20 nucleotides.

        At each position:
        1. Sample a codon weighted by RSCU frequency
        2. Check if the tail of the CDS built so far contains a hairpin
        3. If yes, resample this codon (up to max_codon_attempts times)
        4. If still failing after max attempts, accept best and move on

        The window size (20 nt) covers the minimum hairpin length so we catch
        problems as soon as they're introduced rather than after full assembly.
        """
        window_nt = 20          # nucleotide window to check for local hairpins
        max_codon_attempts = 10 # resampling attempts per codon position

        codons = []

        for aa in peptide:
            options = self.codonTable.get(aa)
            if not options:
                codons.append("NNN")
                continue

            codons_list = [c for c, _ in options]
            weights = [w for _, w in options]

            best_codon = None
            for attempt in range(max_codon_attempts):
                candidate = self.rng.choices(codons_list, weights=weights, k=1)[0]

                # Check local window: tail of built sequence + this candidate codon
                built_so_far = ''.join(codons) + candidate
                local_window = built_so_far[-window_nt:]

                ok, _ = hairpin_checker(local_window)
                if ok:
                    best_codon = candidate
                    break
                elif best_codon is None:
                    # Keep first attempt as fallback in case none pass
                    best_codon = candidate

            codons.append(best_codon)

        codons.append("TAA")
        return codons

    def _count_failures(self, codons: list[str]) -> int:
        """Runs all checkers and returns the number that failed (0 = perfect)."""
        cds = ''.join(codons)
        failures = 0
        if not self.forbiddenChecker.run(cds)[0]:
            failures += 1
        if not hairpin_checker(cds)[0]:
            failures += 1
        if not self.promoterChecker.run(cds)[0]:
            failures += 1
        if not self.codonChecker.run(codons)[0]:
            failures += 1
        if not self.internalRBSChecker.run(cds)[0]:
            failures += 1
        return failures

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide to a codon-optimized DNA sequence and selects an RBS.

        Strategy:
        1. Build codon sequence using sliding window with inline hairpin checking
        2. Run all checkers on the complete CDS
        3. If passing, return immediately; otherwise track best and retry
        4. Return the best candidate after max_attempts
        """
        max_attempts = 100
        best_codons = None
        best_failures = float('inf')

        for _ in range(max_attempts):
            codons = self._sliding_window_sample(peptide)
            failures = self._count_failures(codons)

            if failures < best_failures:
                best_failures = failures
                best_codons = codons

            if best_failures == 0:
                break

        cds = ''.join(best_codons)
        selectedRBS = self.rbsChooser.run(cds, ignores)
        return Transcript(selectedRBS, peptide, best_codons)


if __name__ == "__main__":
    peptide = "MYPFIRTARMTV"
    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()
    transcript = designer.run(peptide, ignores)
    print(transcript)