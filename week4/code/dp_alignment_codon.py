import time


def read_fasta(path):
    """Read a FASTA file and return only the sequence."""
    with open(path) as f:
        lines = f.readlines()
    # Skip the first line if it starts with '>'
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    return "".join(line.strip() for line in lines if line.strip())


from typing import Dict  # Codon supports Python typing names (List/Tuple/Dict/etc.)

def parse_multi_fasta(filepath):
    """Parse a multi-sequence FASTA file into a dict like {'q1': 'AG', ...}"""
    sequences: Dict[str, str] = {}
    current_id = ""  # Use sentinel string instead of None
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                sequences[current_id] = ""
            else:
                sequences[current_id] += line
    return sequences




class Align:
    # Declare fields/types for Codon
    match: int
    mismatch: int
    gaps: int
    gap_open: int
    gap_extension: int
    last_score: int  # use int sentinel instead of None

    def __init__(self, match=3, mismatch=-3, gaps=-2, gap_open=-5, gap_extension=-1):
        # all stuff here

        # track scoring parameters
        self.match = match
        self.mismatch = mismatch
        self.gaps = gaps
        self.gap_open = gap_open
        self.gap_extension = gap_extension

        # most recent alignment score
        self.last_score = -1  # sentinel (Codon: avoid None)

    def _match_score(self, base1, base2):
        """Return score for comparing two bases."""
        return self.match if base1 == base2 else self.mismatch

    def global_alignment(self, s1, s2):
        '''implements Needleman Wunsch for global alignment
        only scoring matrix and final score needed, no full alignment traceback
        '''
        n, m = len(s1), len(s2)

        # initialize scoring matrix
        dp = [[0] * (m + 1) for _ in range(n + 1)]

        # add gap penalties
        for i in range(1, n + 1):
            dp[i][0] = dp[i - 1][0] + self.gaps
        for j in range(1, m + 1):
            dp[0][j] = dp[0][j - 1] + self.gaps

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if s1[i - 1] == s2[j - 1]:
                    match = dp[i - 1][j - 1] + self.match
                else:
                    match = dp[i - 1][j - 1] + self.mismatch

                delete = dp[i - 1][j] + self.gaps
                insert = dp[i][j - 1] + self.gaps

                dp[i][j] = max(match, delete, insert)

            #     print(
            #     f"i={i}, j={j}, compare {s1[i-1]} vs {s2[j-1]}, "
            #     f"match={match}, delete={delete}, insert={insert}, "
            #     f"chosen={dp[i][j]}"
            # )

        # final score
        self.last_score = dp[n][m]
        return self.last_score

    def local_alignment(self, s1, s2):
        '''implements Smith–Waterman for local alignment

        '''
        n, m = len(s1), len(s2)

        # initialize scoring matrix (all zeros)
        dp = [[0] * (m + 1) for _ in range(n + 1)]

        max_score = 0  # keep track of the highest score found

        # fill matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if s1[i - 1] == s2[j - 1]:
                    match = dp[i - 1][j - 1] + self.match
                else:
                    match = dp[i - 1][j - 1] + self.mismatch

                delete = dp[i - 1][j] + self.gaps
                insert = dp[i][j - 1] + self.gaps

                # local alignment rule: no negative values
                dp[i][j] = max(0, match, delete, insert)

                # track the maximum score across the matrix
                if dp[i][j] > max_score:
                    max_score = dp[i][j]

        self.last_score = max_score
        return self.last_score

    def semi_global_alignment(self, s1, s2):
        '''implements semi-global (fitting) alignment

        - gaps at beginning of short sequence are penalized
        - gaps at beginning of long sequence are free
        - best alignment score in last column of query alignment
        ex.
            ATCCTC
            --CC
        '''

        if len(s1) >= len(s2):
            ref = s1  # longer
            q = s2  # shorter
        else:
            ref = s2
            q = s1

        n, m = len(ref), len(q)
        dp = [[0] * (m + 1) for _ in range(n + 1)]

        # penalize gaps in query q
        for j in range(1, m + 1):
            dp[0][j] = dp[0][j - 1] + self.gaps

        # Initialize first column free gaps in ref
        for i in range(1, n + 1):
            dp[i][0] = 0

        # fill
        for i in range(1, n + 1):

            for j in range(1, m + 1):
                if ref[i - 1] == q[j - 1]:
                    diag = dp[i - 1][j - 1] + self.match
                else:
                    diag = dp[i - 1][j - 1] + self.mismatch

                up = dp[i - 1][j] + self.gaps
                left = dp[i][j - 1] + self.gaps

                dp[i][j] = max(diag, up, left)

        # The best alignment ends anywhere along the last column - fully aligned query
        # best alignment score occurs anywhere along last column
        best_score = dp[0][m]
        for i in range(1, n + 1):
            if dp[i][m] > best_score:
                best_score = dp[i][m]

        self.last_score = best_score
        return self.last_score



    def affine_gap_penalty_global_alignment(self, s1, s2):
        '''
        implements affine-gap Needleman Wunsch
        pay extra for opening a gap and less for extending it
        '''

        n, m = len(s1), len(s2)
        NEG_INF = -10**12  # use large negative int to keep matrices homogeneous

        # Matrices: M = match/mismatch, X = gap in s2 (vertical), Y = gap in s1 (horizontal)
        M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
        X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
        Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

        # Initialization
        M[0][0] = 0
        for i in range(1, n + 1):
            X[i][0] = self.gap_open + (i - 1) * self.gap_extension
            M[i][0] = X[i][0]
        for j in range(1, m + 1):
            Y[0][j] = self.gap_open + (j - 1) * self.gap_extension
            M[0][j] = Y[0][j]

        # Fill matrices
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                score = self.match if s1[i - 1] == s2[j - 1] else self.mismatch

                # Match/mismatch
                best_prev = M[i - 1][j - 1]
                if X[i - 1][j - 1] > best_prev:
                    best_prev = X[i - 1][j - 1]
                if Y[i - 1][j - 1] > best_prev:
                    best_prev = Y[i - 1][j - 1]
                M[i][j] = best_prev + score

                # Gap in s2 (vertical gap)
                open_gap = M[i - 1][j] + self.gap_open
                extend_gap = X[i - 1][j] + self.gap_extension
                X[i][j] = open_gap if open_gap > extend_gap else extend_gap

                # Gap in s1 (horizontal gap)
                open_gap_h = M[i][j - 1] + self.gap_open
                extend_gap_h = Y[i][j - 1] + self.gap_extension
                Y[i][j] = open_gap_h if open_gap_h > extend_gap_h else extend_gap_h

        # Final score: best among M, X, Y at bottom-right
        best_final = M[n][m]
        if X[n][m] > best_final:
            best_final = X[n][m]
        if Y[n][m] > best_final:
            best_final = Y[n][m]
        self.last_score = best_final
        return self.last_score




if __name__ == "__main__":

    score = Align(match=3, mismatch=-3, gaps=-2, gap_open=-5, gap_extension=-1)

    # x = score.global_alignment("AG", "AGC")
    # print(score.last_score)  # shows the final score of that alignment

    # -----------------------

    # only remove first line
    MT_human = read_fasta("../MT-human.fa")
    MT_orang = read_fasta("../MT-orang.fa")

    # parse and remove > for q and t fa files
    # assing sequence to each q1,2,3,4,5 and same for t
    qs = parse_multi_fasta("../q1.fa")
    ts = parse_multi_fasta("../t1.fa")

    # test output
    # print(qs["q2"])
    # print(ts["t5"])

    # redundant
    start = time.time()
    result = score.global_alignment(MT_human, MT_orang)

    elapsed = time.time() - start
    elapsed_ms = elapsed * 1000
    # main output 
    print(f"{'global' + '-' + 'mt_human':<20} {'codon':<12} {elapsed_ms:7.2f}")

    start = time.time()
    result = score.local_alignment(MT_human, MT_orang)

    elapsed = time.time() - start
    elapsed_ms = elapsed * 1000
    # main output 
    print(f"{'local' + '-' + 'mt_human':<20} {'codon':<12} {elapsed_ms:7.2f}")

    start = time.time()
    result = score.semi_global_alignment(MT_human, MT_orang)

    elapsed = time.time() - start
    elapsed_ms = elapsed * 1000
    # main output 
    print(f"{'semi-global' + '-' + 'mt_human':<20} {'codon':<12} {elapsed_ms:7.2f}")

    tests = {
        "q1": (qs["q1"], ts["t1"]),
        "q2": (qs["q2"], ts["t2"]),
        "q3": (qs["q3"], ts["t3"]),
        "q4": (qs["q4"], ts["t4"]),
        "q5": (qs["q5"], ts["t5"])}

    # methods to test
    methods = ["global", "local", "semi-global", "affine"]

    # run all tests
    for name, (s1, s2) in tests.items():

        for method_name in methods:

            # super redundant but I won't do anything about it 
            if method_name == "global": 
                start = time.time()
                result = score.global_alignment(s1, s2)
                elapsed_ms = (time.time() - start) * 1000

            elif method_name == "local": 
                start = time.time()
                result = score.local_alignment(s1, s2)
                elapsed_ms = (time.time() - start) * 1000
            
            elif method_name == "semi-global": 
                start = time.time()
                result = score.semi_global_alignment(s1, s2)
                elapsed_ms = (time.time() - start) * 1000

            elif method_name == "affine": 
                start = time.time()
                result = score.affine_gap_penalty_global_alignment(s1, s2)
                elapsed_ms = (time.time() - start) * 1000


            # main output 
            print(f"{method_name + '-' + name:<20} {'codon':<12} {elapsed_ms:7.2f}")


            # optional debug score output
            # print(f"{method_name}: {result}")


