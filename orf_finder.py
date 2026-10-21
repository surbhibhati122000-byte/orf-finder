def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


bases = "TCAG"
amino = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
i = 0
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            CODON_TABLE[b1 + b2 + b3] = amino[i]
            i += 1

STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate(seq):
    return "".join(CODON_TABLE.get(seq[i:i+3], "X") for i in range(0, len(seq) - 2, 3))

def find_orfs(seq, frame, min_len_aa=100):
    orfs = []
    start = None
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if start is None and codon == "ATG":
            start = i
        elif start is not None and codon in STOP_CODONS:
            length_aa = (i - start) // 3
            if length_aa >= min_len_aa:
                orfs.append((start, i + 3, length_aa))
            start = None
    return orfs


def find_all_orfs(seq, min_len_aa=100):
    L = len(seq)
    results = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            for start, end, n_aa in find_orfs(s, frame, min_len_aa):
                protein = translate(s[start:end]).rstrip("*")
                if strand == "+":
                    f_start, f_end = start + 1, end
                else:
                    f_start, f_end = L - end + 1, L - start
                results.append({
                    "strand": strand, "frame": frame + 1,
                    "start": f_start, "end": f_end,
                    "length_aa": n_aa, "protein": protein,
                })
    return results


def read_fasta(path):
    parts = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                parts.append(line.strip().upper())
    return "".join(parts)
if __name__ == "__main__":
    print(reverse_complement("ATGC"))
    print(translate("ATGAAATAG"))
    print(find_all_orfs("ATGAAATAG", min_len_aa=1))