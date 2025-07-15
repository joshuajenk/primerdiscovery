import math
from primer_config import Default_PRIMER_CONC
from primer_config import Default_NA
from primer_config import Default_MG
from primer_config import Default_dNTP


from Bio.Seq import reverse_complement

# nearest-neighbor approx. for delH & delS
NN_PARAMS = {
    'AA':(-7.9, -22.2), 'AT':(-7.2, -20.4), 'TA':(-7.2, -21.3), 'CA':(-8.5, -22.7),
    'GT':(-8.4, -22.4), 'CT':(-7.8, -21.0), 'GA':(-8.2, -22.2), 'CG':(-10.6,-27.2),
    'GC':(-9.8, -24.4), 'GG':(-8.0, -19.9)
}

# NOTE: dangling end considerations for primer positioning to amplicon
# init & symmetry; dangling end bonus
INIT = {'G_or_C': (0.1, -2.8), 'A_or_T': (2.3, 4.1)}
SYMMETRY_CORRECTION = -0.4  # kcal/mol

R = 1.987 # Universal Gas Constant (cal / (k*mol))

import re

# clean seq
def clean_sequence(seq):
    # rem whitespace, digits, and non-letters
    seq = re.sub(r'[^A-Za-z]', '', seq.upper())
    # drop all headers
    if seq.startswith('>'):
        seq = ''
    return seq

# calculate advanced melt temperature
# SantaLucia Method (1998)
def calc_tm_nn(seq, primer_conc = Default_PRIMER_CONC, def_na = Default_NA, def_mg = Default_MG, def_dntp = Default_dNTP, temp_C = None):
    seq = seq.upper()
    H, S = 0.0, 0.0

    # 5' and 3' initialization
    for end in (seq[0], seq[-1]):
        k = 'A_or_T' if end in 'AT' else 'G_or_C'
        h0, s0 = INIT[k]
        H += h0; S += s0

    # nearest-neighbor approx
    for i in range(len(seq)-1):
        nn = seq[i:i+2]
        if nn in NN_PARAMS:
            dh, ds = NN_PARAMS[nn]
            H += dh; S += ds

    # symmetry correction
    if seq == reverse_complement(seq):
        H += SYMMETRY_CORRECTION
        S -= 1.4 # symmetry entropy corr

    # mM to Molar calc
    na_m = def_na / 1000.0
    mg_m = max((def_mg - 0.9 * def_dntp), 0.01) / 1000.0

    # Salt correction; SantaLucia (1998)
    eff_salt = na_m + 120 * math.sqrt(mg_m)
    salt_corr = 0.368 * (len(seq) - 1) * math.log10(eff_salt)
    S += salt_corr

    # delS (def = Tm estimation)
    if temp_C is not None:
        T = temp_C + 273.15
    else:
        T = (1000 * H) / (S + R * math.log(primer_conc/4.0)) #init. estimate

    # Tm calc
    tm_K = (1000 * H) / (S + R * math.log(primer_conc/4.0))
    tm_C = tm_K - 273.15
    return round(tm_C, 1)

def deltaG(seq1, seq2=None, na_conc = Default_NA, mg_conc = Default_MG):
    # duplex formation
    # if seq-2 none --> computes self-delG (hairpin/dimer potential)
    # align seq
    targ = seq2 if seq2 else seq1
    seq1 = seq1.upper(); targ = targ.upper()


    # delG arrange for scoring
    best_dG = 0
    for offset in range(-len(targ) + 1, len(seq1)):
        dG = 0
        for i in range(len(targ)):
            j = i + offset
            if 0 <= j < len(seq1):
                pair = targ[i] + seq1[j]
                if pair in NN_PARAMS:
                    dh, _ = NN_PARAMS[pair]
                    dG += dh
        best_dG = min(best_dG, dG)
    return round(best_dG, 1)

# reverse complement join
def reverse_complement(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G'}
    return ''.join(comp.get(b,'N') for b in seq[::-1])
