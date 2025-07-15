from primer_config import Def_Melt_Tm_LOWER, Def_Melt_Tm_UPPER, Def_GC_LOWER, Def_GC_UPPER

# ASCII codes
RED = '\033[91m'
YELLOW = '\033[93m'
GREEN = '\033[92m'
RESET = '\033[0m'
CYAN = '\033[96m'
BOLD = '\033[1m'
END = '\033[0m'
BOLD_GREEN = '\033[1m\033[92m'
BOLD_BLUE = '\033[1m\033[94m'

# color GC
def color_gc(gc):
    if gc < Def_GC_LOWER:
        return f"{YELLOW}{gc}% (LOW){RESET}"
    elif gc > Def_GC_UPPER:
        return f"{RED}{gc}% (HIGH){RESET}"
    else:
        return f"{GREEN}{gc}%{RESET}"

# color Tm
def color_tm(tm):
    if tm < Def_Melt_Tm_LOWER:
        return f"{YELLOW}{tm}°C (LOW){RESET}"
    elif tm > Def_Melt_Tm_UPPER:
        return f"{RED}{tm}°C (HIGH){RESET}"
    else:
        return f"{GREEN}{tm}°C{RESET}"

# color score
def color_score(score):
    if score >= 90:
        return f"{GREEN}{score} / 100 (EXCELLENT){RESET}"
    elif score >= 70:
        return f"{YELLOW}{score} / 100 (MODERATE){RESET}"
    else:
        return f"{RED}{score} / 100 (LOW){RESET}"

# rank QC
def qc_rank(score):
    if score >= 90:
        return "Ideal"
    elif score >= 70:
        return "Moderate"
    else:
        return "Risky"

# output for each pair
def print_optimized_primers(pairs, highlight=True):
    for idx, pair in enumerate(pairs, 1):
        print()
        print("="*60)
        score = pair.get('pair_score', 0.0)
        score_colored = color_score(score)
        rank = qc_rank(score)
        print(f"{BOLD}OPTIMIZED PRIMER PAIR#{idx}{END} - SCORE: {score_colored} — {rank}")
        print("-" * 60)

        # GC% color filters
        f_gc_colored = color_gc(pair['forward']['gc_content'])
        r_gc_colored = color_gc(pair['reverse']['gc_content'])

        # colored Tm values
        f_tm_colored = color_tm(pair['forward']['tm'])
        r_tm_colored = color_tm(pair['reverse']['tm'])

        # print statements w/ color
        print(f"Forward Primer: {pair['forward']['sequence']} (Tm: {f_tm_colored}, GC: {f_gc_colored})")
        print(f"Reverse Primer: {pair['reverse']['sequence']} (Tm: {r_tm_colored}, GC: {r_gc_colored})")

        try:
            f_start = int(pair['forward']['position'].split('-')[0])
            r_end = int(pair['reverse']['position'].split('-')[1])
            size = abs(f_start - r_end)
            if f_start is not None and r_end is not None:
                print(f"Amplicon size: {abs(f_start - r_end)} bp\n")
            else:
                print("Amplicon size: N/A\n")
        except (KeyError, ValueError, IndexError):
            print("Amplicon Size:   N/A")
        tm_diff = pair.get('tm_difference', 'N/A')
        print(f"Tm Difference:   {tm_diff}°C")
        print(f"Hetero-Dimers:   {pair.get('hetero_dimers', 'N/A')}")
        print("=" * 60)

# comparison of primer pairs to optimized genes
def visualize_amplicon(seq, f_start, f_end, r_start, r_end, wrap=60, idx=None):
    print(f"\n{CYAN}[Visualization for Optimized Pair #{idx}]{END}")
    print(f"{CYAN}Amplicon visualization (5' → 3'){END}")

    for i in range(0, len(seq), wrap):
        line = ''
        for j in range(i, min(i + wrap, len(seq))):
            base = seq[j]
            if f_start <= j < f_end:
                line += f"{BOLD_GREEN}{base}{END}"
            elif r_start <= j < r_end:
                line += f"{BOLD_BLUE}{base}{END}"
            else:
                line += base
        print(line)

    print(f"{BOLD_GREEN}← Forward Primer:{END} {f_start+1}-{f_end}")
    print(f"{BOLD_BLUE}→ Reverse Primer:{END} {r_start+1}-{r_end}")
    print(f"Amplicon size: {abs(f_start - r_end)} bp\n")