# src/qc.py

# filter primer3 results
def filter_primer3_results(results):
    num = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
    if num == 0:
        print("[QC] Primer3 returned 0 pairs → skipping QC")
        return []

    # store pairs
    pairs = []
    for i in range(num):
        fwd = {
            'SEQUENCE': results[f'PRIMER_LEFT_{i}_SEQUENCE'],
            'Tm':       results[f'PRIMER_LEFT_{i}_TM'],
            'GC':       results[f'PRIMER_LEFT_{i}_GC_PERCENT'],
        }
        rev = {
            'SEQUENCE': results[f'PRIMER_RIGHT_{i}_SEQUENCE'],
            'Tm':       results[f'PRIMER_RIGHT_{i}_TM'],
            'GC':       results[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
        }
        pairs.append({'forward': fwd, 'reverse': rev})

    print(f"[QC] Passing through all {num} Primer3 pairs (no filtering)")
    return pairs