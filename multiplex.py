##!! do comments on this code and multiplex code in primer_designer.py

import itertools
from typing import List, Dict, Tuple
from primer_config import MINIMUM_SPACING, MELT_TM_TOLERANCE, MAX_CROSS_DIMERS
from src.thermo import deltaG

def _parse_pos(pos: str) -> Tuple[int,int]:
    a, b = map(int, pos.split('-'))
    return a-1, b

def _check_tm(pairs: List[Dict]) -> bool:
    tms = [p['forward']['tm'] for p in pairs] + [p['reverse']['tm'] for p in pairs]
    return max(tms) - min(tms) <= MELT_TM_TOLERANCE

def _check_spacing(pairs: List[Dict]) -> bool:
    pairs_sorted = sorted(pairs, key=lambda p: p['target'])
    for target, group in itertools.groupby(pairs_sorted, key=lambda p: p['target']):
        recs = []
        for p in group:
            fs, fe = _parse_pos(p['forward']['position'])
            rs, re = _parse_pos(p['reverse']['position'])
            recs.append((min(fs,rs), max(fe,re)))
        recs.sort()
        for (s1,e1),(s2,e2) in zip(recs, recs[1:]):
            if s2 - e1 < MINIMUM_SPACING:
                return False
    return True

def _check_dimers(pairs: List[Dict], cross_dimer_score) -> bool:
    for a, b in itertools.combinations(pairs, 2):
        d1 = deltaG(a['forward']['sequence'], b['reverse']['sequence'])
        d2 = deltaG(a['reverse']['sequence'], b['forward']['sequence'])
        if abs(d1) > cross_dimer_score or abs(d2) > cross_dimer_score:
            print(f"   cross-dimer fail: |d1|={abs(d1):.1f}, |d2|={abs(d2):.1f} > {cross_dimer_score}")
            return False
    return True

def select_multiplex_sets(all_pairs: List[Dict],
                          cross_dimer_score,
                          max_panels: int = 3) -> List[List[Dict]]:
    candidates = []
    for size in range(2, min(6, len(all_pairs) + 1)):
        for combo in itertools.combinations(all_pairs, size):
            # ensure exactly one primer‐pair per target in each panel
            if len({p['target'] for p in combo}) != size:
                continue

            if _check_tm(combo) and \
               _check_spacing(combo) and \
               _check_dimers(combo, cross_dimer_score):
                score = sum(p['pair_score'] for p in combo) / size
                candidates.append((list(combo), score))

    candidates.sort(key=lambda x: x[1], reverse=True)
    return [panel for panel, _ in candidates[:max_panels]]