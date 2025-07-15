# src/primer3_wrapper.py
import primer3
import re
from primer_config import (
  Primer3_PRODUCT_SIZE_RANGE,
  Primer3_PRIMER_SIZE_RANGE,
  Def_Melt_Tm_LOWER, Def_Melt_Tm_UPPER,
  Def_GC_LOWER, Def_GC_UPPER
)

# Primer3.py search query
# implemented as first round
def pick_with_primer3(seq, na_conc, mg_conc):
    # mask non-ACGT bases; N.A. pref
    seq = re.sub(r'[^ACGTN]', 'N', seq.upper())
    L = len(seq)

    # FLAG?
    # adapt min product size to template length
    min_prod, max_prod = Primer3_PRODUCT_SIZE_RANGE[0]
    if L < min_prod:
        # allow half‐amplicon or 30 bp
        min_prod = max(1, min(L//2, min_prod))
    product_range = [(min_prod, max_prod)]

    # params, primer3.py recommendations
    params = {
        'PRIMER_OPT_TM': (Def_Melt_Tm_LOWER + Def_Melt_Tm_UPPER)//2,
        'PRIMER_MIN_TM': Def_Melt_Tm_LOWER,
        'PRIMER_MAX_TM': Def_Melt_Tm_UPPER,
        'PRIMER_OPT_GC_PERCENT': (Def_GC_LOWER + Def_GC_UPPER)//2,
        'PRIMER_MIN_GC': Def_GC_LOWER,
        'PRIMER_MAX_GC': Def_GC_UPPER,
        'PRIMER_PRODUCT_SIZE_RANGE': product_range,
        'PRIMER_MIN_SIZE': Primer3_PRIMER_SIZE_RANGE[0],
        'PRIMER_MAX_SIZE': Primer3_PRIMER_SIZE_RANGE[1],
        'PRIMER_SALT_MONOVALENT': na_conc,
        'PRIMER_SALT_DIVALENT': mg_conc,
        # warnings per primer3.py
        'PRIMER_MAX_POLY_X': 3,
        'PRIMER_MAX_SELF_ANY': 8,
        'PRIMER_MAX_SELF_END': 3
    }

    try:
        return primer3.bindings.design_primers(
            {'SEQUENCE_TEMPLATE': seq},
            params
        )

    # run OS error for small templates
    except OSError as e:
        msg = str(e)
        if 'SEQUENCE_INCLUDED_REGION' in msg:
            print(f"\nERROR: Template length ({L} bp) is smaller than min product size ({min_prod} bp).")
            return {}   # or return
        raise