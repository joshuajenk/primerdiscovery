# imports
from thefuzz import process
from Bio import Entrez
from src.primer_output import color_tm, color_gc, color_score
import pandas as pd
from src.multiplex import select_multiplex_sets
from src.thermo import calc_tm_nn, deltaG
from primer_config import MAX_CROSS_DIMERS

# customs import
import re

# ANSI palettes:
BLUE  = '\033[94m'
GREEN = '\033[92m'
RESET = '\033[0m'

# clean seq
def clean_sequence(seq):
    seq = re.sub(r'[^ACGTacgt]', '', seq)
    return seq.upper()

# configs --> view primer_config.py to change parameters
from primer_config import (
    Def_Melt_Tm_LOWER,
    Def_Melt_Tm_UPPER,
    Def_GC_LOWER,
    Def_GC_UPPER,
    Def_Hairpin_Min,
    Def_SelfDimer_Min,
    Default_MG,
    Default_NA,
    Default_dNTP
)

# PRIMER3.PY: import tool for selecting optimal primers for DNA amplification, ensuring specificity and efficiency in PCR experiments.
try:
    import primer3
    PRIMER3_AVAILABLE = True
except ImportError:
    PRIMER3_AVAILABLE = False

# Bio.Blast: BLAST for Searching... handles XML and tabular inputs
try:
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True

# Fallback Mode
except ImportError:
    BIOPYTHON_AVAILABLE = False
    NCBIWWW = None
    NCBIXML = None

# main class
class PrimerDesigner:
    def __init__(self, pathogen_sequences, ss_to_ls, fast_mode=True):
        self.fast_mode = fast_mode
        self.pathogen_sequences = pathogen_sequences
        self.ss_to_ls = ss_to_ls
        # define compliments
        self.complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

        # initialize NCBI
        self.ncbi_databases = {
            'nt': 'Nucleotide collection (nr/nt)',
            'refseq_rna': 'NCBI Transcript Reference Sequences',
            'refseq_genomic': 'NCBI Genomic Reference Sequences',
            'chromosome': 'Chromosome sequences'
        }

        # NCBI BLAST endpts
        self.blast_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
        self.blast_results_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

        # problematic sequences
        self.common_sequences = self._load_common_sequences() #understand more about this

    def run_cli(self):
        mode = input("\nSelect design mode [single(s) / multiplex(m)]: ").strip().lower()
        if mode in ('multi', 'm', 'multiplex'):
            mode = 'multiplex'
        else:
            mode = 'singleplex'
        print(f"\nSelected {mode} assay.")

        def get_sequence():
            print("\n - Available Pathogens:")
            for n in sorted(self.pathogen_sequences):
                print(f" - {n}")
            print(" - Custom")
            print(" - NCBI Accession")
            choice = input("Enter Pathogen Name / Shorthand, Custom, or NCBI Accession: ").strip()

            if choice.lower().startswith('custom'):
                raw = input("Enter custom sequence (50–1500 bp; DNA or AA): ").strip()
                if not raw:
                    print("No sequence entered. Exiting…")
                    return None, None
                if self.is_amino_acid_sequence(raw):
                    seq = self.amino_acid_to_dna(raw)
                else:
                    seq = clean_sequence(raw)
                if len(seq) < 50 or len(seq) > 1500:
                    print("sequence length must be 50–1500 bp.")
                    return None, None
                return "custom sequence", seq

            if re.match(r'^[A-Z]{1,2}_\d+(\.\d+)?$', choice, re.I):
                try:
                    seq, hdr = self.fetch_ncbi_sequence(accession=choice)
                    return f"ncbi {choice}", seq
                except Exception as e:
                    print(f"error fetching {choice}: {e}")
                    return None, None

            full = self.ss_to_ls.get(choice)
            if not full:
                full_name = self.fuzzy_match_pathogen(choice)  # score ≥80%
            if not full:
                full_name = choice

            seq = self.pathogen_sequences.get(full)
            if not seq:
                print(f"no sequence found for '{choice}'.")
                return
            return full, clean_sequence(seq)

        names, seqs = [], []
        if mode == 'multiplex':
            count = int(input("how many sequences to multiplex?: ").strip())
            for i in range(count):
                name, seq = get_sequence()
                if seq is None:
                    return
                names.append(name)
                seqs.append(seq)
        else:
            name, seq = get_sequence()
            if seq is None:
                return
            names.append(name)
            seqs.append(seq)

        all_pairs = []
        for name, seq in zip(names, seqs):
            print(f"\nsequence for {name} loaded. ({len(seq)} bp)")
            res = self.design_primers(seq)
            if mode == 'singleplex':
                return res

            # tag and collect optimized pairs
            optimized = res.get('optimized_pairs', [])
            for p in optimized:
                p['target'] = name
            all_pairs.extend(optimized)


        if mode == 'multiplex':
            if not all_pairs:
                print("no optimized primer pairs available for multiplexing.")
                return

            panels = select_multiplex_sets(
                all_pairs=all_pairs,
                cross_dimer_score=MAX_CROSS_DIMERS,
                max_panels=5
            )

            if not panels:
                print("No compatible multiplex panels found.")
                return

            # multiplex design
            print("\nMultiplex Panel Designs:")
            for idx, panel in enumerate(panels, 1):
                # compute panel metrics
                tm_vals = [(p['forward']['tm'] + p['reverse']['tm']) / 2 for p in panel]
                gc_vals = [(p['forward']['gc_content'] + p['reverse']['gc_content']) / 2 for p in panel]
                # tm and GC% of ea. panel
                panel_tm = round(sum(tm_vals) / len(tm_vals), 1)
                panel_gc = round(sum(gc_vals) / len(gc_vals), 1)
                panel_score = sum(p['pair_score'] for p in panel) / len(panel)

                # header
                print(f"\n{BLUE}Panel {idx}{RESET} — Score: {color_score(panel_score)}")
                # each primer‐pair
                for p in panel:
                    f = p['forward']
                    r = p['reverse']
                    # color coded fwd / rev
                    print(f"  • {BLUE}FWD{RESET} {BLUE}{f['sequence']}{RESET}"
                          f" (Tm:{color_tm(f['tm'])} GC:{color_gc(f['gc_content'])} Pos:{f['position']})")
                    print(f"    {GREEN}REV{RESET} {GREEN}{r['sequence']}{RESET}"
                          f" (Tm:{color_tm(r['tm'])} GC:{color_gc(r['gc_content'])} Pos:{r['position']})")
                # panel summary
                print(f"  Avg Tm: {panel_tm}°C | Avg GC: {panel_gc}%\n")

        # GC% filtering and warnings for target
                gc = round((seq.count('G')+seq.count('C'))/len(seq)*100, 1)
                if gc < 70 or gc > 20:
                    gc = round((seq.count('G') + seq.count('C')) / len(seq) * 100, 1)
                    print(f"GC% Content: {gc}%")
                else:
                    gc = round((seq.count('G') + seq.count('C')) / len(seq) * 100, 1)
                    print(f"GC% Content -- {gc}% -- Usual Result, Check Sequence")

        # aa to na
        if self.is_amino_acid_sequence(seq):
            print("Amino acid sequence detected—converting to DNA.")
            sequence = self.amino_acid_to_dna(seq)
            print(f"Converted sequence length: {len(sequence)} bp")
        return None

    # keyword check for ss codes; must be in SECOND (B) column of sheet
    def fuzzy_match_pathogen(self, user_input, threshold=80):
        if user_input in self.ss_to_ls:
            return self.ss_to_ls[user_input]
        choices = list(self.pathogen_sequences.keys())
        match, score = process.extractOne(user_input,choices)
        if match in self.pathogen_sequences and score >= threshold:
            return match
        return None

    # fetch NCBI accession codes
    def fetch_ncbi_sequence(self, accession, email="biologydept@vanderbilt.edu"):
        Entrez.email = email
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            fasta = handle.read()
            handle.close()
            # clean sequence
            lines = fasta.splitlines()
            header = lines[0].strip() if lines else "UNKNOWN GENE"
            seq = ''.join(line for line in lines if not line.startswith('>'))
            return clean_sequence(seq), header
        except Exception as e:
            raise RuntimeError(f"Failed to fetch sequence: {e}")

        # ***FALLBACK MODE***
        # if NCBI down, use fallback as problematic seq:
    def _load_common_sequences(self):
        # common seq: recognize and note any significant seq. in tx
        common_seqs = [
            # Ribosomal RNA sequences
            # keep in storage
            "GTGCCAGCAGCCGCGGTAAT",  # 16S rRNA
            "GGCTACCTTGTTACGACTT",   # 16S rRNA
            "CCTACGGGAGGCAGCAGT",    # 16S rRNA
            "ACTCCTACGGGAGGCAGCAGT", # 16S rRNA
            "GTATTACCGCGGCTGCTGGCAC", # 16S rRNA

            # 18S rRNA sequences
            # keep in storage
            "ACCTGGTTGATCCTGCCAG",   # 18S rRNA
            "TGATCCTTCTGCAGGTTCACCTAC", # 18S rRNA
            "CTTGTTACGACTTTTACTTCC",    # 18S rRNA

            # Common Repeat Sequences
            # keep and flag
            "ATATATATATATATATATATAT",  # AT repeats
            "GCGCGCGCGCGCGCGCGCGCGC",  # GC repeats
            "CACACACACACACACACACACA",  # CA repeats
            "TGTGTGTGTGTGTGTGTGTGTG",  # TG repeats
            "AAAAAAAAAAAAAAAAAAAAAAA", # Poly-A
            "TTTTTTTTTTTTTTTTTTTTTT",  # Poly-T
            "GGGGGGGGGGGGGGGGGGGGGG",  # Poly-G
            "CCCCCCCCCCCCCCCCCCCCCC",  # Poly-C

            # Mitochondrial sequences
            # keep in storage
            "TAAACTTCAGGGTGACCAAAAAATCA", # mt-CO1
            "CCGGTTTGAACTCAGATCATGT",     # mt-CO1
            "GCTCGTGTTCTTATTAATAATCTAT",  # mt-CO1

            # Common vector sequences
            # keep in storage
            "TCGAGGTGAAGATGAAAGCGC",      # pUC origin
            "CTTGGCGCGCCAAGCTTGCGC",      # Multiple cloning site
            "GGATCCTCTAGAGTCGACCTGCAG",   # Common restriction sites
        ]
        return common_seqs

    # run NCBI for fallback mode
    def ncbi_blast_search(self, sequence, database='nt', organism='', max_hits=50,
                         expect_threshold=1000, word_size=7, timeout=300):

        # Bio-Python is 2nd motion; fallback will initiate if 1) not avail. or 2) no results from 1-2
        if not BIOPYTHON_AVAILABLE:
            print("Biopython not available: Fallback")
            return self._fallback_specificity_check(sequence)

        try:
            # NCBI BLAST via Biopython
            print(f"Searching NCBI {database} database for primer: {sequence[:15]}...")

            # Parameters for single-seq search
            # all parameters per NCBI Entrez query
            blast_params = {
                'database': database,
                'sequence': sequence,
                'program': 'blastn',
                'expect': expect_threshold,
                'word_size': word_size,
                'hitlist_size': max_hits,
                'filter': 'L',  # Low complexity filter
                'format_type': 'XML'
            }

            # further spec. for Entrez query
            if organism:
                blast_params['entrez_query'] = f'"{organism}"[Organism]'

            # submit BLAST search
            result_handle = NCBIWWW.qblast(**blast_params)
            blast_records = NCBIXML.parse(result_handle)

            # parse results
            hits = []
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        # match percentage: considers NCBI match percentage for sequences
                        match_percent = (hsp.identities / hsp.align_length) * 100

                        # consider significant matches above 80%
                        if match_percent >= 80:  # 80% id
                            hits.append({
                                'title': alignment.title[:100],  # truncate long titles
                                'accession': alignment.accession,
                                'length': alignment.length,
                                'e_value': hsp.expect,
                                'identity': hsp.identities,
                                'align_length': hsp.align_length,
                                'match_percent': round(match_percent, 1),
                                'query_start': hsp.query_start,
                                'query_end': hsp.query_end,
                                'subject_start': hsp.sbjct_start,
                                'subject_end': hsp.sbjct_end
                            })

            result_handle.close()

            # match percentage and e-value
            hits.sort(key=lambda x: (-x['match_percent'], x['e_value']))

            # returns the top hits; not overpopulating
            return hits[:max_hits]

        # uh-oh... BLAST failed; the sequence has failed the fallback
        except Exception as e:
            print(f"NCBI BLAST search failed: {str(e)}")
            print("Using fallback specificity check.")
            return self._fallback_specificity_check(sequence)

    def _fallback_specificity_check(self, sequence):
        issues = []

        # checking for match % errors
        for common_seq in self.common_sequences:
            # forward direction
            for i in range(len(common_seq) - len(sequence) + 1):
                target = common_seq[i:i + len(sequence)]
                match_count = sum(1 for a, b in zip(sequence, target) if a == b)
                match_percent = (match_count / len(sequence)) * 100

                if match_percent >= 80:
                    issues.append({
                        'title': f'Common sequence match: {common_seq[:30]}...',
                        'match_percent': round(match_percent, 1),
                        'type': 'local_database'
                    })

            # reverse compliment
            seq_rc = self.reverse_complement(sequence)
            for i in range(len(common_seq) - len(seq_rc) + 1):
                target = common_seq[i:i + len(seq_rc)]
                match_count = sum(1 for a, b in zip(seq_rc, target) if a == b)
                match_percent = (match_count / len(seq_rc)) * 100

                if match_percent >= 80:
                    issues.append({
                        'title': f'Common sequence match (RC): {common_seq[:30]}...',
                        'match_percent': round(match_percent, 1),
                        'type': 'local_database'
                    })

        return issues
    # def for human seq. for RefSeq RNA db
    def check_primer_specificity_ncbi(self, sequence, organism='Homo sapiens',
                                    database='refseq_rna', max_hits=20):

        print(f"Checking primer specificity against NCBI {database}...")

        # BLAST search
        blast_hits = self.ncbi_blast_search(
            sequence=sequence,
            database=database,
            organism=organism,
            max_hits=max_hits,
            expect_threshold=1000,
            word_size=7
        )

        # analyze results for specificity issues
        specificity_issues = []

        for hit in blast_hits:
            # flag! high-identity matches that might cause non-specific amplification
            if hit.get('match_percent', 0) >= 85:
                specificity_issues.append({
                    'database': database,
                    'title': hit.get('title', 'Unknown'),
                    'accession': hit.get('accession', ''),
                    'match_percent': hit.get('match_percent', 0),
                    'e_value': hit.get('e_value', 0),
                    'align_length': hit.get('align_length', 0),
                    'risk_level': self._assess_specificity_risk(hit)
                })

        return {
            'total_hits': len(blast_hits),
            'high_risk_hits': len([h for h in specificity_issues if h['risk_level'] == 'HIGH']),
            'medium_risk_hits': len([h for h in specificity_issues if h['risk_level'] == 'MEDIUM']),
            'low_risk_hits': len([h for h in specificity_issues if h['risk_level'] == 'LOW']),
            'hits': specificity_issues[:10],  # Top 10 hits
            'overall_risk': self._calculate_overall_specificity_risk(specificity_issues)
        }

    def _assess_specificity_risk(self, hit):
        # risk for primer specificity on BLAST
        match_percent = hit.get('match_percent', 0)
        align_length = hit.get('align_length', 0)

        # high: >95% identity over primer length
        if match_percent >= 95 and align_length >= len(hit.get('query_sequence', '')) * 0.8:
            return 'HIGH'
        # medium: >90% identity over portion
        elif match_percent >= 90 and align_length >= len(hit.get('query_sequence', '')) * 0.6:
            return 'MEDIUM'
        # low: >85% identity but shorter alignment
        elif match_percent >= 85:
            return 'LOW'
        else:
            return 'NEGLIGIBLE'

    # overall specificity risk
    def _calculate_overall_specificity_risk(self, specificity_issues):
        if not specificity_issues:
            return 'LOW'

        high_risk_count = len([h for h in specificity_issues if h['risk_level'] == 'HIGH'])
        medium_risk_count = len([h for h in specificity_issues if h['risk_level'] == 'MEDIUM'])

        # spec. risk score
        if high_risk_count >= 3:
            return 'HIGH'
        elif high_risk_count >= 1 or medium_risk_count >= 5:
            return 'MEDIUM'
        else:
            return 'LOW'

    # rev compliment specificity risk
    def reverse_complement(self, seq):
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(comp.get(base.upper(), 'N') for base in reversed(seq))

    # HAIRPINS**
    # checks hairpin secondary structures
    def check_hairpin(self, seq, min_loop=3, min_stem=4):
        hairpins = []
        seq_len = len(seq)

        # for each sequence, check for compliments
        for i in range(seq_len - min_stem * 2 - min_loop):
            for j in range(i + min_stem + min_loop, seq_len - min_stem + 1):
                left_arm = seq[i:i + min_stem]
                right_arm = seq[j:j + min_stem]
                right_arm_rc = self.reverse_complement(right_arm)

                # analyzing hairpin structures
                matches = sum(1 for a, b in zip(left_arm, right_arm_rc) if a == b) #comp L/R arms
                if matches >= min_stem - 1:  # allow one mismatch
                    loop_size = j - i - min_stem  #gap
                    stability = matches * 2 - (min_stem - matches) * 0.5 #inc stab on match pairs, penalize mismatch
                    hairpins.append({
                        'start': i,
                        'end': j + min_stem,
                        'stem_length': min_stem,
                        'loop_size': loop_size,
                        'stability': stability,
                        'sequence': seq[i:j + min_stem]
                    })

        return hairpins

    def visualize_hairpin(self, seq, start, end):
        left = seq[:start]
        loop = seq[start:end]
        right = seq[end:]

        pointer_line = ' ' * len(left) + '↑' * len(loop)

        print("\n🧷 Detected Hairpin Structure:")
        print(f"5′ {seq} 3′")
        print(f"   {pointer_line}")
        print(f"   {' ' * len(left)}Loop ({len(loop)} nt)")

    # PRIMER & DIMER STABILITY AND RECOGNITION***
    # store stability
    def check_dimers(self, seq1, seq2=None, min_match=4, max_stability=None):
        if max_stability is None:
            max_stability = 5
        if seq2 is None:
            seq2 = seq1  # self-dimer

        dimers = []

        # Check all possible alignments
        for offset in range(-len(seq1) + min_match, len(seq2) - min_match + 1):
            matches = 0
            match_positions = []

        # check complimentary bases seq1/seq2
            for i in range(len(seq1)):
                j = i + offset
                if 0 <= j < len(seq2):
                    if seq1[i] == self.complement_map.get(seq2[j], ''):
                        matches += 1
                        match_positions.append((i, j))

            if matches >= min_match:
                # 3' end involvement
                end_matches = sum(1 for pos_i, pos_j in match_positions[-3:] if pos_i >= len(seq1)- 3 or pos_j >= len(seq2) - 3)

                stability = max(0, matches + end_matches * 2)
                if stability <= max_stability:

                    dimers.append({
                        'matches': matches,
                        'offset': offset,
                        'end_matches': end_matches,
                        'stability': matches + end_matches * 2,
                        'positions': match_positions
                    })

        return sorted(dimers, key=lambda x: x['stability'], reverse=True)

    # enhanced specificity with Entrez
    def check_specificity(self, seq, organism='Homo sapiens', database='refseq_rna', max_matches=2):
        # NCBI database check
        if BIOPYTHON_AVAILABLE and not self.fast_mode:
            return self.check_primer_specificity_ncbi(seq, organism, database)
        else:
            # fallback to local check
            issues = self._fallback_specificity_check(seq)
            return {
                'total_hits': len(issues),
                'high_risk_hits': len([h for h in issues if h.get('match_percent', 0) >= 95]),
                'medium_risk_hits': len([h for h in issues if 85 <= h.get('match_percent', 0) < 95]),
                'low_risk_hits': len([h for h in issues if h.get('match_percent', 0) < 85]),
                'hits': issues,
                'overall_risk': 'MEDIUM' if len(issues) > 3 else 'LOW'
            }

    # primer qc with Entrez
    def analyze_primer_quality(self, seq, organism='Homo sapiens', database='refseq_rna'):
        quality = {
            'sequence': seq,
            'length': len(seq),
            'tm': round(calc_tm_nn(seq, Default_NA, Default_MG, Default_dNTP), 1),
            'gc_content':  round((seq.count('G') + seq.count('C')) / len(seq) * 100, 1),
            'hairpins': self.check_hairpin(seq),
            'self_dimers': self.check_dimers(seq),
            'warnings': []

        }
        # specificity analysis: full NCBI
        if not self.fast_mode and BIOPYTHON_AVAILABLE:
            spec = self.check_primer_specificity_ncbi(seq, organism, database)
        else:
            issues = self._fallback_specificity_check(seq)
            spec = {
                'total_hits': len(issues),
                'high_risk_hits': len([h for h in issues if h['match_percent'] >= 95]),
                'medium_risk_hits': len([h for h in issues if 85 <= h['match_percent'] < 95]),
                'low_risk_hits': len([h for h in issues if h['match_percent'] < 85]),
                'hits': issues,
                'overall_risk': 'MEDIUM' if len(issues) > 3 else 'LOW'
            }

        quality['specificity_analysis'] = spec
        quality['specificity_issues']     = spec.get('hits', [])

        # warnings on analysis
        if quality['tm'] < Def_Melt_Tm_LOWER or quality['tm'] > Def_Melt_Tm_UPPER:
            quality['warnings'].append(f"Tm ({quality['tm']}°C) outside optimal range (55-65°C)")

        # warning on GC content
        if quality['gc_content'] < Def_GC_LOWER or quality['gc_content'] > Def_GC_UPPER:
            quality['warnings'].append(f"GC content ({quality['gc_content']}%) outside optimal range (40-60%)")

        # warning on hairpins
        if quality['hairpins']:
            strong_hairpins = [h for h in quality['hairpins'] if h['stability'] > Def_Hairpin_Min]
            if strong_hairpins:
                quality['warnings'].append(f"Strong hairpin structures detected ({len(strong_hairpins)})")
                for h in strong_hairpins:
                    print(
                        f"🧷 Hairpin start: {h['start']} | end: {h['end']} | loop: {h['loop_size']} nt | stability: {h['stability']:.2f}")
                    self.visualize_hairpin(seq, h['start'], h['end'])

        # warning on self-dimers
        if quality['self_dimers']:
            strong_dimers = [d for d in quality['self_dimers'] if d['stability'] > Def_SelfDimer_Min]
            if strong_dimers:
                quality['warnings'].append(f"Strong self-dimer potential detected")

        # specificity warnings
        spec_analysis = quality['specificity_analysis']
        if spec_analysis['overall_risk'] == 'HIGH':
            quality['warnings'].append(f"HIGH specificity risk - {spec_analysis['high_risk_hits']} high-risk matches")
        elif spec_analysis['overall_risk'] == 'MEDIUM':
            quality['warnings'].append(f"MEDIUM specificity risk - check target specificity")

        # calculate overall qcs with enhanced scoring
        score = 100
        score -= len(quality['warnings']) * 8
        score -= len(quality['hairpins']) * 2
        score -= len([d for d in quality['self_dimers'] if d['stability'] > 3]) * 5

        # specificity penalty
        if spec_analysis['overall_risk'] == 'HIGH':
            score -= 30
        elif spec_analysis['overall_risk'] == 'MEDIUM':
            score -= 15

        # scoring guidelines per spec analysis
        score -= spec_analysis['high_risk_hits'] * 10
        score -= spec_analysis['medium_risk_hits'] * 5

        quality['quality_score'] = max(0, score)

        return quality

    # optimizing primer pairs based on dimer formation
    def optimize_primer_pair(self, forward_primers, reverse_primers):
        best_pairs = []
        for f_primer in forward_primers:
            for r_primer in reverse_primers:
                # hetero-dimers
                hetero_dimers = self.check_dimers(f_primer['sequence'], r_primer['sequence'])

                # pair score
                tm_diff = abs(f_primer['tm'] - r_primer['tm'])
                pair_score = (f_primer.get('quality_score', 50) + r_primer.get('quality_score', 50)) / 2

                # large Tm differences and strong hetero-dimers
                pair_score -= tm_diff * 5
                if hetero_dimers:
                    strongest_dimer = max(hetero_dimers, key=lambda x: x['stability'])
                    if strongest_dimer['stability'] > 5:
                        pair_score -= strongest_dimer['stability'] * 10

                # store best pairs and information
                best_pairs.append({
                    'forward': f_primer,
                    'reverse': r_primer,
                    'tm_difference': round(tm_diff, 1),
                    'hetero_dimers': len([d for d in hetero_dimers if d['stability'] > 3]),
                    'pair_score': max(0, pair_score)
                })

        return sorted(best_pairs, key=lambda x: x['pair_score'], reverse=True)

    # AA to NA
    def is_amino_acid_sequence(self, sequence):
        aa_only = set('DEFHIKLMNPQRSVWY')
        return any(char in aa_only for char in sequence.upper())

    # AA to NA
    def amino_acid_to_dna(self, aa_seq):
        codon_table = {
            'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
            'Q': 'CAG', 'E': 'GAG', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
            'L': 'CTG', 'K': 'AAG', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
            'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT'
        }
        dna_seq = ''
        for aa in aa_seq.upper():
            if aa in codon_table:
                dna_seq += codon_table[aa]
            elif aa not in 'ATGC':
                continue
        return dna_seq

    # find primers based on query
    def find_primers(self, sequence, target_tm=60, product_size=500, max_primers=10,
                    organism='Homo sapiens', database='refseq_rna'):
        # convert AA seq
        if self.is_amino_acid_sequence(sequence):
            print("Amino acid sequence detected - converting to DNA...")
            sequence = self.amino_acid_to_dna(sequence)
            print(f"Converted DNA sequence: {sequence[:50]}...")

        sequence = re.sub(r'[^ATGC]', '', sequence.upper())
        if len(sequence) < 100:
            return None

        print(f"Analyzing primers with NCBI {database} database for {organism}...")

        forward_primers = []
        reverse_primers = []

        # forward primers with enhanced qc
        print("Designing forward primers...")
        for length in range(18, 26):
            for start in range(min(100, len(sequence) - length)):
                primer = sequence[start:start + length]
                quality = self.analyze_primer_quality(primer, organism, database)

                # analyzing qc
                if Def_Melt_Tm_LOWER <= quality['tm'] <= Def_Melt_Tm_UPPER and Def_GC_LOWER <= quality['gc_content'] <= Def_GC_UPPER:
                    quality.update({
                        'position': f"{start+1}-{start+length}",
                        'direction': 'forward'
                    })
                    forward_primers.append(quality)

                    if len(forward_primers) >= 3:
                        break

        # reverse primers with enhanced qc
        print("Designing reverse primers...")
        rev_start = max(0, len(sequence) - 200)
        rev_end = len(sequence)

        for length in range(18, 26):
            for start in range(rev_start, rev_end - length):
                if start + length > len(sequence):
                    continue
                primer_seq = sequence[start:start + length]
                primer = self.reverse_complement(primer_seq)
                quality = self.analyze_primer_quality(primer, organism, database)

                if Def_Melt_Tm_LOWER <= quality['tm'] <= Def_Melt_Tm_UPPER and Def_GC_LOWER <= quality['gc_content'] <= Def_GC_UPPER:
                    quality.update({
                        'position': f"{start+1}-{start+length}",
                        'direction': 'reverse'
                    })
                    reverse_primers.append(quality)

        # quality score and target Tm
        forward_primers.sort(key=lambda x: (x['quality_score'], -abs(x['tm'] - target_tm)), reverse=True)
        reverse_primers.sort(key=lambda x: (x['quality_score'], -abs(x['tm'] - target_tm)), reverse=True)

        # primer pairs
        if forward_primers and reverse_primers:
            print("Optimizing primer pairs...")
            top_forward = forward_primers[:min(5, len(forward_primers))]
            top_reverse = reverse_primers[:min(5, len(reverse_primers))]
            optimized_pairs = self.optimize_primer_pair(top_forward, top_reverse)

            return {
                'forward': forward_primers[:max_primers],
                'reverse': reverse_primers[:max_primers],
                'optimized_pairs': optimized_pairs[:5]
            }

        return {
            'forward': forward_primers[:max_primers],
            'reverse': reverse_primers[:max_primers],
            'optimized_pairs': []
        }

    # main output framework
    def design_primers(self, sequence, detailed=True, verbose=False):
        print("ENHANCED PRIMER DESIGNER -- Designed by Joshua Jenkelowitz")
        print(f"Primer3-py available: {PRIMER3_AVAILABLE}")
        print(f"BIOPYTHON available: {BIOPYTHON_AVAILABLE}")
        print()

        from src.primer3_wrapper import pick_with_primer3
        from src.qc import filter_primer3_results
        from primer_config import Default_NA, Default_MG

        # to primer3_wrapper.py
        raw = pick_with_primer3(sequence, Default_NA, Default_MG)
        pairs_raw = filter_primer3_results(raw)
        if not pairs_raw:
            print("No suitable primers found")
            return {"forward": [], "reverse": [], "optimized_pairs": []}

        # pairs storage
        pairs = []
        for pr in pairs_raw[:5]:
            f = pr['forward']
            r = pr['reverse']
            f_dict = {
                'sequence': f['SEQUENCE'],
                'tm': round(f['Tm'], 1),
                'gc_content': round(f['GC'], 1),
                'position': '',
                'direction': 'forward',
                'warnings': []
            }
            r_dict = {
                'sequence': r['SEQUENCE'],
                'tm': round(r['Tm'], 1),
                'gc_content': round(r['GC'], 1),
                'position': '',
                'direction': 'reverse',
                'warnings': []
            }

            # dimer storage
            hetero_dimers = self.check_dimers(f_dict['sequence'], r_dict['sequence'])
            strong_dimer_count = len([d for d in hetero_dimers if d['stability'] > 3])

            # map fwd primer
            f_index = sequence.find(f_dict['sequence'])
            if f_index != -1:
                f_dict['position'] = f"{f_index + 1}-{f_index + len(f_dict['sequence'])}"

            # map rev primer
            r_seq_rc = self.reverse_complement(r_dict['sequence'])
            r_index = sequence.find(r_seq_rc)
            if r_index != -1:
                r_dict['position'] = f"{r_index + 1}-{r_index + len(r_dict['sequence'])}"

            # tm difference and scoring
            tm_diff = abs(f_dict['tm'] - r_dict['tm'])
            pair_score = 100.0 - tm_diff * 2.5 - strong_dimer_count * 5.5
            pairs.append({
                'forward': f_dict,
                'reverse': r_dict,
                'tm_difference': round(tm_diff, 1),
                'pair_score': max(0, round(pair_score, 1)),
                'hetero_dimers': strong_dimer_count
            })

        # storage accession
        forward_primers = [p['forward'] for p in pairs]
        reverse_primers = [p['reverse'] for p in pairs]
        optimized_pairs = pairs

        print(f"Forward primers found: {len(forward_primers)}")
        print(f"Reverse primers found: {len(reverse_primers)}")

        from src.primer_output import print_optimized_primers, visualize_amplicon
        top_pairs = optimized_pairs[:5]

        # index top pairs for top 5 results
        if top_pairs:
            print_optimized_primers(top_pairs)
            for idx, pair in enumerate(top_pairs, 1):
                f_pos = pair['forward']['position']
                r_pos = pair['reverse']['position']
                f_start, f_end = map(int, f_pos.split('-')) if '-' in f_pos else (None, None)
                r_start, r_end = map(int, r_pos.split('-')) if '-' in r_pos else (None, None)

                # visualize amplicon sequence
                if f_start and r_start:
                    f_start -= 1
                    r_start -= 1
                    visualize_amplicon(sequence, f_start, f_end, r_start, r_end, wrap=60, idx=idx)

        # output framework
        rows = []
        for p in forward_primers + reverse_primers:
            rows.append({
                'Direction': p['direction'],
                'Sequence': p['sequence'],
                'Quality Score': p.get('quality_score', ''),
                'Tm (°C)': p['tm'],
                'GC (%)': p['gc_content'],
                'Position': p['position'],
                'Warnings': ";".join(p['warnings']) if p['warnings'] else ""
            })
        df = pd.DataFrame(rows)

        # export_utils.py followup code... export to excel
        from export_utils import export_top_primers_to_excel
        response = input("\nExport all primers, top pairs only, or skip export? (all / top / none): ").strip().lower()
        if response == 'all':
            export_top_primers_to_excel(
                df=df,
                full_sequence=sequence,
                output_folder="output"
            )
        elif response == 'top':
            top_df = pd.DataFrame([
                {
                    'Forward Sequence': pair['forward']['sequence'],
                    'Reverse Sequence': pair['reverse']['sequence'],
                    'Forward Tm (°C)': pair['forward']['tm'],
                    'Reverse Tm (°C)': pair['reverse']['tm'],
                    'Forward GC (%)': pair['forward']['gc_content'],
                    'Reverse GC (%)': pair['reverse']['gc_content'],
                    'Forward Pos': pair['forward']['position'],
                    'Reverse Pos': pair['reverse']['position'],
                    'Pair Tm Δ (°C)': pair['tm_difference'],
                    'Pair Score': pair['pair_score'],
                    'QC Rank': (
                        "Ideal" if pair['pair_score'] >= 85
                        else "Moderate" if pair['pair_score'] >= 70
                        else "Risky"
                    )
                }
                for pair in top_pairs
            ])
            float_cols = [col for col in top_df.columns if "Tm" in col or "GC" in col or "Δ" in col or "Score" in col]
            top_df[float_cols] = top_df[float_cols].round(1)
            export_top_primers_to_excel(
                df=top_df,
                full_sequence=sequence,
                output_folder="output"
            )
        else:
            print("No Excel file will be created.")

        return {"forward": forward_primers, "reverse": reverse_primers, "optimized_pairs": optimized_pairs}

