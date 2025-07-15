import pandas as pd

# file loader from .xlsx: BATCH POSSIBLE
def load_sequences_from_excel(file_path):
    df = pd.read_excel(file_path)
    pathogen_sequences = {}
    ss_to_ls = {}
    gene_labels = {}

    # ls or ss conversion
    for _, row in df.iterrows():
        long_name = str(row.get("Pathogen Name", "")).strip()
        short_code = str(row.get("SS", "")).strip()
        gene = str(row.get("Target Gene", "")).strip()
        sequence = str(row.get("Upper Sequence", "")).strip().upper()

        if sequence:
            pathogen_sequences[long_name] = sequence
            gene_labels[long_name] = gene
            if short_code:
                ss_to_ls[short_code] = long_name

    return pathogen_sequences, ss_to_ls, gene_labels
