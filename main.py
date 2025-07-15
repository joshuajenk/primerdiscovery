import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from src.sequence_loader import load_sequences_from_excel
from src.primer_designer import PrimerDesigner

def main():
    filepath = "data/blood_borne_pathogens.xlsx" #change this based on desired xlsx file
    pathogen_sequences, ss_to_ls, _ = load_sequences_from_excel(filepath)
    designer = PrimerDesigner(
        pathogen_sequences=pathogen_sequences,
        ss_to_ls=ss_to_ls
    )
    designer.run_cli()

# cmd line interface run
def run_cli(self):
    print("\n Available pathogens:")
    for name in sorted(self.pathogen_sequences):
        print(f" - {name}")

    # prompt
    choice = input("\nEnter a pathogen name or shorthand code: ").strip()
    full_name = self.ss_to_ls.get(choice, choice)

    # y/n seq
    sequence = self.pathogen_sequences.get(full_name)
    if not sequence:
        print(f"No sequence found for '{choice}'.")
        return

    # parameters for target
    print(f"\nSequence for {full_name} loaded.")
    print(f"Length: {len(sequence)} bp")
    print(f"GC% Content: {gc_content(sequence)}%")
    print(f"Melt Temperature: {melt_tm(sequence)} °C")
    self.design_primers(sequence)

if __name__ == "__main__":
    main()

