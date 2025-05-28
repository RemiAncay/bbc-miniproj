from pathlib import Path
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


genetic_code = {
    'ATG': 'M',
    'TTT': 'F', 'TTC': 'F',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y',
    'CAT': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C',
    'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'TAA': '*', 'TAG': '*', 'TGA': '*'
}

def translate_dna_to_protein(dna_sequence: str, start_codon: str = 'ATG') -> str:
    """
    Traduit une séquence ADN en protéine à partir du premier codon start ATG
    jusqu’au premier codon stop (non inclus). Les codons inconnus deviennent 'X'.
    """
    dna_sequence = dna_sequence.upper()
    start_index = dna_sequence.find(start_codon)
    if start_index == -1:
        return ''                                # pas de start → pas de protéine
    
    protein = []
    for i in range(start_index, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) < 3:
            break                               # codon incomplet en fin de chaîne
        #if codon not in ['TAA', 'TAG', 'TGA']:      # codon stop
        protein.append(genetic_code.get(codon, ''))
    return ''.join(protein)


def translate_fasta_file(nucl_fasta: Path) -> Path:
    """
    Traduit toutes les séquences d'un fichier FASTA nucléotidique.
    Renvoie le chemin du fichier FASTA protéique produit.
    """
    output_path = nucl_fasta.with_name(nucl_fasta.stem + "_Protein.fasta")
    translated_records = []

    for record in SeqIO.parse(nucl_fasta, "fasta"):
        prot_seq = translate_dna_to_protein(str(record.seq))
        # On conserve l'identifiant et la description
        translated_records.append(
            SeqRecord(Seq(prot_seq),
                      id=record.id,
                      description="translated_from_" + record.id)
        )

    SeqIO.write(translated_records, output_path, "fasta")
    return output_path


if __name__ == "__main__":
 
    fasta_path = Path("data/assemblies.fasta")

    out_path = translate_fasta_file(fasta_path)
    print(f"Traduction terminée : {out_path}")
