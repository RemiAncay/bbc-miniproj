import os
import sys
import glob
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def assemble_sequences_in_file(fasta_path: Path) -> SeqRecord | None:
    """Renvoie un SeqRecord résultant de la concaténation
       de toutes les séquences d’un fichier FASTA."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:                       # fichier vide
        return None

    merged_seq = Seq("".join(str(r.seq) for r in records))
    # Identifiant = nom de fichier sans extension
    seq_id = fasta_path.stem
    # Optionnel : conservez le 1er descripteur si vous préférez
    description = f"assembled_from_{fasta_path.name}"
    return SeqRecord(merged_seq, id=seq_id, description=description)

def main(folder: Path) -> None:
    # Extensions considérées ; adaptez si besoin
    patterns = ["*.fa", "*.fasta", "*.fna", "*.fas"]
    fasta_files = [p for pattern in patterns for p in folder.glob(pattern)]

    if not fasta_files:
        sys.exit(f"Aucun fichier FASTA trouvé dans {folder}")

    assembled_records = []
    for f in sorted(fasta_files):         # tri pour reproductibilité
        rec = assemble_sequences_in_file(f)
        if rec:
            assembled_records.append(rec)
            print(f"{f.name:30}  ->  ajouté ({len(rec.seq):,} nt)")
        else:
            print(f"{f.name:30}  ->  ignoré (vide)")

    output_path = folder / "assemblies.fasta"
    SeqIO.write(assembled_records, output_path, "fasta")
    print(f"Assemblage terminé : {len(assembled_records)} séquences écrites dans {output_path}")

if __name__ == "__main__":
    target_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data")
    main(target_dir)
