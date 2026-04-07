import sys
import pandas as pd
from Bio import SeqIO # Requires biopython in your environment

fasta_in = sys.argv[1]
csv_out = sys.argv[2]

metadata = []

# Parse the UniProt FASTA downloaded in our earlier test rule
for record in SeqIO.parse(fasta_in, "fasta"):
    # We assign 'True' to is_characterized because we specifically downloaded 
    # 'reviewed' (Swiss-Prot) lectins which are all known ("light") proteins.
    metadata.append({
        "protein_id": record.id,
        "is_characterized": True, 
        "has_interpro": True,
        "description": record.description
    })

df = pd.DataFrame(metadata)
df.to_csv(csv_out, index=False)