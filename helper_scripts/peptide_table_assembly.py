import re
import pandas

with open("./new_library.fasta") as f:
    peptides = f.read()
    peptides = peptides.split("\n")

peptide_name = []
peptide_oligo = []
for i in range(0,len(peptides)):
    if i % 2 == 0:
        peptide_name.append(peptides[i].replace(">", ""))
    else:
        peptide_oligo.append(peptides[i])

peptide_df = pd.DataFrame(zip(peptide_name, peptide_oligo), columns=["name", "oligo"])

peptide_df['name_breakdown'] = peptide_df['name'].str.split("|")
peptide_df['name_breakdown_count'] = peptide_df['name_breakdown'].apply(len)

# Define the regex pattern
pattern = re.compile(r'^\d+-\d+$')

# Apply the match function to each list in the column
def extract_range_or_second(lst):
    for item in lst:
        if len(lst) <= 2:
            if pattern.match(item):
                return item
        elif len(lst) >= 3:
            if pattern.match(item):
                return item
            elif "CTERM" in item:
                return item
    return lst[1] if len(lst) > 1 else None

peptide_df['protein_range'] = peptide_df['name_breakdown'].apply(extract_range_or_second)

peptide_df['full_name'] = peptide_df['name']
peptide_df['gene'] = peptide_df['name_breakdown'].str[0]

#Peptide Table reorg
peptide_df.drop(columns=['name_breakdown', 'name_breakdown_count'])
peptide_df = peptide_df.loc[:, ['gene', 'full_name', 'oligo', 'protein_range']]

peptide_df.to_csv("peptide_table.csv", index_label='peptide_id')