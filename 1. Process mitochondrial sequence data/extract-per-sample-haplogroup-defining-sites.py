import pandas as pd
import re

df = pd.read_csv('haplogroups_extended.csv', sep='\t')
print(f"Loaded {len(df)} samples")

def parse_positions(poly_string):
    if pd.isna(poly_string) or str(poly_string).strip() == '':
        return []
    positions = []
    for token in poly_string.split():
        clean = re.sub(r'\(.*?\)', '', token).strip()
        m = re.match(r'^(\d+)[A-Za-z!d]+$', clean)
        if m:
            positions.append(int(m.group(1)))
    return positions

records = []
for _, row in df.iterrows():
    for pos in parse_positions(row['Found_Polys']):
        records.append({'SampleID': row['SampleID'], 'Position': pos})

haplo_df = pd.DataFrame(records)
haplo_df.to_csv('per_sample_haplodefining_positions.tsv', sep='\t', index=False)
print(f"Written: {len(haplo_df)} rows")
print(f"Mean haplogroup-defining sites per individual: {haplo_df.groupby('SampleID').size().mean():.1f}")