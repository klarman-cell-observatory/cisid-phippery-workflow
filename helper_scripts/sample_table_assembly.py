import pandas as pd
from google.cloud import storage

BUCKET_NAME     = "fc-2286e7bb-9ae4-4814-a0d4-28fa504d3d59"
BASEPATH        = "PhIP_Seq/raw/fastq_files/250515_VH00997_413_AACVG7HHV_fastq/250515_VH00997_413_AACVG7HHV_fastqs/sample_fastqs"
SAMPLE_TABLE_IN = "/Users/chene/Downloads/Sepsis_PhIPseq_batch3_sample_table.csv"

storage_client = storage.Client()
bucket = storage_client.bucket(BUCKET_NAME)

all_blobs = list(bucket.list_blobs(prefix=BASEPATH))

r1_uris = [f"gs://{BUCKET_NAME}/{b.name}" for b in all_blobs if "_R1_" in b.name and b.name.endswith(".fastq.gz")]
r2_uris = [f"gs://{BUCKET_NAME}/{b.name}" for b in all_blobs if "_R2_" in b.name and b.name.endswith(".fastq.gz")]

print(f"Found {len(r1_uris)} R1 files and {len(r2_uris)} R2 files under {BASEPATH}")


def find_uri(sample_id, uri_list):
    """Return the GCS URI whose filename contains /{sample_id}_."""
    return next((u for u in uri_list if f"/{sample_id}_" in u), None)


sample_table = pd.read_csv(SAMPLE_TABLE_IN)

sample_table["R1_cloud_filepath"] = sample_table["sample_ID"].apply(lambda sid: find_uri(sid, r1_uris))
sample_table["R2_cloud_filepath"] = sample_table["sample_ID"].apply(lambda sid: find_uri(sid, r2_uris))

# Drop columns that the new WDL no longer uses
sample_table = sample_table.drop(columns=["cloud_filepath", "fastq_filepath"], errors="ignore")

# Report any samples that couldn't be matched
missing = sample_table[sample_table["R1_cloud_filepath"].isna() | sample_table["R2_cloud_filepath"].isna()]
if not missing.empty:
    print("WARNING: could not find R1/R2 for the following samples:")
    for _, row in missing.iterrows():
        print(f"  {row['sample_ID']}  R1={row['R1_cloud_filepath']}  R2={row['R2_cloud_filepath']}")

sample_table.to_csv(SAMPLE_TABLE_IN, index=False)
print(f"Wrote {len(sample_table)} rows to {SAMPLE_TABLE_IN}")
