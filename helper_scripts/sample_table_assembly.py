from google.cloud import storage

def token_match(sid, path_list):
    return next((p for p in path_list if f'/{sid}_' in p), None)

storage_client= storage.Client()
bucket= storage_client.bucket("fc-2286e7bb-9ae4-4814-a0d4-28fa504d3d59")
basepath = "PhIP_Seq/raw/fastq_files/250515_VH00997_413_AACVG7HHV_fastq/250515_VH00997_413_AACVG7HHV_fastqs/sample_fastqs"
blobs= bucket.list_blobs(prefix=basepath)

blob_list=[]
for i in blobs:
    blob_list.append(i)

reformat_blob_list= [str(x)[7:(len(str(x))-19)].replace(", ","/") for x in blob_list if "_R1" in str(x)]

def token_match(sid, path_list):
    return next((p for p in path_list if f'/{sid}_' in p), None)

sample_table = pd.read_csv("/Users/chene/Downloads/Sepsis_PhIPseq_batch3_sample_table.csv")

sample_table['cloud_filepath'] = sample_table['sample_ID'].apply(
    lambda sid: "gs://" + token_match(sid, reformat_blob_list)
)
sample_table['fastq_filepath'] = sample_table['cloud_filepath'].apply(lambda fq: "data/seq/" + fq.split("/")[-1])

sample_table.to_csv("/Users/chene/Downloads/Sepsis_PhIPseq_batch3_sample_table.csv")