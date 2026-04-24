import gzip, csv, os

finngen_gz = "c:/ProjectTEDGWAS/TrackA_MR/data/finngen_R12_GRAVES_OPHT.gz"
output_csv = "c:/ProjectTEDGWAS/TrackA_MR/data/finngen_TSHR_region.csv"

if not os.path.exists(finngen_gz):
    print("FinnGen GZ file not found.")
    exit(1)

print(f"Extracting TSHR region (chr14:80.5M-81.7M) from {finngen_gz}...")
with gzip.open(finngen_gz, "rt", encoding="utf-8") as f_in, \
     open(output_csv, "w", newline="", encoding="utf-8") as f_out:
    
    reader = csv.reader(f_in, delimiter="\t")
    writer = csv.writer(f_out, delimiter=",")
    
    header = next(reader)
    writer.writerow(header)
    
    count = 0
    try:
        for row in reader:
            if not row: continue
            chrom = row[0].replace('chr', '')
            try:
                pos = int(row[1])
            except ValueError:
                continue
                
            if chrom == "14" and 80500000 <= pos <= 81700000:
                writer.writerow(row)
                count += 1
    except Exception as e:
        print(f"Read error: {e}")

print(f"Extracted {count} SNPs to {output_csv}")
