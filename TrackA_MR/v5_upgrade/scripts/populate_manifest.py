import csv
import hashlib
import os
import datetime

manifest_in = r"c:\ProjectTEDGWAS\Week1_file_manifest_v1.csv"
manifest_out = r"c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\Week1_file_manifest_v1.csv"
base_dir = r"c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade"

with open(manifest_in, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    rows = list(reader)

for row in rows:
    file_path = row['file_path']
    abs_path = os.path.join(base_dir, file_path.replace('/', '\\'))
        
    if os.path.exists(abs_path):
        # file size
        row['file_size_bytes'] = str(os.path.getsize(abs_path))
        # md5
        md5 = hashlib.md5()
        with open(abs_path, 'rb') as f_in:
            for chunk in iter(lambda: f_in.read(4096), b""):
                md5.update(chunk)
        row['md5_hash'] = md5.hexdigest()
        # created at
        mtime = os.path.getmtime(abs_path)
        row['created_at'] = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')
    else:
        if 'TBD' in row['md5_hash']:
            row['md5_hash'] = 'FILE_NOT_FOUND'

with open(manifest_out, 'w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
print("Manifest populated successfully!")
