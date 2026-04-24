import subprocess

res1 = subprocess.run(["C:\\R\\R-4.3.3\\bin\\Rscript.exe", "c:/ProjectTEDGWAS/TrackA_MR/scripts/validate_sensitivity.R"], capture_output=True)
res2 = subprocess.run(["python", "c:/ProjectTEDGWAS/TrackC_Offtarget/scripts/rebuild_fig4_enrichr.py"], capture_output=True)

with open("c:/ProjectTEDGWAS/Validation_Full.txt", "wb") as f:
    f.write(b"=== MR SENSITIVITY ===\n")
    f.write(res1.stdout)
    f.write(b"\n=== ENRICHR & FIG 4 ===\n")
    f.write(res2.stdout)
