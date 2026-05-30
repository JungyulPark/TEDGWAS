# PowerShell script to download and extract MAGMA reference panels from Zenodo

$data_dir = "c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\data"
New-Item -ItemType Directory -Force -Path $data_dir

# 1. Download and Extract EUR Reference
$eur_zip = "$data_dir\g1000_eur.zip"
$eur_dest = "$data_dir\1000G_EUR_phase3"
Write-Output "Downloading EUR Reference from Zenodo..."
curl.exe -L -o $eur_zip "https://zenodo.org/records/14142263/files/g1000_eur.zip?download=1"

Write-Output "Extracting EUR Reference..."
New-Item -ItemType Directory -Force -Path $eur_dest
Expand-Archive -Path $eur_zip -DestinationPath $eur_dest -Force
Remove-Item -Path $eur_zip
Write-Output "EUR Reference installation complete."

# 2. Download and Extract EAS Reference
$eas_zip = "$data_dir\g1000_eas.zip"
$eas_dest = "$data_dir\1000G_EAS_phase3"
Write-Output "Downloading EAS Reference from Zenodo..."
curl.exe -L -o $eas_zip "https://zenodo.org/records/14108873/files/g1000_eas.zip?download=1"

Write-Output "Extracting EAS Reference..."
New-Item -ItemType Directory -Force -Path $eas_dest
Expand-Archive -Path $eas_zip -DestinationPath $eas_dest -Force
Remove-Item -Path $eas_zip
Write-Output "EAS Reference installation complete."

Write-Output "All panels downloaded and extracted successfully."
