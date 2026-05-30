# PowerShell script to start BITS transfers for 1000G panels

$data_dir = "c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\data"
if (!(Test-Path $data_dir)) {
    New-Item -ItemType Directory -Force -Path $data_dir
}

# Remove partial files
Remove-Item -Path "$data_dir\g1000_eur.zip" -ErrorAction SilentlyContinue
Remove-Item -Path "$data_dir\g1000_eas.zip" -ErrorAction SilentlyContinue

Write-Output "Starting background BITS transfer for EUR panel..."
Start-BitsTransfer -Source "https://zenodo.org/records/14142263/files/g1000_eur.zip?download=1" -Destination "$data_dir\g1000_eur.zip" -Asynchronous

Write-Output "Starting background BITS transfer for EAS panel..."
Start-BitsTransfer -Source "https://zenodo.org/records/14108873/files/g1000_eas.zip?download=1" -Destination "$data_dir\g1000_eas.zip" -Asynchronous

# List active BITS transfers
Get-BitsTransfer
