# PowerShell script to download official MAGMA 1000G hg19 panels from SurfSara

$data_dir = "c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\data"
if (!(Test-Path $data_dir)) {
    New-Item -ItemType Directory -Force -Path $data_dir
}

Write-Output "Canceling all active BITS transfers..."
Get-BitsTransfer | Remove-BitsTransfer -ErrorAction SilentlyContinue

# Clean up any existing partial files
Remove-Item -Path "$data_dir\g1000_eur.zip" -ErrorAction SilentlyContinue
Remove-Item -Path "$data_dir\g1000_eas.zip" -ErrorAction SilentlyContinue

Write-Output "Starting background BITS transfer for official EUR panel (SurfSara)..."
Start-BitsTransfer -Source "https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download" -Destination "$data_dir\g1000_eur.zip" -Asynchronous

Write-Output "Starting background BITS transfer for official EAS panel (SurfSara)..."
Start-BitsTransfer -Source "https://vu.data.surfsara.nl/index.php/s/dz6PYdKOi3xVqHn/download" -Destination "$data_dir\g1000_eas.zip" -Asynchronous

# List active BITS transfers
Get-BitsTransfer
