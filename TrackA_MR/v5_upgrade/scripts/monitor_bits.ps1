# PowerShell script to monitor BITS downloads, extract them, and run PLINK pre-processing

$data_dir = "c:\ProjectTEDGWAS\TrackA_MR\v5_upgrade\data"
$plink_bin = "c:\ProjectTEDGWAS\tools\plink.exe"

# 1. Poll BITS transfers until complete
Write-Output "Monitoring active BITS jobs..."
while ($true) {
    $jobs = Get-BitsTransfer | Where-Object { $_.JobState -ne "Acknowledged" }
    if ($jobs.Count -eq 0) {
        Write-Output "No active BITS jobs found."
        break
    }
    
    $all_done = $true
    foreach ($job in $jobs) {
        if ($job.JobState -eq "Transferred") {
            $job | Complete-BitsTransfer
            Write-Output "BITS Job completed and finalized: $($job.DisplayName) (ID: $($job.JobId))"
        } elseif ($job.JobState -eq "Error") {
            Write-Output "ERROR: BITS Job $($job.JobId) failed. Details: $($job.Description)"
            # Try to resume or restart if failed
            $job | Resume-BitsTransfer -ErrorAction SilentlyContinue
            $all_done = $false
        } else {
            $all_done = $false
            $pct = if ($job.BytesTotal -gt 0) { [math]::Round(($job.BytesTransferred / $job.BytesTotal) * 100, 1) } else { 0 }
            Write-Output "Job $($job.JobId) ($($job.DisplayName)): $pct% ($($job.BytesTransferred) / $($job.BytesTotal) bytes) - Status: $($job.JobState)"
        }
    }
    
    if ($all_done) {
        Write-Output "All transfers completed successfully!"
        break
    }
    
    Start-Sleep -Seconds 15
}

# 2. Extract EUR reference if zip exists
$eur_zip = "$data_dir\g1000_eur.zip"
$eur_dest = "$data_dir\1000G_EUR_phase3"
if (Test-Path $eur_zip) {
    Write-Output "Extracting EUR Reference Zip..."
    New-Item -ItemType Directory -Force -Path $eur_dest
    Expand-Archive -Path $eur_zip -DestinationPath $eur_dest -Force
    Remove-Item -Path $eur_zip
    Write-Output "EUR Reference extracted."
}

# 3. Extract EAS reference if zip exists
$eas_zip = "$data_dir\g1000_eas.zip"
$eas_dest = "$data_dir\1000G_EAS_phase3"
if (Test-Path $eas_zip) {
    Write-Output "Extracting EAS Reference Zip..."
    New-Item -ItemType Directory -Force -Path $eas_dest
    Expand-Archive -Path $eas_zip -DestinationPath $eas_dest -Force
    Remove-Item -Path $eas_zip
    Write-Output "EAS Reference extracted."
}

# 4. PLINK Pre-processing: Extract Chromosome 14 and calculate frequencies
Write-Output "Running PLINK pre-processing for EUR..."
if (Test-Path "$eur_dest\g1000_eur.bed") {
    & $plink_bin --bfile "$eur_dest\g1000_eur" --chr 14 --make-bed --out "$eur_dest\EUR_chr14"
    & $plink_bin --bfile "$eur_dest\EUR_chr14" --freq --out "$eur_dest\EUR_chr14_freq"
    Write-Output "EUR PLINK pre-processing complete."
} else {
    Write-Output "WARNING: EUR g1000_eur.bed not found at $eur_dest"
}

Write-Output "Running PLINK pre-processing for EAS..."
if (Test-Path "$eas_dest\g1000_eas.bed") {
    & $plink_bin --bfile "$eas_dest\g1000_eas" --chr 14 --make-bed --out "$eas_dest\EAS_chr14"
    & $plink_bin --bfile "$eas_dest\EAS_chr14" --freq --out "$eas_dest\EAS_chr14_freq"
    Write-Output "EAS PLINK pre-processing complete."
} else {
    Write-Output "WARNING: EAS g1000_eas.bed not found at $eas_dest"
}

Write-Output "Pre-processing and setup finished!"
