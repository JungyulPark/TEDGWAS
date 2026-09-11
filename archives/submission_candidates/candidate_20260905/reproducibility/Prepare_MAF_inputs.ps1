param(
    [string]$OutputDirectory = (Join-Path $PSScriptRoot 'maf_inputs'),
    [string]$Rscript = 'C:\R\R-4.3.3\bin\Rscript.exe'
)
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
$taskOutput = (Resolve-Path -LiteralPath $OutputDirectory).Path
$afName = '2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz'
$afPath = Join-Path $taskOutput $afName
$afPartial = $afPath + '.partial'
$afUrl = 'https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/' + $afName
Write-Output 'This prepares inputs only; it does not perform a sensitivity analysis.'
if (-not (Test-Path -LiteralPath $afPath)) {
    Invoke-WebRequest -Uri $afUrl -OutFile $afPartial -UseBasicParsing -ErrorAction Stop
    $stream = [System.IO.File]::OpenRead($afPartial)
    try {
        if (($stream.ReadByte() -ne 31) -or ($stream.ReadByte() -ne 139)) {
            throw 'The response is not a gzip file. Do not use it as allele-frequency data.'
        }
    } finally { $stream.Dispose() }
    Move-Item -LiteralPath $afPartial -Destination $afPath
}
Get-FileHash -LiteralPath $afPath -Algorithm SHA256 | Format-List | Out-String | Set-Content -LiteralPath (Join-Path $taskOutput 'AF_SHA256.txt') -Encoding UTF8
if (Test-Path -LiteralPath $Rscript) {
    $probePath = Join-Path $taskOutput 'check_R_packages.R'
    @'
cat(R.version.string, "\n")
for (p in c("data.table", "TwoSampleMR", "coloc")) {
  ok <- requireNamespace(p, quietly=TRUE)
  cat(p, if (ok) as.character(packageVersion(p)) else "MISSING", "\n")
}
'@ | Set-Content -LiteralPath $probePath -Encoding UTF8
    & $Rscript --vanilla $probePath 2>&1 | Out-String | Set-Content -LiteralPath (Join-Path $taskOutput 'R_package_check.txt') -Encoding UTF8
}
Write-Output "Prepared files in: $taskOutput"
Write-Output 'If HTTPS certificate validation fails, stop and use the official eQTLGen download page. Do not bypass certificate checks.'
