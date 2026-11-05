# Bookkeeping for the Windows CI cache of conda package archives (ADR-0078), run from the
# workspace, whose conda_pkgs folder is CONDA_PKGS_DIRS.
#
#   -Phase before-install   lists the archives the cache restored (conda_restored.txt) and says
#                           how many and how large.
#   -Phase after-install    records the installed packages (conda_contents.txt, the sorted URLs of
#                           `conda list --explicit --md5`), reports how many archives were reused
#                           from the cache, downloaded in this job, or restored and left unused, and
#                           the versions of the packages the build depends on most.  With -Prune it
#                           deletes every archive not installed, so a saved cache holds this
#                           environment only; with -KeyPrefix it writes `key=<prefix>-<SHA-256 of
#                           conda_contents.txt>` to GITHUB_OUTPUT.

param(
    [Parameter(Mandatory)][ValidateSet('before-install', 'after-install')][string]$Phase,
    [string]$Environment = 'b2-windows',
    [string]$KeyPrefix = '',
    [switch]$Prune
)

$ErrorActionPreference = 'Stop'
$pkgsDir = 'conda_pkgs'

function Get-Archives {
    if (-not (Test-Path $pkgsDir)) { return @() }
    Get-ChildItem $pkgsDir -File | Where-Object { $_.Name -like '*.conda' -or $_.Name -like '*.tar.bz2' }
}

function Format-MB([long]$bytes) { '{0:N1} MB' -f ($bytes / 1MB) }

function Write-Group([string]$title, $archives) {
    $list = @($archives)
    $bytes = ($list | Measure-Object Length -Sum).Sum
    Write-Host ("{0}: {1} archives, {2}" -f $title, $list.Count, (Format-MB $bytes))
    foreach ($archive in ($list | Sort-Object Name | Select-Object -First 40)) {
        Write-Host ("    {0}  ({1})" -f $archive.Name, (Format-MB $archive.Length))
    }
    if ($list.Count -gt 40) { Write-Host "    ... and $($list.Count - 40) more" }
}

if ($Phase -eq 'before-install') {
    $restored = @(Get-Archives)
    $restored | ForEach-Object { "{0}`t{1}" -f $_.Name, $_.LastWriteTimeUtc.Ticks } | Set-Content conda_restored.txt
    $bytes = ($restored | Measure-Object Length -Sum).Sum
    Write-Host ("restored from the cache: {0} archives, {1}" -f $restored.Count, (Format-MB $bytes))
    return
}

conda list -n $Environment --explicit --md5 | Where-Object { $_ -match '^https://' } | Sort-Object | Set-Content conda_contents.txt
$installed = @(Get-Content conda_contents.txt | ForEach-Object { ($_ -split '#')[0].Split('/')[-1] })
$restoredTicks = @{}
if (Test-Path conda_restored.txt) {
    foreach ($line in Get-Content conda_restored.txt) { $name, $ticks = $line -split "`t"; $restoredTicks[$name] = [long]$ticks }
}
$present = @(Get-Archives)
function Test-AsRestored($archive) {
    $restoredTicks.ContainsKey($archive.Name) -and $restoredTicks[$archive.Name] -eq $archive.LastWriteTimeUtc.Ticks
}

Write-Host "installed packages: $($installed.Count)"
Write-Group 'reused from the cache' ($present | Where-Object { (Test-AsRestored $_) -and $installed -contains $_.Name })
Write-Group 'downloaded in this job (new, or rewritten after the restore)' ($present | Where-Object { -not (Test-AsRestored $_) })
$unused = @($present | Where-Object { $installed -notcontains $_.Name })
Write-Group 'present but not installed' $unused
$missing = @($installed | Where-Object { -not (Test-Path (Join-Path $pkgsDir $_)) })
Write-Host "installed with no archive in $pkgsDir`: $($missing.Count)"
foreach ($name in ($missing | Select-Object -First 20)) { Write-Host "    $name" }

Write-Host 'versions:'
conda list -n $Environment '^(python|clang|libboost-devel|eigen|eigenpy|numpy|scipy|libblas|cmake|ninja|ccache|msmpi|vs2022_win-64)$'

if ($Prune -and $unused.Count -gt 0) {
    $unused | Remove-Item
    Write-Host "pruned $($unused.Count) archives that are not installed"
}

if ($KeyPrefix) {
    $hash = (Get-FileHash conda_contents.txt -Algorithm SHA256).Hash
    "key=$KeyPrefix-$hash" >> $env:GITHUB_OUTPUT
    Write-Host "cache key: $KeyPrefix-$hash"
}
