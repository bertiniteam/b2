# Finds the link tools of the Visual Studio whose STL clang-cl compiles against, for the Windows CI
# builds (ADR-0078): link.exe and lib.exe from its default MSVC toolset, rc.exe and mt.exe from the
# newest Windows 10/11 SDK.  Returns a hashtable of CMake variable -> full path.
#
# clang-cl finds Visual Studio by itself, but the conda environment's activation does not run
# Visual Studio's setup script in PowerShell, so none of these is on PATH (and Git's link.exe,
# which is not a linker, may be).  Left to itself, CMake takes the runner's own LLVM tools from
# C:\Program Files\LLVM, which nothing pins.

$ErrorActionPreference = 'Stop'

$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
$vsPath = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vsPath) { throw "vswhere found no Visual Studio with the MSVC x64 tools" }

$toolsetVersion = (Get-Content (Join-Path $vsPath 'VC\Auxiliary\Build\Microsoft.VCToolsVersion.default.txt')).Trim()
$msvcBin = Join-Path $vsPath "VC\Tools\MSVC\$toolsetVersion\bin\Hostx64\x64"

$sdkRoot = Join-Path ${env:ProgramFiles(x86)} 'Windows Kits\10\bin'
$sdkBin = Get-ChildItem $sdkRoot -Directory |
    Where-Object { $_.Name -match '^10\.' -and (Test-Path (Join-Path $_.FullName 'x64\rc.exe')) } |
    Sort-Object { [version]$_.Name } |
    Select-Object -Last 1 |
    ForEach-Object { Join-Path $_.FullName 'x64' }
if (-not $sdkBin) { throw "no Windows SDK with x64\rc.exe under $sdkRoot" }

$tools = [ordered]@{
    CMAKE_LINKER      = Join-Path $msvcBin 'link.exe'
    CMAKE_AR          = Join-Path $msvcBin 'lib.exe'
    CMAKE_RC_COMPILER = Join-Path $sdkBin 'rc.exe'
    CMAKE_MT          = Join-Path $sdkBin 'mt.exe'
}
foreach ($entry in $tools.GetEnumerator()) {
    if (-not (Test-Path $entry.Value)) { throw "$($entry.Key): not found at $($entry.Value)" }
    Write-Host "$($entry.Key) = $($entry.Value)"
}
$tools
