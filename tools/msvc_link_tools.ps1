# Enters the developer environment of the Visual Studio whose STL clang-cl compiles against, for the
# Windows CI builds (ADR-0078), and returns its link tools: a hashtable of CMake variable -> full path
# for link.exe and lib.exe (its MSVC toolset) and rc.exe and mt.exe (its Windows SDK).
#
# clang-cl finds Visual Studio by itself, and lld-link finds the Windows libraries by itself, but
# link.exe needs the developer environment (LIB above all), and the conda activation does not run
# Visual Studio's setup in PowerShell.  Left to itself, CMake takes the runner's own LLVM tools from
# C:\Program Files\LLVM, which nothing pins.  Run it with `&` in the step that configures: the
# environment it sets lasts for that step's process.  Full paths, because Git's link.exe, which is not
# a linker, can come first on PATH.

$ErrorActionPreference = 'Stop'

$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
$vsPath = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vsPath) { throw "vswhere found no Visual Studio with the MSVC x64 tools" }
Write-Host "Visual Studio: $vsPath"

# The developer shell also puts Visual Studio's own CMake, Ninja and tools first on PATH; the build
# must keep conda's cmake, ninja, ccache and clang-cl, so PATH goes back as it was.  What stays is
# LIB, LIBPATH and INCLUDE (the libraries and headers of one MSVC toolset and one SDK) and the
# variables naming them; the tools below are passed to CMake by full path.
$pathBefore = $env:PATH
& (Join-Path $vsPath 'Common7\Tools\Launch-VsDevShell.ps1') -Arch amd64 -HostArch amd64 -SkipAutomaticLocation | Out-Host
$env:PATH = $pathBefore
foreach ($name in 'VCToolsInstallDir', 'WindowsSdkVerBinPath', 'LIB') {
    if (-not (Get-Item "Env:$name" -ErrorAction SilentlyContinue)) { throw "the developer shell did not set $name" }
}
Write-Host "LIB = $env:LIB"
Write-Host "INCLUDE = $env:INCLUDE"

$msvcBin = Join-Path $env:VCToolsInstallDir 'bin\Hostx64\x64'
$sdkBin = Join-Path $env:WindowsSdkVerBinPath 'x64'
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
