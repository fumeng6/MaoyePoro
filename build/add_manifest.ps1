param(
    [Parameter(Mandatory = $true)]
    [string[]]$Executable,
    [string]$Manifest
)

$ErrorActionPreference = 'Stop'
if (-not $Manifest) { $Manifest = Join-Path $PSScriptRoot 'windows.manifest' }
if (-not (Test-Path -LiteralPath $Manifest -PathType Leaf)) {
    throw "Manifest not found: $Manifest"
}
$manifestBytes = [System.IO.File]::ReadAllBytes((Resolve-Path -LiteralPath $Manifest).ProviderPath)

if (-not ('MaoyePoroManifestResource' -as [type])) {
    Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;

public static class MaoyePoroManifestResource
{
    [DllImport("kernel32.dll", CharSet = CharSet.Unicode, ExactSpelling = true, SetLastError = true)]
    public static extern IntPtr BeginUpdateResourceW(string fileName,
        [MarshalAs(UnmanagedType.Bool)] bool deleteExistingResources);

    [DllImport("kernel32.dll", CharSet = CharSet.Unicode, ExactSpelling = true, SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool UpdateResourceW(IntPtr handle, IntPtr type,
        IntPtr name, ushort language, byte[] data, uint size);

    [DllImport("kernel32.dll", CharSet = CharSet.Unicode, ExactSpelling = true, SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool EndUpdateResourceW(IntPtr handle,
        [MarshalAs(UnmanagedType.Bool)] bool discard);
}
'@
}

foreach ($file in $Executable) {
    if (-not (Test-Path -LiteralPath $file -PathType Leaf)) {
        throw "Executable not found: $file"
    }
    $path = (Resolve-Path -LiteralPath $file).ProviderPath
    $handle = [MaoyePoroManifestResource]::BeginUpdateResourceW($path, $false)
    if ($handle -eq [IntPtr]::Zero) {
        throw [System.ComponentModel.Win32Exception]::new(
            [Runtime.InteropServices.Marshal]::GetLastWin32Error(),
            "Cannot open executable resources: $path")
    }
    try {
        # RT_MANIFEST=24, CREATEPROCESS_MANIFEST_RESOURCE_ID=1, neutral language.
        $updated = [MaoyePoroManifestResource]::UpdateResourceW(
            $handle, [IntPtr]24, [IntPtr]1, 0, $manifestBytes, [uint32]$manifestBytes.Length)
        if (-not $updated) {
            throw [System.ComponentModel.Win32Exception]::new(
                [Runtime.InteropServices.Marshal]::GetLastWin32Error(),
                "Cannot update executable manifest: $path")
        }
        $committed = [MaoyePoroManifestResource]::EndUpdateResourceW($handle, $false)
        $errorCode = [Runtime.InteropServices.Marshal]::GetLastWin32Error()
        $handle = [IntPtr]::Zero
        if (-not $committed) {
            throw [System.ComponentModel.Win32Exception]::new(
                $errorCode, "Cannot save executable manifest: $path")
        }
    }
    finally {
        if ($handle -ne [IntPtr]::Zero) {
            $discarded = [MaoyePoroManifestResource]::EndUpdateResourceW($handle, $true)
            if (-not $discarded) {
                Write-Warning "Cannot discard pending executable resource update: $path"
            }
        }
    }
    Write-Output "Embedded UTF-8 process manifest: $path"
}
