$ErrorActionPreference = 'Stop'
$bundleRoot = [IO.Path]::GetFullPath($PSScriptRoot).TrimEnd('\')

foreach ($program in 'myporo.exe', 'mycfs.exe') {
    $programPath = Join-Path $bundleRoot "bin\$program"
    if (-not (Test-Path -LiteralPath $programPath -PathType Leaf)) {
        throw "未找到 $programPath。请先完整解压程序包，再运行本脚本。"
    }
}

$userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
$userEntries = @($userPath -split ';' | ForEach-Object {
    [Environment]::ExpandEnvironmentVariables($_.Trim().Trim('"')).TrimEnd('\')
})
if ($userEntries -notcontains $bundleRoot) {
    $newUserPath = if ([string]::IsNullOrWhiteSpace($userPath)) {
        $bundleRoot
    } else {
        $userPath.TrimEnd(';') + ';' + $bundleRoot
    }
    [Environment]::SetEnvironmentVariable('Path', $newUserPath, 'User')
}

$processEntries = @($env:Path -split ';' | ForEach-Object {
    [Environment]::ExpandEnvironmentVariables($_.Trim().Trim('"')).TrimEnd('\')
})
if ($processEntries -notcontains $bundleRoot) {
    $env:Path = $env:Path.TrimEnd(';') + ';' + $bundleRoot
}

Write-Host "已将程序目录加入当前用户 PATH：$bundleRoot"
Write-Host '请重新打开终端，运行 myporo --help。库仑应力命令为 mycfs。'
Write-Host '请保留整个程序文件夹；移动文件夹后需从新位置重新运行本脚本。'
