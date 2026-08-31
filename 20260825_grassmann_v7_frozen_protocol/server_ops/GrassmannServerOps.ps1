[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [ValidateSet('Plan', 'UploadV7', 'Bootstrap', 'Audit', 'T00', 'GpuTest', 'T03')]
    [string]$Action = 'Plan',

    [string]$RemoteHost = '166.111.94.174',
    [string]$RemoteUser = 'tanyuxiao',
    [string]$RemoteRoot = '/data1/home/tanyuxiao/Grassmann_model',
    [string]$ReleaseId = (Get-Date -Format 'yyyyMMddTHHmmss'),
    [string]$BundlePath = '',
    [string]$RemotePython = '',
    [string]$T03OutputId = ''
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

if (-not $BundlePath) {
    $BundlePath = Join-Path $PSScriptRoot '..\..\v7_p0_server_bundle_20260825.tar.gz'
}

if ($RemoteRoot -ne '/data1/home/tanyuxiao/Grassmann_model') {
    throw "RemoteRoot must remain /data1/home/tanyuxiao/Grassmann_model"
}
if ($RemoteHost -notmatch '^[A-Za-z0-9.-]+$') { throw 'Unsafe RemoteHost.' }
if ($RemoteUser -notmatch '^[A-Za-z0-9._-]+$') { throw 'Unsafe RemoteUser.' }
if ($ReleaseId -notmatch '^[A-Za-z0-9._-]+$') { throw 'Unsafe ReleaseId.' }
if ($T03OutputId -and $T03OutputId -notmatch '^[A-Za-z0-9._-]+$') { throw 'Unsafe T03OutputId.' }

$SshExe = Join-Path $env:WINDIR 'System32\OpenSSH\ssh.exe'
$ScpExe = Join-Path $env:WINDIR 'System32\OpenSSH\scp.exe'
$RemoteTarget = "${RemoteUser}@${RemoteHost}"
$ReleaseDir = "${RemoteRoot}/v7/code/releases/${ReleaseId}"
$IncomingFile = "${RemoteRoot}/incoming/v7_${ReleaseId}.tar.gz"
if (-not $RemotePython) {
    $RemotePython = "${ReleaseDir}/.venv/bin/python"
}
if (-not $T03OutputId) {
    $T03OutputId = "t03_profile_${ReleaseId}"
}

function Invoke-RemoteCommand {
    param([Parameter(Mandatory)][string]$Command)
    Write-Host "REMOTE> $Command"
    if ($PSCmdlet.ShouldProcess($RemoteTarget, $Command)) {
        & $SshExe -o StrictHostKeyChecking=yes $RemoteTarget $Command
        if ($LASTEXITCODE -ne 0) {
            throw "SSH command failed with exit code $LASTEXITCODE"
        }
    }
}

function Invoke-UploadV7 {
    $ResolvedBundle = (Resolve-Path -LiteralPath $BundlePath).Path
    $BundleHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $ResolvedBundle).Hash.ToLowerInvariant()
    $Prepare = "umask 077; mkdir -p '${RemoteRoot}/incoming' '${RemoteRoot}/v7/code/releases'; " +
        "if [ -e '${ReleaseDir}' ]; then echo 'Release already exists; refusing overwrite' >&2; exit 41; fi"
    Invoke-RemoteCommand $Prepare

    Write-Host "UPLOAD> $ResolvedBundle -> ${RemoteTarget}:$IncomingFile"
    if ($PSCmdlet.ShouldProcess("${RemoteTarget}:$IncomingFile", 'Upload immutable v7 bundle')) {
        & $ScpExe -o StrictHostKeyChecking=yes -- $ResolvedBundle ($RemoteTarget + ':' + $IncomingFile)
        if ($LASTEXITCODE -ne 0) {
            throw "SCP failed with exit code $LASTEXITCODE"
        }
    }

    $Extract = "set -e; actual=`$(sha256sum '${IncomingFile}' | cut -d ' ' -f1); " +
        "test `"`$actual`" = '${BundleHash}' || { echo 'Bundle hash mismatch' >&2; exit 42; }; " +
        "mkdir '${ReleaseDir}'; tar -xzf '${IncomingFile}' --strip-components=1 -C '${ReleaseDir}'; " +
        "cd '${ReleaseDir}'; python3 validate_freeze.py; python3 p0/write_implementation_manifest.py; " +
        "sha256sum -c p0/P0_IMPLEMENTATION.sha256"
    Invoke-RemoteCommand $Extract
}

function Invoke-Bootstrap {
    Invoke-RemoteCommand "bash '${ReleaseDir}/server_ops/bootstrap_v6_v7.sh' '${RemoteRoot}'"
}

function Invoke-Audit {
    $AuditId = "audit_${ReleaseId}"
    $AuditPath = "${RemoteRoot}/v7/resources/audits/${AuditId}.json"
    Invoke-RemoteCommand (
        "umask 077; mkdir -p '${RemoteRoot}/v7/resources/audits'; " +
        "python3 '${ReleaseDir}/p0/audit_server.py' --output '${AuditPath}'"
    )
}

function Invoke-T00 {
    Invoke-RemoteCommand "cd '${ReleaseDir}' && bash p0/run_t00_nonintrusive.sh"
}

function Invoke-GpuTest {
    Invoke-RemoteCommand (
        "bash '${ReleaseDir}/server_ops/run_gpu_test_nonintrusive.sh' " +
        "--root '${RemoteRoot}' --python '${RemotePython}'"
    )
}

function Invoke-T03 {
    $Output = "${RemoteRoot}/v7/results/${T03OutputId}"
    Invoke-RemoteCommand "cd '${ReleaseDir}' && bash p0/run_t03_nonintrusive.sh '${Output}'"
}

switch ($Action) {
    'Plan' {
        Write-Host 'No server connection was made. Prepared values:'
        [pscustomobject]@{
            RemoteTarget = $RemoteTarget
            RemoteRoot = $RemoteRoot
            ReleaseId = $ReleaseId
            ReleaseDir = $ReleaseDir
            BundlePath = $BundlePath
            RemotePython = $RemotePython
            T03Output = "${RemoteRoot}/v7/results/${T03OutputId}"
        } | Format-List
        Write-Host 'Recommended order: UploadV7 -> Bootstrap -> Audit -> T00 -> GpuTest -> T03'
        Write-Host 'Audit/T00/GpuTest/T03 all retain records; GPU actions refuse to start when no eligible idle GPU exists.'
    }
    'UploadV7' { Invoke-UploadV7 }
    'Bootstrap' { Invoke-Bootstrap }
    'Audit' { Invoke-Audit }
    'T00' { Invoke-T00 }
    'GpuTest' { Invoke-GpuTest }
    'T03' { Invoke-T03 }
}
