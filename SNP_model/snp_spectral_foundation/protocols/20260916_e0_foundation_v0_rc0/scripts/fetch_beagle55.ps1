$ErrorActionPreference = 'Stop'

$destination = 'F:\Yuxiao Tan\Document\Grassmann_model_external\beagle_5.5_27Feb25.75f'
New-Item -ItemType Directory -Path $destination -Force | Out-Null

$assets = @(
    @{
        Url = 'https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar'
        Name = 'beagle.27Feb25.75f.jar'
    },
    @{
        Url = 'https://faculty.washington.edu/browning/beagle/beagle_5.5_17Dec24.pdf'
        Name = 'beagle_5.5_17Dec24.pdf'
    },
    @{
        Url = 'https://faculty.washington.edu/browning/beagle/run.beagle.27Feb25.75f.example'
        Name = 'run.beagle.27Feb25.75f.example'
    }
)

foreach ($asset in $assets) {
    $target = Join-Path $destination $asset.Name
    Invoke-WebRequest -Uri $asset.Url -OutFile $target
}
