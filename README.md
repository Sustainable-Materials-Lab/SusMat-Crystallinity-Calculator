# WAXS crystallinity calculator

## Description
This python script can be used to plot the fitted WAXS data output from TOPAS academic and determine the cellulose crystallinity using the Ruland method. 

## Installation
Download the repo and install with pip.

## Usage
This requires a .txt file saved from TOPAS (academic) as input.

### Local usage
To plot the output from one fitting:
```powershell
sm-cryst --svg filename.txt
```

To plot the output from a folder full of files using powershell:

```powershell
Get-ChildItem -Filter *.txt | ForEach-Object -Process {sm-cryst --svg $_.BaseName}
```

### Remote usage (iRODS/ManGO)
Calculate the crystallinity for a collection (folder) containing many .txt files in iRODS. Ensure that the txt files and dat files have the same name in order to copy metadata onto the resulting image files. Enter the full path to the collection in iRODS:

```powershell
irods-cryst --svg /set/home/SusMat/$project_ID/$initials/SWAXS/$target_collection
```

## Support
Contact @u0092172

## Roadmap
No planned feature upgrades, reported bugs will be fixed.

## Contributing
Contributions are welcome, contact @u0092172

## Authors and acknowledgment
Authored by @u0092172

**AI Assistance**: Code refactoring and architecture improvements contributed by GitHub Copilot (October 2025)

## License
Currently not published.
