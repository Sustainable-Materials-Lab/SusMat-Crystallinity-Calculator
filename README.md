# WAXS crystallinity calculator

## Description
This python script can be used to plot the fitted WAXS data output from TOPAS academic and determine the cellulose crystallinity using the Ruland method. 

## Installation
Requires python 3, matplotlib, scipy and numpy. In windows, installation of Anaconda python 3 distribution will cover all dependencies

## Usage
This requires a .txt file saved from TOPAS (academic) as input.

To plot the output from one fitting:
```python
python cell_cryst.py --xye --svg <filename without extension>
```

To plot the output from a folder full of files using powershell:

```powershell
Get-ChildItem -Filter *.txt | ForEach-Object -Process {python cell_cryst.py --xye --svg $_.BaseName}
```

## Support
Contact @u0092172

## Roadmap
No planned feature upgrades, reported bugs will be fixed.

## Contributing
Contributions are welcome, contact @u0092172

## Authors and acknowledgment
Authored by @u0092172

## License
Currently not published.
