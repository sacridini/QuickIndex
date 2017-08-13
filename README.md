# QuickIndex

![28606761-031da9b8-71af-11e7-8e4a-3a716e8a9886](https://user-images.githubusercontent.com/7756611/29248089-40c76e2c-7fe5-11e7-8c3f-73400c6f2c67.jpg)

## Dependences are:
- Numpy
- Gdal
- Rasterio

## Usage
QuickIndex(indices = None, green = None, red = None, nir = None, swir = None)

. indices: List of indexes like: ["ndvi", "tvi", "savi"]
  possible choices: NDVI; SR; TVI; TTVI; SAVI; MSAVI; MSAVI2, RVI, DVI, CTVI, GEMI, NRVI, WDVI
  
. green: Character or integer. Green band.

. red: Character or integer. Red band.

. nir: Character or integer. Near-infrared band

. swir: Character or integer. Short-wave-infrared band (1400-1800nm).

### Example
```
idx = ["ndvi"]
QuickIndex(idx, red = red_file, nir = nir_file)
```
