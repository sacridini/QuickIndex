# QuickIndex

Dependences are:
Numpy
Gdal
Rasterio

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
