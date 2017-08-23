################################# [ quickIndex.py ] ################################
#	Author: Eduardo (sacridini) Lacerda
#	e-mail: eduardolacerdageo@gmail.com
#	Version: 0.1.7.2
#
#	Calculate multispectral indices such as NDVI, SAVI, TVI, etc.
#	Dependencies: numpy, gdal and rasterio
# 
#################################################################################

import numpy as np
import rasterio

class QuickIndex(object):

	# SR - Simple Ratio Vegetation Index
	# nir/red
	def sr(self, red, nir, kwargs):
		sr = (self.nir.astype(float) / self.red.astype(float))
		with rasterio.open('sr_specidx.tif', 'w', **kwargs) as dst_sr:
			dst_sr.write_band(1, sr.astype(rasterio.float32))

	# NDVI - Normalised Difference Vegetation Index
	# (nir - red)/(nir + red)
	def ndvi(self, red, nir, kwargs):
		ndvi = (self.nir.astype(float) - self.red.astype(float)) / (self.nir + self.red)
		with rasterio.open('ndvi_specidx.tif', 'w', **kwargs) as dst_ndvi:
			dst_ndvi.write_band(1, ndvi.astype(rasterio.float32))

	# TVI - Transformed Vegetation Index
	# sqrt((nir - red)/(nir + red) + 0.5)
	def tvi(self, red, nir, kwargs):
		tvi = np.sqrt((self.nir.astype(float) - self.red.astype(float)) / (self.nir + self.red) + 0.5)
		with rasterio.open('tvi_specidx.tif', 'w', **kwargs) as dst_tvi:
			dst_tvi.write_band(1, tvi.astype(rasterio.float32))
				
	# TTVI - Thiam's Transformed Vegetation Index
	# sqrt(abs((nir - red)/(nir + red) + 0.5))
	def ttvi(self, red, nir, kwargs):
		ttvi = np.sqrt(np.absolute((self.nir.astype(float) - self.red.astype(float)) / (self.nir + self.red) + 0.5))
		with rasterio.open('ttvi_specidx.tif', 'w', **kwargs) as dst_ttvi:
			dst_ttvi.write_band(1, ttvi.astype(rasterio.float32))
	
	# NRVI - Normalised Ratio Vegetation Index
	# (red/nir - 1)/(red/nir + 1)
	def nrvi(self, red, nir, kwargs):
		nrvi = (self.red.astype(float) / self.nir.astype(float) - 1) / ((self.red / self.nir) + 1) 
		with rasterio.open('nrvi_specidx.tif', 'w', **kwargs) as dst_nrvi:
			dst_nrvi.write_band(1, nrvi.astype(rasterio.float32))

	# RVI - Ratio Vegetation Index
	# red / nir
	def rvi(self, red, nir, kwargs):	
		rvi = (self.red.astype(float) / self.nir.astype(float))
		with rasterio.open('rvi_specidx.tif', 'w', **kwargs) as dst_rvi:
			dst_rvi.write_band(1, rvi.astype(rasterio.float32))		

	# SAVI - Soil Adjusted Vegetation Index
	# (nir - red) * (1 + L)/(nir + red + L)
	def savi(self, red, nir, kwargs):
		savi = (self.nir.astype(float) - self.red.astype(float) * (1 + 0.5)) / (self.nir + self.red + 0.5)
		with rasterio.open('savi_specidx.tif', 'w', **kwargs) as dst_savi:
			dst_savi.write_band(1, savi.astype(rasterio.float32))

	# MSAVI - Modified Soil Adjusted Vegetation Index
	# nir + 0.5 - (0.5 * sqrt((2 * nir + 1)^2 - 8 * (nir - (2 * red))))
	def msavi(self, red, nir, kwargs):
		msavi = self.nir.astype(float) + 0.5 - (0.5 * np.sqrt((2 * self.nir.astype(float) + 1)**2 - 8 * (self.nir.astype(float) - (2 * self.red.astype(float)))))
		with rasterio.open('msavi_specidx.tif', 'w', **kwargs) as dst_msavi:
			dst_msavi.write_band(1, msavi.astype(rasterio.float32))		

	# Modified Soil Adjusted Vegetation Index 2
	# (2 * (nir + 1) - sqrt((2 * nir + 1)^2 - 8 * (nir - red)))/2
	def msavi2(self, red, nir, kwargs):
		msavi2 = (2 * (self.nir.astype(float) + 1) - np.sqrt((2 * self.nir.astype(float) + 1)**2 - 8 * (self.nir.astype(float) - self.red.astype(float))))/2
		with rasterio.open('msavi2_specidx.tif', 'w', **kwargs) as dst_msavi2:
			dst_msavi2.write_band(1, msavi2.astype(rasterio.float32))	
	# Global Environmental Monitoring Index
	# (((nir^2 - red^2) * 2 + (nir * 1.5) + (red * 0.5))/(nir + red + 0.5)) * (1 - ((((nir^2 - red^2) * 2 + (nir * 1.5) + (red * 0.5))/(nir + red + 0.5)) * 0.25)) - ((red - 0.125)/(1 - red))
	def gemi(self, red, nir, kwargs):
		gemi = (((self.nir.astype(float)**2 - self.red.astype(float)**2) * 2 + (self.nir * 1.5) + (self.red * 0.5))/(self.nir + self.red + 0.5)) * (1 - ((((self.nir**2 - self.red**2) * 2 + (self.nir * 1.5) + (self.red * 0.5))/(self.nir + self.red + 0.5)) * 0.25)) - ((self.red - 0.125)/(1 - self.red))
		with rasterio.open('gemi_specidx.tif', 'w', **kwargs) as dst_gemi:
			dst_gemi.write_band(1, gemi.astype(rasterio.float32))

	# Corrected Transformed Vegetation Index
	# (NDVI + 0.5)/sqrt(abs(NDVI + 0.5))
	def ctvi(self, red, nir, kwargs):
		ndvi = (self.nir.astype(float) - self.red.astype(float)) / (self.nir + self.red)
		ctvi = (ndvi + 0.5) / np.sqrt(np.absolute(ndvi + 0.5))
		with rasterio.open('ctvi_specidx.tif', 'w', **kwargs) as dst_ctvi:
			dst_ctvi.write_band(1, ctvi.astype(rasterio.float32))		

	# Difference Vegetation Index
	# s * nir - red
	def dvi(self, red, nir, kwargs):
		s = self.nir.astype(float) / self.red.astype(float)
		dvi = s * self.nir.astype(float) - self.red.astype(float)
		with rasterio.open('dvi_specidx.tif', 'w', **kwargs) as dst_dvi:
			dst_dvi.write_band(1, dvi.astype(rasterio.float32))		
	 
	# Weighted Difference Vegetation Index
	# nir - s * red
	def wdvi(self, red, nir, kwargs):
		s = self.nir.astype(float) / self.red.astype(float)
		wdvi = self.nir.astype(float) - s * self.red.astype(float)
		with rasterio.open('wdvi_specidx.tif', 'w', **kwargs) as dst_wdvi:
			dst_wdvi.write_band(1, wdvi.astype(rasterio.float32))

	# Green Normalised Difference Vegetation Index 
	# (nir - green)/(nir + green)
	def gndvi(self, green, nir, kwargs):
		gndvi = (self.nir.astype(float) - self.green.astype(float)) / (self.nir + self.green)
		with rasterio.open('gndvi_specidx.tif', 'w', **kwargs) as dst_gndvi:
			dst_gndvi.write_band(1, gndvi.astype(rasterio.float32))

	# Normalised Difference Water Index
	# (green - nir)/(green + nir)
	def ndwi(self, green, nir, kwargs):
		ndwi = (self.green.astype(float) - self.nir.astype(float)) / (self.green + self.nir)
		with rasterio.open('ndwi_specidx.tif', 'w', **kwargs) as dst_ndwi:
			dst_ndwi.write_band(1, ndwi.astype(rasterio.float32))

	# Normalised Difference Water Index 2
	# (nir - swir2)/(nir + swir2)
	def ndwi2(self, nir, swir, kwargs):
		ndwi2 = (self.nir.astype(float) - self.swir.astype(float)) / (self.nir + self.swir)
		with rasterio.open('ndwi2_specidx.tif', 'w', **kwargs) as dst_ndwi2:
			dst_ndwi2.write_band(1, ndwi2.astype(rasterio.float32))

	# Modified Normalised Difference Water Index 
	# (green - swir2)/(green + swir2)
	def mndwi(self, green, swir, kwargs):
		mndwi = (self.green.astype(float) - self.swir.astype(float)) / (self.green + self.swir)
		with rasterio.open('mndwi_specidx.tif', 'w', **kwargs) as dst_mndwi:
			dst_mndwi.write_band(1, mndwi.astype(rasterio.float32))

	# Corrected Normalised Difference Vegetation Index
	# (nir - red)/(nir + red) * (1 - ((swir2 - swir2ccc)/(swir2coc - swir2ccc)))
	def ndvic(self, red, nir, swir, kwargs):
		pass

	# Specific Leaf Area Vegetation Index
	# nir/(red + swir2) 
	def slavi(self, red, nir, swir, kwargs):
		slavi = self.nir.astype(float) / (self.red.astype(float) + self.swir.astype(float))
		with rasterio.open('slavi_specidx.tif', 'w', **kwargs) as dst_slavi:
			dst_slavi.write_band(1, slavi.astype(rasterio.float32))


	# Create all indices that use the RED and the NIR bands
	def genAllRedNir(self, red, nir, kwargs):
		self.savi(self.red, self.nir, self.kwargs)
		self.msavi(self.red, self.nir, self.kwargs)
		self.msavi2(self.red, self.nir, self.kwargs)
		self.gemi(self.red, self.nir, self.kwargs)
		self.ctvi(self.red, self.nir, self.kwargs)
		self.sr(self.red, self.nir, self.kwargs)
		self.dvi(self.red, self.nir, self.kwargs)
		self.rvi(self.red, self.nir, self.kwargs)
		self.tvi(self.red, self.nir, self.kwargs)
		self.ttvi(self.red, self.nir, self.kwargs)
		self.nrvi(self.red, self.nir, self.kwargs)
		self.wdvi(self.red, self.nir, self.kwargs)		

	# Create all indices that use the GREEN and the NIR bands
	def genAllGreenNir(self, green, nir, kwargs):
		self.gndvi(self.green, self.nir, self.kwargs)
		self.ndwi(self.green, self.nir, self.kwargs)
	
	# Create all indices that use the GREEN and the SWIR bands
	def genAllGreenSwir(self, green, swir, kwargs):
		self.mndwi(self.green, self.swir, self.kwargs)

	# Create all indices that use the RED, NIR and the SWIR bands
	def getAllRedNirSwir(self, red, nir, swir, kwargs):
		self.slavi(self.red, self.nir, self.swir, self.kwargs)

	# Create all indices that use the NIR and the SWIR bands
	def genAllNirSwir(self, nir, swir, kwargs):
		self.ndwi2(self.nir, self.swir, self.kwargs)

	# Constructor
	def __init__(self, indices = None, green = None, red = None, nir = None, swir = None):
		
		self.green_filepath = green
		self.red_filepath = red
		self.nir_filepath = nir
		self.swir_filepath = swir

		if green is not None and red is not None and nir is not None and swir is not None:
			with rasterio.open(green) as src_green:
				self.green = src_green.read(1)
			with rasterio.open(red) as src_red:
				self.red = src_red.read(1)
			with rasterio.open(nir) as src_nir:
				self.nir = src_nir.read(1)
			with rasterio.open(swir) as src_swir:
				self.swir = src_swir.read(1)
			self.kwargs = src_red.meta
			self.kwargs.update(
			dtype=rasterio.float32,
				count = 1)
			if indices == None:
				self.getAllRedNir(self.red, self.nir, self.kwargs)
				self.getAllRedNirSwir(self.red, self.nir, self.swir, self.kwargs)
				self.genAllGreenNir(self.green, self.nir, self.kwargs)
				self.genAllGreenSwir(self.green, self.swir, self.kwargs)
			else:
				print "Error"

		elif red is not None and nir is not None and swir is not None:
			with rasterio.open(red) as src_red:
				self.red = src_red.read(1)
			with rasterio.open(nir) as src_nir:
				self.nir = src_nir.read(1)
			with rasterio.open(swir) as src_swir:
				self.swir = src_swir.read(1)
			self.kwargs = src_red.meta
			self.kwargs.update(
				dtype=rasterio.float32,
				count = 1)

			# Check if any specific index was requested
			# otherwise, it will generate all possible indexes
			if indices == None:
				self.getAllRedNir(self.red, self.nir, self.kwargs)
				self.getAllRedNirSwir(self.red, self.nir, self.swir, self.kwargs)
			else:
				for idx in indices:
					if idx == 'slavi':
						self.slavi(self.red, self.nir, self.swir, self.kwargs)
					else:
						print "Error: Index [" + idx + "] doesnt exist"

		elif green is not None and nir is not None:
			with rasterio.open(green) as src_green:
				self.green = src_green.read(1)
			with rasterio.open(nir) as src_nir:
				self.nir = src_nir.read(1)
			self.kwargs = src_nir.meta
			self.kwargs.update(
				dtype=rasterio.float32,
				count = 1)

			# Check if any specific index was requested
			# otherwise, it will generate all possible indexes
			if indices == None:
				self.genAllGreenNir(self.green, self.nir, self.kwargs)
			else:
				for idx in indices:
					if idx == 'gndvi':
						self.gndvi(self.green, self.nir, self.kwargs)
					elif idx == 'ndwi':
						self.ndwi(self.green, self.nir, self.kwargs)
					else:
						print "Error: Index [" + idx + "] doesnt exist"


		elif green is not None and swir is not None:
			with rasterio.open(green) as src_green:
				self.green = src_green.read(1)
			with rasterio.open(swir) as src_swir:
				self.swir = src_swir.read(1)
			self.kwargs = src_green.meta
			self.kwargs.update(
				dtype=rasterio.float32,
				count = 1)

			# Check if any specific index was requested
			# otherwise, it will generate all possible indexes			
			if indices == None:
				self.genAllGreenSwir(self.green, self.swir, self.kwargs)
			else:
				for idx in indices:
					if idx == 'mndwi':
						self.mndwi(self.green, self.swir, self.kwargs)
					else:
						print "Error: Index [" + idx + "] doesnt exist"


		elif red is not None and nir is not None:
			with rasterio.open(red) as src_red:
				self.red = src_red.read(1)
			with rasterio.open(nir) as src_nir:
					self.nir = src_nir.read(1)
			self.kwargs = src_red.meta
			self.kwargs.update(
				dtype=rasterio.float32,
				count = 1)

			# Check if any specific index was requested
			# otherwise, it will generate all possible indexes
			if indices == None:
				self.genAllRedNir(self.red, self.nir, self.kwargs)
			else:
				for idx in indices:
					if idx == 'ndvi':
						self.ndvi(self.red, self.nir, self.kwargs)
					elif idx == 'sr':
						self.sr(self.red, self.nir, self.kwargs)
					elif idx == 'tvi':
						self.tvi(self.red, self.nir, self.kwargs)
					elif idx == 'ttvi':
						self.ttvi(self.red, self.nir, self.kwargs)
					elif idx == 'nrvi':
						self.nrvi(self.red, self.nir, self.kwargs)
					elif idx == 'rvi':
						self.rvi(self.red, self.nir, self.kwargs)
					elif idx == 'savi':
						self.savi(self.red, self.nir, self.kwargs)
					elif idx == 'msavi':
						self.msavi(self.red, self.nir, self.kwargs)
					elif idx == 'msavi2':
						self.msavi2(self.red, self.nir, self.kwargs)
					elif idx == 'gemi':
						self.msavi2(self.red, self.nir, self.kwargs)
					elif idx == 'dvi':
						self.msavi2(self.red, self.nir, self.kwargs)
					elif idx == 'wdvi':
						self.msavi2(self.red, self.nir, self.kwargs)	
					elif idx == 'ctvi':
						self.msavi2(self.red, self.nir, self.kwargs)																		
					else:
						print "Error: Index [" + idx + "] cant be created"
						print_one = "Try one of these: sr, ndvi, rvi, tvi, ttvi, nrvi,"
						print_two = " savi, msavi, msavi2, gemi, dvi, wdvi or ctvi"
						print print_one + print_two

		elif nir is not None and swir is not None:
			with rasterio.open(nir) as src_nir:
				nir = src_nir.read(1)
			with rasterio.open(swir) as src_swir:
				swir = src_swir.read(1)
				self.kwargs = src_red.meta
				self.kwargs.update(
					dtype=rasterio.float32,
					count = 1)
			if indices == None:
				self.genAllNirSwir(self.nir, self.swir, self.kwargs)
			else:
				for idx in indices:
					if idx == 'ndwi2':
						self.ndwi2(self.nir, self.swir, self.kwargs)
					else:
						print "Error: Index [" + idx + "] cant be created"
		else:
			print "------------- [ Welcome to QuickIndex ] -------------"
			print

# UNIX
# green_file = '/home/eduardo/Documents/images/LT05_L1TP_217076_20110813_20161007_01_T1_B2.TIF'
# red_file = '/home/eduardo/Documents/images/LT05_L1TP_217076_20110813_20161007_01_T1_B3.TIF'
# nir_file = '/home/eduardo/Documents/images/LT05_L1TP_217076_20110813_20161007_01_T1_B4.TIF'

# WIN
# green_file = 'C:\\Users\\eduardo\\Documents\\LT05_L1TP_217076_20110813_20161007_01_T1_B2.TIF'
# red_file = 'C:\\Users\\eduardo\\Documents\\LT05_L1TP_217076_20110813_20161007_01_T1_B3.TIF'
# nir_file = 'C:\\Users\\eduardo\\Documents\\LT05_L1TP_217076_20110813_20161007_01_T1_B4.TIF'

# idx = ["gndvi"]
# QuickIndex(idx, green = green_file, nir = nir_file)