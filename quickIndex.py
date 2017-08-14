################################# [ quickIndex.py ] ################################
#	Author: Eduardo (sacridini) Lacerda
#	e-mail: eduardolacerdageo@gmail.com
#	Version: 0.1.2
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
		msavi = self.nir.astype(float) + 0.5 - (0.5 * np.sqrt((2 * self.nir.astype(float) + 1)**2 - 8 * (nir - (2 * self.red.astype(float)))))
		with rasterio.open('msavi_specidx.tif', 'w', **kwargs) as dst_msavi:
			dst_msavi.write_band(1, msavi.astype(rasterio.float32))		

	# Modified Soil Adjusted Vegetation Index 2
	# (2 * (nir + 1) - sqrt((2 * nir + 1)^2 - 8 * (nir - red)))/2
	def msavi2(self):
		pass # TODO	

	# Global Environmental Monitoring Index
	# (((nir^2 - red^2) * 2 + (nir * 1.5) + (red * 0.5))/(nir + red + 0.5)) * (1 - ((((nir^2 - red^2) * 2 + (nir * 1.5) + (red * 0.5))/(nir + red + 0.5)) * 0.25)) - ((red - 0.125)/(1 - red))
	def gemi(self):
		pass # TODO

	# Corrected Transformed Vegetation Index
	# (NDVI + 0.5)/sqrt(abs(NDVI + 0.5))
	def ctvi(self):
		pass # TODO

	# Difference Vegetation Index
	# s * nir - red
	def dvi(self):
		pass # TODO

	# Weighted Difference Vegetation Index
	# nir - s * red
	def wdvi(self):
		pass # TODO

	# Create all indices that use the RED and the NIR bands
	def genAllRedNir(self, red, nir, kwargs):
		self.savi(self.red, self.nir, self.kwargs)
		self.msavi(self.red, self.nir, self.kwargs)
		# self.msavi2(self.red, self.nir, self.kwargs)
		# self.gemi(self.red, self.nir, self.kwargs)
		# self.ctvi(self.red, self.nir, self.kwargs)
		self.sr(self.red, self.nir, self.kwargs)
		# self.dvi(self.red, self.nir, self.kwargs)
		self.rvi(self.red, self.nir, self.kwargs)
		self.tvi(self.red, self.nir, self.kwargs)
		self.ttvi(self.red, self.nir, self.kwargs)
		self.nrvi(self.red, self.nir, self.kwargs)
		# self.wdvi(self.red, self.nir, self.kwargs)		


	def __init__(self, indices = None, green = None, red = None, nir = None, swir = None):
		
		self.green_filepath = green
		self.red_filepath = red
		self.nir_filepath = nir
		self.swir_filepath = swir

		if green is not None and red is not None and nir is not None and swir is not None:
		 	with rasterio.open(green_filepath) as src_green:
				green = src_green.read(1)
			with rasterio.open(red_filepath) as src_red:
				red = src_red.read(1)
			with rasterio.open(nir_filepath) as src_nir:
				nir = src_nir.read(1)
			with rasterio.open(swir_filepath) as src_swir:
				swir = src_swir.read(1)
				# TODO 

		elif red is not None and nir is not None and swir is not None:
			with rasterio.open(red_filepath) as src_red:
				red = src_red.read(1)
			with rasterio.open(nir_filepath) as src_nir:
				nir = src_nir.read(1)
			with rasterio.open(swir_filepath) as src_swir:
				swir = src_swir.read(1)
				# TODO 

		elif green is not None and nir is not None:
			with rasterio.open(green_filepath) as src_green:
				green = src_green.read(1)
			with rasterio.open(nir_filepath) as src_nir:
				nir = src_nir.read(1)
				# TODO 

		elif green is not None and swir is not None:
			with rasterio.open(green_filepath) as src_green:
				green = src_green.read(1)
			with rasterio.open(swir_filepath) as src_swir:
				swir = src_swir.read(1)
				# TODO 

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
			with rasterio.open(nir_filepath) as src_nir:
				nir = src_nir.read(1)
			with rasterio.open(swir_filepath) as src_swir:
				swir = src_swir.read(1)
				# TODO 
		else:
			pass


red_file = 'C:\\Users\\eduardo\\Documents\\LT05_L1TP_217076_20110813_20161007_01_T1_B3.TIF'
nir_file = 'C:\\Users\\eduardo\\Documents\\LT05_L1TP_217076_20110813_20161007_01_T1_B4.TIF'
idx = ["ndvi"]
QuickIndex(idx, red = red_file, nir = nir_file)