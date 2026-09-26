import os
import numpy as np
from qgis.core import QgsProject, QgsRasterLayer
from qgis import processing
from osgeo import gdal  # needed for Ubuntu 24

# define the output
output_file = os.path.expanduser("~/Soutpansberg_1m_Animal_DSM_Final.tif")

# match layer names
dtm_name = "DTM_5m_Limpopo"
chm_name = "meta_tree_height"

dtm_layer = QgsProject.instance().mapLayersByName(dtm_name)[0]
chm_layer = QgsProject.instance().mapLayersByName(chm_name)[0]

# open high res 
chm_ds = gdal.Open(chm_layer.source())
chm_band = chm_ds.GetRasterBand(1)
chm_array = chm_band.ReadAsArray().astype(np.float32)

# Extract spatial information from the 1m master matrix
geo_transform = chm_ds.GetGeoTransform()
projection = chm_ds.GetProjection()
cols = chm_ds.RasterXSize
rows = chm_ds.RasterYSize

# load ground blocks and align
dtm_ds = gdal.Open(dtm_layer.source())
mem_driver = gdal.GetDriverByName('MEM')
tmp_ds = mem_driver.Create('', cols, rows, 1, gdal.GDT_Float32)
tmp_ds.SetGeoTransform(geo_transform)
tmp_ds.SetProjection(projection)

# Warp the ground matrix grids to match the canopy 1m footprint via bilinear smoothing
gdal.ReprojectImage(dtm_ds, tmp_ds, dtm_ds.GetProjection(), projection, gdal.GRA_Bilinear)
dtm_array = tmp_ds.GetRasterBand(1).ReadAsArray().astype(np.float32)


# strip out corruption
# Safely handles any pixel sitting on background spatial padding (>10000 or <-10000)
dtm_array[(dtm_array > 10000) | (dtm_array < -10000)] = np.nan
chm_array[(chm_array > 10000) | (chm_array < -10000)] = np.nan

# Calculate the actual structural top obstacle surface
dsm_array = dtm_array + chm_array

# Re-map uninitialized margins to a standard clean null marker (-9999)
dsm_array[np.isnan(dsm_array)] = -9999


# explort geotiff into storage
tiff_driver = gdal.GetDriverByName('GTiff')
out_ds = tiff_driver.Create(output_file, cols, rows, 1, gdal.GDT_Float32, options=['COMPRESS=LZW'])
out_ds.SetGeoTransform(geo_transform)
out_ds.SetProjection(projection)
out_band = out_ds.GetRasterBand(1)
out_band.WriteArray(dsm_array)
out_band.SetNoDataValue(-9999)

out_band.FlushCache()
out_ds = None
tmp_ds = None
chm_ds = None
dtm_ds = None

# add to map canvas
final_dsm = QgsRasterLayer(output_file, "Soutpansberg_1m_Animal_DSM_Final")
QgsProject.instance().addMapLayer(final_dsm)
