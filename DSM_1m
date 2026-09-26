from qgis.core import QgsProject, QgsRasterLayer
from qgis import processing

# reference layers
dtm_name = "DTM_5m_Limpopo"
chm_name = "meta_tree_height"

dtm_layer = QgsProject.instance().mapLayersByName(dtm_name)[0]
chm_layer = QgsProject.instance().mapLayersByName(chm_name)[0]


# both at 1-m layout
calc_params = {
    'EXPRESSION': f'"{dtm_name}@1" + "{chm_name}@1"',
    'LAYERS': [dtm_layer, chm_layer],
    'CELLSIZE': 1.0,               # Enforce 1-meter target layout
    'EXTENT': chm_layer.extent(),  # Lock boundary dimensions to the high-res canopy
    'REFERENCE_LAYER': chm_layer,
    'OUTPUT': 'TEMPORARY_OUTPUT'
}

# interpolating calculation
result = processing.run("qgis:rastercalculator", calc_params)
output_path = result['OUTPUT']

# convert to raster
dsm_layer = QgsRasterLayer(output_path, "Soutpansberg_1m_Animal_DSM")

if dsm_layer.isValid():
    QgsProject.instance().addMapLayer(dsm_layer)
