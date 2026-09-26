var fabdem = ee.ImageCollection("projects/sat-io/open-datasets/FABDEM");
var aoi = ee.Geometry.Rectangle([29.38449, -23.05484, 29.50343, -22.99472]);

Map.centerObject(aoi, 12);
Map.addLayer(fabdem.mosaic().clip(aoi), {min: 800, max: 1600}, 'FABDEM');

Export.image.toDrive({
  image: fabdem.mosaic().clip(aoi),
  description: 'fabdem_soutpansberg',
  scale: 30,
  region: aoi,
  crs: 'EPSG:32735',
  maxPixels: 1e9
});
