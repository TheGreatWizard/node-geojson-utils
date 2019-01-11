# node-geojson-utils
Utility function and spatial operations with geojson geometries

## UNDER CONSTRUCTION

The purpose of this object is to provide a common set of spatial operation functions
for 2D geometries on the earth surface.

1. Geographic coordinates in WGS84 datum.
2. Geojson format
3. Precision 1mm
4. Efficient calculations  

### Methods

1. Conversion - e.g. convert polygon to linestring
2. Spatial relation - e.g. is point inside polygon. are linestring inretect
3. Morphology - e.g. calculate the intersection polygon of two polygons
4. Utility - e.g. calculate polygon extent or convex-hull

##### Supported Geometries
1. Point
2. Polygon
3. LineString

## Example

Calculate LineString length:

    var gsl = require('geojsonlib');
    
    var line = {
        type: "LineString",
        coordinates: [[70.72897, 123.046875],
            [63.074865, 79.453125],
            [47.517204, 45.3515625],
            [46.316584, 5.9765625]]
    };
   
    console.log('The length of the LineString is', gsl.lineStringLength(line),'m');
        
        `
