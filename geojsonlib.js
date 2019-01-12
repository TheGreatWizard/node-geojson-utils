var prl = require('./privatelib');
var geojsonlib = function () {

};


geojsonlib.validate = function (geometry) {
    if (geometry.type === undefined) return false;
    if (geometry.coordinates === undefined) return false;
    if (!Array.isArray(geometry.coordinates)) return false;

    if (!["Point", "MultiPoint", "LineString", "MultiLineString",
        "Polygon", "MultiPolygon", "GeometryCollection"].includes(geometry.type)) return false;

    return true;
};

geojsonlib.isPolygon = function (geometry) {
    return geometry.type === "Polygon";
};
geojsonlib.isLineString = function (geometry) {
    return geometry.type === "LineString";
};
geojsonlib.isPoint = function (geometry) {
    return geometry.type === "Point";
};

// properties
geojsonlib.extentArea = function (extent) {
    var a = 6378137;
    var b = 6356752.3142;
    var f = 1 / 298.257223563;
    var e = 0.0818191908426215;
    var e2 = Math.pow(0.0818191908426215, 2);
    var phi1 = extent[1] * 0.017453292519943295769236908;
    var phi2 = extent[3] * 0.017453292519943295769236908;
    var lam1 = extent[0] * 0.017453292519943295769236908;
    var lam2 = extent[2] * 0.017453292519943295769236908;
    var p1 = Math.sin(phi1) / (1 - e2 * Math.sin(phi1) * Math.sin(phi1)) + Math.log((1 + e * Math.sin(phi1)) / (1 - e * Math.sin(phi1))) / (2 * e);
    var p2 = Math.sin(phi2) / (1 - e2 * Math.sin(phi2) * Math.sin(phi2)) + Math.log((1 + e * Math.sin(phi2)) / (1 - e * Math.sin(phi2))) / (2 * e);
    var p = 0.5 * a * a * (1 - e2) * (lam2 - lam1) * (p2 - p1);
    return Math.round(p * 1000) / 1000;
};

/**
 * returns the area under edge(a geodesic) and parallel defined by the given latitude value
 */
geojsonlib.areaUnderEdge = function (edge, latitude = 0, precision = 1000000) {
    var areaUnderGeodesic = function (phi0, phi1, lam1, phi2, lam2, alpha1, l, p) {
        //console.log((lam2 - lam1) * 180 / Math.PI);
        var a = 6378137;
        var b = 6356752.3142;
        var maxArea = prl.extentArea([lam1, phi0, lam2, phi2]);
        var minArea = prl.extentArea([lam1, phi0, lam2, phi1]);
        if (maxArea < minArea) {
            throw new Error('calculation error area');
        }
        if (lam2 < lam1) {
            throw new Error('calculation error lam');
        }
        var deltaArea = maxArea - minArea;
        //(deltaArea / 3 < p) ||
        if (
            ( (lam2 - lam1) < precision * Math.PI / 180)) {
            //var bestArea = (maxArea + minArea) / 2;

            var acosA2 = Math.pow(a * Math.cos(phi1), 2);
            var bsinA2 = Math.pow(b * Math.sin(phi1), 2);
            var NA = Math.pow(a, 2) / Math.sqrt(acosA2 + bsinA2);
            var acosB2 = Math.pow(a * Math.cos(phi2), 2);
            var bsinB2 = Math.pow(b * Math.sin(phi2), 2);
            var NB = Math.pow(a, 2) / Math.sqrt(acosB2 + bsinB2);
            var bestArea = minArea + deltaArea * NA * Math.cos(phi1) / (NA * Math.cos(phi1) + NB * Math.cos(phi2));


            return bestArea;
        }
        else {
            d = prl.vincentyForward(phi1, lam1, alpha1, l / 2);
            var phiM = d[0];
            var lamM = d[1];
            var alpha2 = d[2];
            if (!((phi1 < phiM) && (phiM < phi2))) {
                throw new Error('calculation error');
            }
            var A1 = areaUnderGeodesic(phi0, phi1, lam1, phiM, lamM, alpha1, l / 2, p );
            var A2 = areaUnderGeodesic(phi0, phiM, lamM, phi2, lam2, alpha2, l / 2, p );
            return A1 + A2;
        }
    };
    var phi1 = edge[0][1] * Math.PI / 180;
    var lam1 = edge[0][0] * Math.PI / 180;
    var phi2 = edge[1][1] * Math.PI / 180;
    var lam2 = edge[1][0] * Math.PI / 180;
    var vb = prl.vincentyBackward(phi1, phi2, lam2 - lam1);
    return areaUnderGeodesic(latitude, phi1, lam1, phi2, lam2, vb[1], vb[0], precision);
}


geojsonlib.polygonArea = function (polygon) {
    throw new Error('method is not implemented');
};

geojsonlib.polygonPerimeter = function (polygon) {
    var s = 0;
    for (var k = 0; k < polygon.coordinates.length; k++) {
        for (var t = 0; t < (polygon.coordinates[k].length - 1); t++) {
            s += this.distanceCoordsToCoords(polygon.coordinates[k][t], polygon.coordinates[k][t + 1]);
        }
    }
    return Math.round(s * 1000) / 1000;
};

geojsonlib.lineStringLength = function (lineString) {
    var s = 0;
    for (var k = 0; k < (lineString.coordinates.length - 1); k++) {
        s += this.distanceCoordsToCoords(lineString.coordinates[k], lineString.coordinates[k + 1]);
    }
    return Math.round(s * 1000) / 1000;
};

geojsonlib.distanceCoordsToCoords = function (coordsA, coordsB) {
    // WGS-84 ellipsiod
    var a = 6378137;
    var b = 6356752.3142;
    var f = 1 / 298.257223563;
    var phi1 = coordsA[1] * 0.017453292519943295769236908;
    var phi2 = coordsB[1] * 0.017453292519943295769236908;
    var L = (coordsB[0] - coordsA[0]) * 0.017453292519943295769236908;

    if (Math.abs(coordsB[1] - coordsA[1]) < 2e-7 && Math.abs(coordsB[0] - coordsA[0]) < 2e-7) {

        var dphi = (coordsB[1] - coordsA[1]) * 0.017453292519943295769236908;

        var phim = (phi1 + phi2) / 2;
        var acos2 = Math.pow(a * Math.cos(phim), 2);
        var bsin2 = Math.pow(b * Math.sin(phim), 2);

        var M1 = Math.pow(a * b, 2) / Math.pow(acos2 + bsin2, 3 / 2);
        var N1 = Math.pow(a, 2) / Math.sqrt(acos2 + bsin2);
        var theta = Math.sqrt(Math.pow(L, 2) + Math.pow(dphi, 2));
        var sinAlpha = L / theta;
        var cosAlpha = dphi / theta;
        var ro = (M1 * N1) / (N1 * Math.pow(cosAlpha, 2) + M1 * Math.pow(sinAlpha, 2));
        return Math.round(ro * theta * 1000) / 1000
    } else {
        var U1 = Math.atan((1 - f) * Math.tan(phi1));
        var U2 = Math.atan((1 - f) * Math.tan(phi2));
        var lambda = L;
        var lambdat = 2 * Math.PI;
        var sinU1 = Math.sin(U1), cosU1 = Math.cos(U1);
        var sinU2 = Math.sin(U2), cosU2 = Math.cos(U2);

        var iterLimit = 20;
        while (Math.abs(lambda - lambdat) > 1e-12 && --iterLimit > 0) {
            var sinLambda = Math.sin(lambda), cosLambda = Math.cos(lambda);
            var sinSigma = Math.sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) +
                (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
            if (sinSigma === 0) return 0;
            var cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
            var sigma = Math.atan2(sinSigma, cosSigma);
            var sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
            var cosSqAlpha = 1 - sinAlpha * sinAlpha;
            var cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
            if (isNaN(cos2SigmaM)) cos2SigmaM = 0;
            var C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
            lambdat = lambda;
            lambda = L + (1 - C) * f * sinAlpha *
                (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
        }
        if (iterLimit === 0) return null;

        var uSq = cosSqAlpha * (a * a - b * b) / (b * b);
        var A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
        var B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
        var deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
            B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
        var s = b * A * (sigma - deltaSigma);

        return Math.round(s * 1000) / 1000;
    }
};


/*
https://www.movable-type.co.uk/scripts/latlong-vincenty.html for reverse
 */

geojsonlib.distanceCoordsToPoint = function (coords, point) {
    return this.distanceCoordsToCoords(coords, point.coordinates);
};
geojsonlib.distanceCoordsToEdge = function (coords, edge) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsToCoordsArray = function (coords, coordsArray) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsToEdgeArray = function (coords, edgeArray) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsToLineString = function (coords, lineString) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsToExtent = function (coords, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsToPolygon = function (coords, polygon) {
    throw new Error('method is not implemented');
};
// distance - from point
geojsonlib.distancePointToPoint = function (pointA, pointB) {
    return this.distanceCoordsToCoords(pointA.coordinates, pointB.coordinates);
};
geojsonlib.distancePointToEdge = function (point, edge) {
    throw new Error('method is not implemented');
};
geojsonlib.distancePointToCoordsArray = function (point, coordsArray) {
    throw new Error('method is not implemented');
};
geojsonlib.distancePointToEdgeArray = function (point, edgeArray) {
    throw new Error('method is not implemented');
};
geojsonlib.distancePointToLineString = function (point, lineString) {
    throw new Error('method is not implemented');
};
geojsonlib.distancePointToExtent = function (point, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.distancePointToPolygon = function (point, polygon) {
    throw new Error('method is not implemented');
};

// distance - from edge
geojsonlib.distanceEdgeToEdge = function (edgeA, edgeB) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceEdgeToCoordsArray = function (edge, coordsArray) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceEdgeToEdgeArray = function (edge, edgeArray) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceEdgeToLineString = function (edge, lineString) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceEdgeToExtent = function (edge, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceEdgeToPolygon = function (edge, polygon) {
    throw new Error('method is not implemented');
};

// distance - from coordsArray

geojsonlib.distanceCoordsArrayToCoordsArray = function (coordsArrayA, coordsArrayB) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsArrayToEdgeArray = function (coordsArray, edgeArray) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsArrayToLineString = function (coordsArray, lineString) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsArrayToExtent = function (coordsArray, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceCoordsArrayToPolygon = function (coordsArray, polygon) {
    throw new Error('method is not implemented');
};

// distance - from edgeArray

geojsonlib.distanceEdgeArrayToEdgeArray = function (edgeArrayA, edgeArrayB) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceEdgeArrayToLineString = function (edgeArray, lineString) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceEdgeArrayToExtent = function (edgeArray, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceEdgeArrayToPolygon = function (edgeArray, polygon) {
    throw new Error('method is not implemented');
};

// distance - from lineString

geojsonlib.distanceLineStringToLineString = function (lineStringA, lineStringB) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceLineStringToExtent = function (lineString, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceLineStringToPolygon = function (lineString, polygon) {
    throw new Error('method is not implemented');
};

// distance - from Extent and polygon

geojsonlib.distanceExtentToExtent = function (extentA, extentB) {
    throw new Error('method is not implemented');
};
geojsonlib.distanceExtentToPolygon = function (extent, polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.distancePolygonToPolygon = function (polygonA, polygonB) {
    throw new Error('method is not implemented');
};


// contains
geojsonlib.isPointInPolygon = function (point, polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.isLineStringInPolygon = function (linestring, polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.isPolygonInPolygon = function (polygonA, polygonB) {
    throw new Error('method is not implemented');
};
geojsonlib.isPointInExtent = function (point, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.isLineStringInExtent = function (linestring, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.isPolygonInExtent = function (polygon, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.isExtentInPolygon = function (extent, polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.isExtentInExtent = function (extentA, extentB) {
    throw new Error('method is not implemented');
};

// intersects
geojsonlib.isLineStringIntersectsLineString = function (linestringA, linestringB) {
    throw new Error('method is not implemented');
};
geojsonlib.isLineStringIntersectsPolygon = function (linestring, polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.isPolygonIntersectsPolygon = function (polygonA, polygonB) {
    throw new Error('method is not implemented');
};
geojsonlib.isLineStringIntersectsExtent = function (linestring, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.isPolygonIntersectsExtent = function (polygon, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.isExtentIntersectsExtent = function (extentA, extentB) {
    throw new Error('method is not implemented');
};

// touches
geojsonlib.isLineStringTouchesLineString = function (linestringA, linestringB, epsilon = 10e-7) {
    throw new Error('method is not implemented');
};
geojsonlib.isLineStringTouchesPolygon = function (linestring, polygon, epsilon = 10e-7) {
    throw new Error('method is not implemented');
};
geojsonlib.isPolygonTouchesPolygon = function (polygonA, polygonB, epsilon = 10e-7) {
    throw new Error('method is not implemented');
};
geojsonlib.isLineStringTouchesExtent = function (linestring, extent, epsilon = 10e-7) {
    throw new Error('method is not implemented');
};
geojsonlib.isPolygonTouchesExtent = function (polygon, extent, epsilon = 10e-7) {
    throw new Error('method is not implemented');
};
geojsonlib.isExtentTouchesExtent = function (extentA, extentB, epsilon = 10e-7) {
    throw new Error('method is not implemented');
};

// conversions
geojsonlib.pointToCoords = function (point) {
    throw new Error('method is not implemented');
};
geojsonlib.lineStringToPolygon = function (linestring) {
    throw new Error('method is not implemented');
};
geojsonlib.lineStringToExtent = function (linestring) {
    throw new Error('method is not implemented');
};
geojsonlib.lineStringToCoordsArray = function (linestring) {
    throw new Error('method is not implemented');
};
geojsonlib.lineStringToEdgeArray = function (linestring) {
    throw new Error('method is not implemented');
};
geojsonlib.lineStringToEdge = function (linestring) {
    throw new Error('method is not implemented');
};


geojsonlib.polygonToLineString = function (polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.polygonToExtent = function (polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.polygonToCoordsArray = function (polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.polygonToEdgeArray = function (polygon) {
    throw new Error('method is not implemented');
};


geojsonlib.coordsArrayToLineString = function (coords) {
    throw new Error('method is not implemented');
};
geojsonlib.coordsArrayToExtent = function (coords) {
    throw new Error('method is not implemented');
};
geojsonlib.coordsArrayToPolygon = function (coords) {
    throw new Error('method is not implemented');
};
geojsonlib.coordsArrayToConvexHull = function (coords) {
    throw new Error('method is not implemented');
};
geojsonlib.coordsArrayToEdgeArray = function (coords) {
    throw new Error('method is not implemented');
};


geojsonlib.edgeArrayToLineString = function (edges) {
    throw new Error('method is not implemented');
};
geojsonlib.edgeArrayToExtent = function (edges) {
    throw new Error('method is not implemented');
};
geojsonlib.edgeArrayToPolygon = function (edges) {
    throw new Error('method is not implemented');
};
geojsonlib.edgeArrayToCoordsArray = function (edges) {
    throw new Error('method is not implemented');
};
geojsonlib.edgeToLineString = function (edge) {
    throw new Error('method is not implemented');
};

// morphology - intersect
geojsonlib.intersectLineStringWithLineString = function (linestringA, linestringB) {
    throw new Error('method is not implemented');
};
geojsonlib.intersectLineStringWithPolygon = function (linestring, polygon) {
    throw new Error('method is not implemented');
};
geojsonlib.intersectPolygonWithPolygon = function (polygonA, polygonB) {
    throw new Error('method is not implemented');
};
geojsonlib.intersectLineStringWithExtent = function (linestring, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.intersectPolygonWithExtent = function (polygon, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.intersectExtentWithExtent = function (extentA, extentB) {
    throw new Error('method is not implemented');
};

// morphology - merge
geojsonlib.mergeLineStringWithLineString = function (linestringA, linestringB) {
    throw new Error('method is not implemented');
};
geojsonlib.mergePolygonWithPolygon = function (polygonA, polygonB) {
    throw new Error('method is not implemented');
};
geojsonlib.mergePolygonWithExtent = function (polygon, extent) {
    throw new Error('method is not implemented');
};
geojsonlib.mergeExtentWithExtent = function (extentA, extentB) {
    throw new Error('method is not implemented');
};

// other
geojsonlib.simplifyPolygon = function (polygon, epsilon = 1e-5) {
    throw new Error('method is not implemented');
};
geojsonlib.simplifyLineString = function (polygon, epsilon = 1e-5) {
    throw new Error('method is not implemented');
};


module.exports = geojsonlib;
