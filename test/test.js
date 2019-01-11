var assert = require('assert');
var gsl = require('../geojsonlib');
describe('validate', function () {
    it('correct point', function () {
        var pnt = {
            type: "Point",
            coordinates: [48.515625, 45.336701]
        };

        assert.equal(gsl.validate(pnt), true);
    });
    it('no type', function () {
        var pnt = {
            tipe: "Point",
            coordinates: [48.515625, 45.336701]
        };

        assert.equal(gsl.validate(pnt), false);
    });

    it('bad type', function () {
        var pnt = {
            type: "SuperGeomerty",
            coordinates: [48.515625, 45.336701]
        };

        assert.equal(gsl.validate(pnt), false);
    });
    it('no coordinates', function () {
        var pnt = {
            type: "Point",
            cordinates: [48.515625, 45.336701]
        };

        assert.equal(gsl.validate(pnt), false);
    });
    it('bad coordinates', function () {
        var pnt = {
            type: "Point",
            coordinates: "[48.515625, 45.336701]"
        };

        assert.equal(gsl.validate(pnt), false);
    });

});
describe('polygonArea', function () {
    it('calculate geodesic', function () {
        var phi1 = 0;
        var lam1 = 0;
        var phi2 = 10 * Math.PI / 180;
        var lam2 = 10 * Math.PI / 180;
        var n = 1000000;

        d2 = gsl.vincentyBackward(phi1, phi2, lam2 - lam1);
        var step = d2[0] / n;
        var fwdAz = d2[1];
        var area = 0;
        var area2 = 0;
        var lamA = lam1 * 180 / Math.PI;
        var phiA = phi1 * 180 / Math.PI;

        for (var k = 0; k < n; k++) {
            d = gsl.vincentyForward(phi1, lam1, fwdAz, step * (k + 1));
            var phiB = d[0] * 180 / Math.PI;
            var lamB = d[1] * 180 / Math.PI;
           // console.log(d[0] * 180 / Math.PI, d[1] * 180 / Math.PI);
            area += gsl.extentArea([lamA, 0, lamB, phiB]);
            area2 += gsl.extentArea([lamA, 0, lamB, phiA]);
            lamA = lamB;
            phiA = phiB;
        }
        //  122473224427.63501 10
        console.log('Area =', area);
        console.log('Area2 =', area2);
        console.log('delta Area =', area - area2);
        //
        // d = gsl.vincentyForward(0, 0, 45 * Math.PI / 180, 2000);
        // var phi2 = d[0];
        // var lambda2 = d[1];
        // var revAz = d[2];
        // d2 = gsl.vincentyBackward(0, phi2, lambda2 - 0);
        // //var s, fwdAz, revAz];
        // var s = d2[0];
        // var fwdAz = d2[1] * 180 / Math.PI;
        // var revAz2 = d2[2] * 180 / Math.PI;
        // console.log(revAz2);
    });
});
describe('extentArea', function () {
    it('calculate extent area', function () {
        assert.equal(gsl.extentArea([-22, -32, 50, 36]), 5.6835331379361180e+13);
    });
});


describe('polygonPerimeter', function () {
    it('calculate polygon perimeter', function () {
        var polygon = {
            type: "Polygon",
            coordinates: [[[-45.703125, 23.241346],
                [-21.796875, 23.241346], [-21.796875, 41.508577],
                [-45.703125, 41.508577], [-45.703125, 23.241346]]]
        };

        assert.equal(gsl.polygonPerimeter(polygon), 8484727.45);
    });
});

describe('lineStringLength', function () {
    it('calculate length', function () {
        var line = {
            type: "LineString",
            coordinates: [[123.046875, 70.72897],
                [79.453125, 63.074865],
                [45.3515625, 47.517204],
                [5.9765625, 46.316584]]
        };

        assert.equal(gsl.lineStringLength(line), 7721785.168);
    });
})
;

describe('distanceCoordsToCoords', function () {
    it('correct type', function () {
        var coords1 = [2.000, 53.000];
        var coords2 = [1.000, 52.000];
        assert.equal(typeof gsl.distanceCoordsToCoords(coords1, coords2), 'number')
    });

    it('correct length 1', function () {
        var coords1 = [2.000, 53.000];
        var coords2 = [1.000, 52.000];
        assert.equal(gsl.distanceCoordsToCoords(coords1, coords2), 130359.286);
    });

    it('correct length 2', function () {
        var coords1 = [13.95846582, 18.13256954];
        var coords2 = [-22.25558497, -36.52895455];
        assert.equal(gsl.distanceCoordsToCoords(coords1, coords2), 7148447.418);
    });

    it('short distance timeout', function () {
        this.timeout(600);
        var epsilon = 1e-7;
        var coords1 = [13.25846582, 18.13256954];
        var coords2 = [coords1[0] + epsilon, coords1[1] + epsilon];
        var start = new Date().getTime();
        for (var i = 0; i < 1e6; ++i) {
            var len = gsl.distanceCoordsToCoords(coords1, coords2);
        }
        var end = new Date().getTime();
        var time = end - start;
        assert.ok(time < 600);
    });

    it('long distance timeout', function () {

        var coords1 = [13.95846582, 18.13256954];
        var coords2 = [-22.25558497, -36.52895455];
        var start = new Date().getTime();
        for (var i = 0; i < 1e6; ++i) {
            var len = gsl.distanceCoordsToCoords(coords1, coords2);
        }
        var end = new Date().getTime();
        var time = end - start;
        assert.ok(time < 1200);
    });
});

describe('distancePointToPoint', function () {
    it('correct length 1', function () {
        var pnt1 = {
            type: "Point",
            coordinates: [31.95841582, -50.52695544]
        };
        var pnt2 = {
            type: "Point",
            coordinates: [104.27778239, 81.52895455]
        };
        assert.equal(gsl.distancePointToPoint(pnt1, pnt2), 15242863.460);
    });

    it('correct length 2', function () {
        var pnt1 = {
            type: "Point",
            coordinates: [82.31958415, 45.05269554]
        };
        var pnt2 = {
            type: "Point",
            coordinates: [0.00104277, 55.81528954]
        };
        assert.equal(gsl.distancePointToPoint(pnt1, pnt2), 5611645.267);
    });
});

describe('distanceCoordsToPoint', function () {
    it('correct length 1', function () {
        var coords = [-73.82319584, -20.0];

        var pnt = {
            type: "Point",
            coordinates: [0.0, 0.0]
        };

        assert.equal(gsl.distanceCoordsToPoint(coords, pnt), 8326509.923);
    });
    it('correct length 2', function () {
        var coords = [51.17934297, 10.89843751];

        var pnt = {
            type: "Point",
            coordinates: [48.92249926, 31.9921875]
        };

        assert.equal(gsl.distanceCoordsToPoint(coords, pnt), 2347247.168);
    });

});


// console.log(Object.keys(gsl).length);
// console.log(Object.keys(gsl));
