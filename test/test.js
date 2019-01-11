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

describe('lineStringLength', function () {
    it('calculate length', function () {
        var line = {
            type: "LineString",
            coordinates: [[70.72897, 123.046875],
                [63.074865, 79.453125],
                [47.517204, 45.3515625],
                [46.316584, 5.9765625]]
        };
        console.log((gsl.distanceCoordsToCoords([70.72897, 123.046875], [63.074865, 79.453125])) +
            (gsl.distanceCoordsToCoords([63.074865, 79.453125], [47.517204, 45.3515625])) +
            (gsl.distanceCoordsToCoords([47.517204, 45.3515625], [46.316584, 5.9765625])));

        console.log('The length of the LineString is', gsl.lineStringLength(line),'m');
        assert.equal(gsl.lineStringLength(line), 7721785.168);
    });
})
;

describe('distanceCoordsToCoords', function () {
    it('correct type', function () {
        var coords1 = [53.000, 2.000];
        var coords2 = [52.000, 1.000];
        assert.equal(typeof gsl.distanceCoordsToCoords(coords1, coords2), 'number')
    });

    it('correct length 1', function () {
        var coords1 = [53.000, 2.000];
        var coords2 = [52.000, 1.000];
        assert.equal(gsl.distanceCoordsToCoords(coords1, coords2), 130359.286);
    });

    it('correct length 2', function () {
        var coords1 = [18.13256954, 13.95846582];
        var coords2 = [-36.52895455, -22.25558497];
        assert.equal(gsl.distanceCoordsToCoords(coords1, coords2), 7148447.418);
    });

    it('short distance timeout', function () {
        this.timeout(600);
        var epsilon = 1e-7;
        var coords1 = [18.13256954, 13.25846582];
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

        var coords1 = [18.13256954, 13.95846582];
        var coords2 = [-36.52895455, -22.25558497];
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
            coordinates: [-50.52695544, 31.95841582]
        };
        var pnt2 = {
            type: "Point",
            coordinates: [81.52895455, 104.27778239]
        };
        assert.equal(gsl.distancePointToPoint(pnt1, pnt2), 15242863.460);
    });

    it('correct length 2', function () {
        var pnt1 = {
            type: "Point",
            coordinates: [45.05269554, 82.31958415]
        };
        var pnt2 = {
            type: "Point",
            coordinates: [55.81528954, 0.00104277]
        };
        assert.equal(gsl.distancePointToPoint(pnt1, pnt2), 5611645.267);
    });
});

describe('distanceCoordsToPoint', function () {
    it('correct length 1', function () {
        var coords = [-20.0, -73.82319584];

        var pnt = {
            type: "Point",
            coordinates: [0.0, 0.0]
        };

        assert.equal(gsl.distanceCoordsToPoint(coords, pnt), 8326509.923);
    });
    it('correct length 2', function () {
        var coords = [10.89843751, 51.17934297];

        var pnt = {
            type: "Point",
            coordinates: [31.9921875, 48.92249926]
        };

        assert.equal(gsl.distanceCoordsToPoint(coords, pnt), 2347247.168);
    });

});


// console.log(Object.keys(gsl).length);
// console.log(Object.keys(gsl));
