var privatelib = function () {

};


privatelib.vincentyBackward = function (phi1, phi2, L) {
    var a = 6378137;
    var b = 6356752.3142;
    var f = 1 / 298.257223563;
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
    var fwdAz = Math.atan2(cosU2 * Math.sin(lambda), cosU1 * sinU2 - sinU1 * cosU2 * Math.cos(lambda));
    var revAz = Math.atan2(cosU1 * Math.sin(lambda), -sinU1 * cosU2 + cosU1 * sinU2 * Math.cos(lambda));
    return [s, fwdAz, revAz];
}

privatelib.vincentyForward = function (phi1, lambda1, fwdAz, s) {
    var a = 6378137;
    var b = 6356752.3142;
    var f = 1 / 298.257223563;

    var sinfwdAz = Math.sin(fwdAz);
    var cosfwdAz = Math.cos(fwdAz);

    var tanU1 = (1 - f) * Math.tan(phi1), cosU1 = 1 / Math.sqrt((1 + tanU1 * tanU1)), sinU1 = tanU1 * cosU1;
    var sigma1 = Math.atan2(tanU1, cosfwdAz);
    var sinAlpha = cosU1 * sinfwdAz;
    var cosSqAlpha = 1 - sinAlpha * sinAlpha;
    var uSq = cosSqAlpha * (a * a - b * b) / (b * b);
    var A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
    var B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));

    var Sigma = s / (b * A), Sigmat;
    do {
        var cos2SigmaM = Math.cos(2 * sigma1 + Sigma);
        var sinSigma = Math.sin(Sigma);
        var cosSigma = Math.cos(Sigma);
        var DSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
            B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
        Sigmat = Sigma;
        Sigma = s / (b * A) + DSigma;
    } while (Math.abs(Sigma - Sigmat) > 1e-12);

    var tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosfwdAz;
    var phi2 = Math.atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosfwdAz, (1 - f) * Math.sqrt(sinAlpha * sinAlpha + tmp * tmp));
    var lambda = Math.atan2(sinSigma * sinfwdAz, cosU1 * cosSigma - sinU1 * sinSigma * cosfwdAz);
    var C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
    var L = lambda - (1 - C) * f * sinAlpha *
        (Sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
    var lambda2 = (lambda1 + L + 3 * Math.PI) % (2 * Math.PI) - Math.PI;  // normalise to -180...+180

    var revAz = Math.atan2(sinAlpha, -tmp);
    return [phi2, lambda2, revAz];
}

privatelib.extentArea = function (extent) {
    var a = 6378137;
    var b = 6356752.3142;
    var f = 1 / 298.257223563;
    var e = 0.0818191908426215;
    var e2 = Math.pow(0.0818191908426215, 2);
    var phi1 = extent[1];
    var phi2 = extent[3];
    var lam1 = extent[0];
    var lam2 = extent[2];
    var p1 =   Math.sin(phi1) / (1 - e2 * Math.sin(phi1) * Math.sin(phi1)) + Math.log((1 + e * Math.sin(phi1)) / (1 - e * Math.sin(phi1))) / (2 * e);
    var p2 =   Math.sin(phi2) / (1 - e2 * Math.sin(phi2) * Math.sin(phi2)) + Math.log((1 + e * Math.sin(phi2)) / (1 - e * Math.sin(phi2))) / (2 * e);
    var p = 0.5 * a * a * (1 - e2) * (lam2 - lam1) * (p2 - p1);
    return Math.round(p * 1000) / 1000;
};


module.exports = privatelib;
