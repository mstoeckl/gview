/* SPDX-License-Identifier: GPL-3.0-only */
#include "General.hh"

#include <QString>
#include <cmath>

QRgb FColor::rgbaRound() const {
    return qRgba(std::round(v[0] * 255.f), std::round(v[1] * 255.f),
                 std::round(v[2] * 255.f), std::round(v[3] * 255.f));
}

static f3 color_srgb_to_rgb(const f3 &srgb) {
    f3 rgb;
    for (int i = 0; i < 3; i++) {
        rgb[i] = srgb[i] > 0.04045 ? std::pow((srgb[i] + 0.055) / 1.055, 2.4)
                                   : srgb[i] / 12.92;
    }
    return rgb;
}

static f3 color_rgb_to_srgb(const f3 &rgb) {
    f3 srgb;
    for (int i = 0; i < 3; i++) {
        srgb[i] = rgb[i] > 0.04045 / 12.92
                      ? std::pow(rgb[i], 1. / 2.4) * 1.055 - 0.055
                      : 12.92 * rgb[i];
    }
    return srgb;
}
static f3 color_rgb_to_xyz(const f3 &rgb) {
    const float mtx[3][3] = {{0.4124, 0.3576, 0.1805},
                             {0.2126, 0.7152, 0.0722},
                             {0.0193, 0.1192, 0.9505}};

    f3 xyz;
    for (int i = 0; i < 3; i++) {
        xyz[i] = mtx[i][0] * rgb[0] + mtx[i][1] * rgb[1] + mtx[i][2] * rgb[2];
    }
    return xyz;
}
static f3 color_xyz_to_rgb(const f3 &xyz) {
    const float mtx[3][3] = {
        {3.240625477320054, -1.5372079722103185, -0.4986285986982478},
        {-0.9689307147293196, 1.8757560608852413, 0.04151752384295395},
        {0.05571012044551063, -0.20402105059848671, 1.0569959422543882}};
    f3 rgb;
    for (int i = 0; i < 3; i++) {
        rgb[i] = mtx[i][0] * xyz[0] + mtx[i][1] * xyz[1] + mtx[i][2] * xyz[2];
    }
    return rgb;
}
static float lab_f(float x) {
    const float d = 6. / 29;
    return x > d * d * d ? std::cbrt(x) : x / (3 * d * d) + 4. / 29.;
}
static float lab_finv(float y) {
    const float d = 6. / 29;
    return y > d ? y * y * y : (3 * d * d) * (y - 4. / 29);
}
static f3 color_xyz_to_lab(const f3 &xyz) {
    const float d65[3] = {95.047, 100.000, 108.883};
    f3 lab;
    lab[0] = 116. * lab_f(xyz[1] / d65[1]) - 16.0;
    lab[1] = 500. * (lab_f(xyz[0] / d65[0]) - lab_f(xyz[1] / d65[1]));
    lab[2] = 200. * (lab_f(xyz[1] / d65[1]) - lab_f(xyz[2] / d65[2]));
    return lab;
}
static f3 color_lab_to_xyz(const f3 &lab) {
    const float d65[3] = {95.047, 100.000, 108.883};
    f3 xyz;
    xyz[0] = d65[0] * lab_finv((lab[0] + 16.) / 116. + lab[1] / 500);
    xyz[1] = d65[1] * lab_finv((lab[0] + 16.) / 116.);
    xyz[2] = d65[2] * lab_finv((lab[0] + 16.) / 116. - lab[2] / 200);
    return xyz;
}
static f3 color_lab_to_hlc(const f3 &lab) {
    f3 hlc;
    hlc[0] = std::atan2(lab[2], lab[1]);
    hlc[1] = lab[0];
    hlc[2] = std::sqrt(lab[1] * lab[1] + lab[2] * lab[2]);
    return hlc;
}
static f3 color_hlc_to_lab(const f3 &hlc) {
    f3 lab;
    lab[0] = hlc[1];
    lab[1] = hlc[2] * std::cos(hlc[0]);
    lab[2] = hlc[2] * std::sin(hlc[0]);
    return lab;
}
static f3 color_srgb_to_hlc(const f3 &srgb) {
    return color_lab_to_hlc(
        color_xyz_to_lab(color_rgb_to_xyz(color_srgb_to_rgb(srgb))));
}
static f3 color_hlc_to_srgb(const f3 &hlc) {
    return color_rgb_to_srgb(
        color_xyz_to_rgb(color_lab_to_xyz(color_hlc_to_lab(hlc))));
}
float color_srgb_to_nhue(const f3 &srgb) {
    f3 hlc = color_srgb_to_hlc(srgb);
    return (hlc[0] + M_PI) / (M_PI * 2);
}
void seek_max_circle() {
    int N = 100, M = 10000;
    float maxl = 0., maxc = 0.;
    for (int ic = 0; ic <= N; ic++) {
        for (int il = 1; il <= N; il++) {
            // 9.0 is max possible lightness; 16 max chroma
            float fl = 4 + 5. * il * 1. / N, fc = 5. + 5. * ic * 1. / N;
            int nfail = 0;
            for (int ih = 0; ih <= M; ih++) {
                float fh = ih * 1. / M;
                f3 hlc(fh * M_PI * 2 - M_PI, fl, fc);
                f3 srgb = color_hlc_to_srgb(hlc);
                if (srgb[0] > 1.f || srgb[1] > 1.f || srgb[2] > 1.f ||
                    srgb[0] < -0.01f || srgb[1] < -0.01f || srgb[2] < -0.01f) {
                    nfail++;
                    break;
                }
            }
            //            qDebug("l=%f c=%f fail=%d", fl, fc, nfail);
            if (!nfail && fc >= maxc) {
                maxc = fc;
                maxl = fl;
            }
        }
    }
    qDebug("l=%f c=%f", maxl, maxc);
}

f3 rainbow_nhue(const float nhue) {
    static bool r = true;
    if (!r) {
        seek_max_circle();
        r = true;
    }

    // todo: pick lightness for which the chroma circle is widest
    f3 hlc(nhue * M_PI * 2 - M_PI, 4.5, 5.75);
    f3 srgb = color_hlc_to_srgb(hlc);
    for (int i = 0; i < 3; i++) {
        // low-side excess is a rounding error, high-side excess is truncation
        srgb[i] = std::max(srgb[i], 0.f);
        if (srgb[i] > 1.f) {
            qDebug("NHue %f saturated: srgb = %f %f %f", nhue, srgb[0], srgb[1],
                   srgb[2]);
            srgb[i] = 1.f;
        }
    }
    return srgb;
}
