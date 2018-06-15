/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include <QPointF>
#include <QRgb>

inline uint32_t randint(uint32_t excl_upper) {
    return ((uint32_t)qrand()) % excl_upper;
}

// Very fast, simple, and compact color handling
class VColor {
public:
    VColor(uint8_t r, uint8_t g, uint8_t b) : pr(r), pg(g), pb(b) {}
    VColor(QRgb rgb) : pr(qRed(rgb)), pg(qGreen(rgb)), pb(qBlue(rgb)) {}
    static VColor fromRgbF(float r, float g, float b) {
        return VColor(255.f * r, 255.f * g, 255.f * b);
    }
    QRgb rgb() const { return qRgb(pr, pg, pb); }
    float redF() const { return pr / 255.f; }
    float greenF() const { return pg / 255.f; }
    float blueF() const { return pb / 255.f; }

private:
    uint8_t pr, pg, pb;
};

// Fast color handling
class FColor {
public:
    FColor() {
        v[0] = 0.f;
        v[1] = 0.f;
        v[2] = 0.f;
        v[3] = 0.f;
    }
    FColor(float r, float g, float b, float a = 1.0f) {
        v[0] = r;
        v[1] = g;
        v[2] = b;
        v[3] = a;
    }
    explicit FColor(QRgb x) {
        v[0] = qRed(x) / 255.f;
        v[1] = qGreen(x) / 255.f;
        v[2] = qBlue(x) / 255.f;
        v[3] = qAlpha(x) / 255.f;
    }

    QRgb rgba() const {
        return qRgba(v[0] * 255.f, v[1] * 255.f, v[2] * 255.f, v[3] * 255.f);
    }
    inline float redF() const { return v[0]; }
    inline float greenF() const { return v[1]; }
    inline float blueF() const { return v[2]; }
    inline float alphaF() const { return v[3]; }
    float magnitude() const {
        return v[3] * std::max(v[2], std::max(v[1], v[0]));
    }
    QString hexName() const;
    static FColor blend(const FColor &a, const FColor &b, float s) {
        float t = 1. - s;
        return FColor(t * a.v[0] + s * b.v[0], t * a.v[1] + s * b.v[1],
                      t * a.v[2] + s * b.v[2], t * a.v[3] + s * b.v[3]);
    }

private:
    float v[4];
};

inline int ceil_div(int v, int s) { return v / s + ((v % s) > 0 ? 1 : 0); }

/**
 * Grid point specification. If we map the screen pixel array to the viewport,
 * provides center points for subsampled pixels. Supersample pixels are square
 * and are aligned to the (0,0) corner.
 */
class GridSpec {
public:
    GridSpec(int iw, int ih, int pixel_size)
        : w(iw), h(ih), s(pixel_size), imind_t_s_t_2(s * 2. / std::min(w, h)),
          mind_t_is_d_2(1. / imind_t_s_t_2),
          off_ts_x(-(0.5 - 0.5 * ceil_div(w, s))),
          off_ts_y(-(0.5 - 0.5 * ceil_div(h, s))),
          off_tv_x(-imind_t_s_t_2 * off_ts_x),
          off_tv_y(-imind_t_s_t_2 * off_ts_y) {}

    // ceil[w/s], ceil[h/s]
    inline int sampleWidth() const { return ceil_div(w, s); }
    inline int sampleHeight() const { return ceil_div(h, s); }
    inline int imageWidth() const { return w; }
    inline int imageHeight() const { return h; }
    inline int pixelSize() const { return s; }

    // the coordinates for a sample point
    // includes barely the [-1,1]x[-1,1] box
    inline QPointF toViewCoord(int x, int y) const {
        // sample point center in image coordinates
        qreal sx = imind_t_s_t_2 * x + off_tv_x;
        qreal sy = imind_t_s_t_2 * y + off_tv_y;
        return QPointF(sx, sy);
    }

    // Map from box incl. [-1,1]x[-1,1] to sample coordinates
    inline QPointF toSampleCoord(QPointF viewCoord) const {
        qreal x = mind_t_is_d_2 * viewCoord.x() + off_ts_x;
        qreal y = mind_t_is_d_2 * viewCoord.y() + off_ts_y;
        return QPointF(x, y);
    }

private:
    int w, h, s;
    qreal imind_t_s_t_2, mind_t_is_d_2;
    qreal off_ts_x, off_ts_y, off_tv_x, off_tv_y;
};
