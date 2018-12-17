/* SPDX-License-Identifier: GPL-3.0-only */
#include "Bezier.hh"

#include <algorithm>
#include <cmath>
#include <stdio.h>

static inline double e2(double x) { return x * x; }
static inline double e3(double x) { return x * x * x; }
static inline double e4(double x) { return (x * x) * (x * x); }
static inline double e5(double x) { return (x * x) * (x * x * x); }
static inline double e6(double x) { return (x * x * x) * (x * x * x); }

/**
 * Fast 1d bezier routine.
 * `ts`: length N-2, the arc length prior
 * `vs`: length N, includes first/last coordinates
 */
static void fast_bezier_1d(int N, const double *ts, const double *vs,
                           double ctrl[4]) {
    double b1 = 0., b2 = 0., a11 = 0., a12 = 0., a22 = 0.;
    for (int i = 1; i < N - 1; i++) {
        double t = ts[i - 1];
        double w0 = (1. - t) * (1. - t) * (1. - t);
        double w1 = 3 * (1. - t) * (1. - t) * t;
        double w2 = 3 * (1. - t) * t * t;
        double w3 = t * t * t;

        double r = vs[i] - w0 * vs[0] - w3 * vs[N - 1];
        b1 += r * w1;
        b2 += r * w2;

        // standard amatrix computation
        a11 += w1 * w1;
        a12 += w1 * w2;
        a22 += w2 * w2;
    }

    double det = a11 * a22 - a12 * a12;
    ctrl[0] = vs[0];
    ctrl[3] = vs[N - 1];
    if (det <= 0.) {
        ctrl[1] = vs[N / 3];
        ctrl[2] = vs[2 * N / 3];
    } else {
        ctrl[1] = (a22 * b1 - a12 * b2) / det;
        ctrl[2] = (-a12 * b1 + a11 * b2) / det;
    }
}

/**
 * Given known control points and t-matches, compute the
 * (squared) L2 prediction error.
 */
static double compute_l2_err(int N, const double *xpts, const double *ypts,
                             const double *ts, const double xctrl[4],
                             const double yctrl[4]) {
    double err = 0.;
    for (int i = 0; i < N; i++) {
        double t = ts[i];
        double w0 = (1. - t) * (1. - t) * (1. - t);
        double w1 = 3 * (1. - t) * (1. - t) * t;
        double w2 = 3 * (1. - t) * t * t;
        double w3 = t * t * t;

        double xerr = xpts[i] - w0 * xctrl[0] - w1 * xctrl[1] - w2 * xctrl[2] -
                      w3 * xctrl[3];
        double yerr = ypts[i] - w0 * yctrl[0] - w1 * yctrl[1] - w2 * yctrl[2] -
                      w3 * yctrl[3];
        err += e2(xerr) + e2(yerr);
    }
    return err;
}

/**
 * Given known control points and t-matches, compute the L-infinity prediction
 * error.
 */
static double compute_linf_err(int N, const double *xpts, const double *ypts,
                               const double *ts, const double xctrl[4],
                               const double yctrl[4]) {
    double err = 0.;
    for (int i = 0; i < N; i++) {
        double t = ts[i];
        double w0 = (1. - t) * (1. - t) * (1. - t);
        double w1 = 3 * (1. - t) * (1. - t) * t;
        double w2 = 3 * (1. - t) * t * t;
        double w3 = t * t * t;

        double xerr = xpts[i] - w0 * xctrl[0] - w1 * xctrl[1] - w2 * xctrl[2] -
                      w3 * xctrl[3];
        double yerr = ypts[i] - w0 * yctrl[0] - w1 * yctrl[1] - w2 * yctrl[2] -
                      w3 * yctrl[3];
        err = std::max(err, e2(xerr) + e2(yerr));
    }
    return std::sqrt(err);
}

/**
 * Modify `seq` to enforce that it be a monotonic sequence inside [minv, maxv].
 * (Done by interpolating the longest increasing subsequence.)
 *
 * `scratchA` needs size N, `scratchB` needs size N+1.
 */
static void force_monotonic(int N, double *seq, double minv, double maxv,
                            int *scratchA, int *scratchB) {
    for (int i = 0; i < N; i++) {
        seq[i] = std::min(maxv, std::max(minv, seq[i]));
    }

    int len = 0;
    for (int i = 0; i < N; i++) {
        int lo = 1;
        int hi = len;
        while (lo <= hi) {
            int mid = (lo + hi + 1) / 2;
            if (seq[scratchB[mid]] < seq[i]) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        int newL = lo;

        scratchA[i] = scratchB[newL - 1];
        scratchB[newL] = i;
        if (newL > len) {
            len = newL;
        }
    }

    /* we reuse scratchB this time */
    int k = scratchB[len];
    for (int i = len; i-- > 0;) {
        scratchB[i] = k;
        k = scratchA[k];
    }

    /* now we interpolate all points not part of the monotonic sequence */
    for (int x = 0; x < scratchB[0]; x++) {
        seq[x] = (seq[scratchB[0]] - minv) * x / (double)scratchB[0] + minv;
    }
    for (int i = 0; i < len - 1; i++) {
        for (int x = scratchB[i] + 1; x < scratchB[i + 1]; x++) {
            seq[x] = (x - scratchB[i]) /
                         (double)(scratchB[i + 1] - scratchB[i]) *
                         (seq[scratchB[i + 1]] - seq[scratchB[i]]) +
                     seq[scratchB[i]];
        }
    }
    for (int x = scratchB[len - 1] + 1; x < N; x++) {
        seq[x] = (maxv - seq[scratchB[len - 1]]) * (x - scratchB[len - 1]) /
                     (double)(N - 1 - scratchB[len - 1]) +
                 seq[scratchB[len - 1]];
    }
}

/**
 * Compute the first partial t-derivative of the error objective function at
 * specific [t,xv], given bezier control [x0,x1,x2,x3].
 */
static inline double compute_dt(double t, double xv, double x0, double x1,
                                double x2, double x3) {
    return 6 * e5(t) * e2(x0 - 3 * x1 + 3 * x2 - x3) -
           30 * e4(t) * (x0 - 2 * x1 + x2) * (x0 - 3 * x1 + 3 * x2 - x3) +
           12 * e3(t) *
               (2 * (x0 - x1) * (x0 - 3 * x1 + 3 * x2 - x3) +
                3 * e2(x0 - 2 * x1 + x2)) -
           6 * e2(t) *
               (9 * (x0 - x1) * (x0 - 2 * x1 + x2) +
                (x0 - xv) * (x0 - 3 * x1 + 3 * x2 - x3)) +
           6 * t *
               (5 * e2(x0) - 10 * x0 * x1 + 2 * x0 * x2 - 2 * x0 * xv +
                3 * e2(x1) + 4 * x1 * xv - 2 * x2 * xv) -
           6 * e2(x0) + 6 * x0 * x1 + 6 * x0 * xv - 6 * x1 * xv;
}
/**
 * Compute the second partial t-derivative of the error objective function at
 * specific [t,xv], given bezier control [x0,x1,x2,x3].
 */
static inline double compute_dt_dt(double t, double xv, double x0, double x1,
                                   double x2, double x3) {
    return 30 * e4(t) * e2(-x0 + 3 * x1 - 3 * x2 + x3) +
           20 * e3(t) * (6 * x0 - 12 * x1 + 6 * x2) *
               (-x0 + 3 * x1 - 3 * x2 + x3) +
           3 * e2(t) *
               (4 * (-6 * x0 + 6 * x1) * (-x0 + 3 * x1 - 3 * x2 + x3) +
                4 * e2(3 * x0 - 6 * x1 + 3 * x2)) +
           2 * t *
               (3 * (-6 * x0 + 6 * x1) * (3 * x0 - 6 * x1 + 3 * x2) +
                3 * (2 * x0 - 2 * xv) * (-x0 + 3 * x1 - 3 * x2 + x3)) +
           30 * e2(x0) - 60 * x0 * x1 + 12 * x0 * x2 - 12 * x0 * xv +
           18 * e2(x1) + 24 * x1 * xv - 12 * x2 * xv;
}
static inline double compute_dx1(double t, double xv, double x0, double x1,
                                 double x2, double x3) {
    return 6 * t *
           (e5(t) * (-x0 + 3 * x1 - 3 * x2 + x3) +
            3 * e4(t) * (x0 - 2 * x1 + x2) +
            2 * e4(t) * (x0 - 3 * x1 + 3 * x2 - x3) +
            e3(t) * (-10 * x0 + 18 * x1 - 9 * x2 + x3) +
            e2(t) * (10 * x0 - 12 * x1 + 3 * x2 - xv) +
            t * (-5 * x0 + 3 * x1 + 2 * xv) + x0 - xv);
}
static inline double compute_dx2(double t, double xv, double x0, double x1,
                                 double x2, double x3) {
    return 6 * e2(t) *
           (e4(t) * (x0 - 3 * x1 + 3 * x2 - x3) +
            3 * e3(t) * (-x0 + 2 * x1 - x2) +
            e3(t) * (-x0 + 3 * x1 - 3 * x2 + x3) +
            3 * e2(t) * (2 * x0 - 3 * x1 + x2) + t * (-4 * x0 + 3 * x1 + xv) +
            x0 - xv);
}
static inline double compute_dx1_dt(double t, double xv, double x0, double x1,
                                    double x2, double x3) {
    return 6 * e5(t) * (-6 * x0 + 18 * x1 - 18 * x2 + 6 * x3) +
           15 * e4(t) * (6 * x0 - 12 * x1 + 6 * x2) -
           60 * e4(t) * (-x0 + 3 * x1 - 3 * x2 + x3) +
           4 * e3(t) * (-60 * x0 + 108 * x1 - 54 * x2 + 6 * x3) +
           3 * e2(t) * (60 * x0 - 72 * x1 + 18 * x2 - 6 * xv) +
           2 * t * (-30 * x0 + 18 * x1 + 12 * xv) + 6 * x0 - 6 * xv;
}
static inline double compute_dx2_dt(double t, double xv, double x0, double x1,
                                    double x2, double x3) {
    return 6 * e5(t) * (6 * x0 - 18 * x1 + 18 * x2 - 6 * x3) -
           15 * e4(t) * (6 * x0 - 12 * x1 + 6 * x2) +
           30 * e4(t) * (-x0 + 3 * x1 - 3 * x2 + x3) +
           4 * e3(t) * (36 * x0 - 54 * x1 + 18 * x2) +
           3 * e2(t) * (-24 * x0 + 18 * x1 + 6 * xv) +
           2 * t * (6 * x0 - 6 * xv);
}

static void compute_A(int N, const double *ts, double *a11, double *a12,
                      double *a22) {
    *a11 = 0.;
    *a12 = 0.;
    *a22 = 0.;
    for (int i = 0; i < N; i++) {
        double t = ts[i];
        double w1 = 3 * (1. - t) * (1. - t) * t;
        double w2 = 3 * (1. - t) * t * t;
        *a11 += w1 * w1;
        *a12 += w1 * w2;
        *a22 += w2 * w2;
    }
}

static void newton_t_step(int N, const double *xpts, const double *ypts,
                          const double xctrl[4], const double yctrl[4],
                          double *told, double *tnew) {
    for (int i = 0; i < N; i++) {
        double xdt = compute_dt(told[i], xpts[i], xctrl[0], xctrl[1], xctrl[2],
                                xctrl[3]);
        double ydt = compute_dt(told[i], ypts[i], yctrl[0], yctrl[1], yctrl[2],
                                yctrl[3]);
        double dt = xdt + ydt;

        double xddt = compute_dt_dt(told[i], xpts[i], xctrl[0], xctrl[1],
                                    xctrl[2], xctrl[3]);
        double yddt = compute_dt_dt(told[i], ypts[i], yctrl[0], yctrl[1],
                                    yctrl[2], yctrl[3]);
        double ddt = xddt + yddt;

        if (ddt != 0.) {
            tnew[i] = told[i] - dt / ddt;
        } else {
            tnew[i] = told[i];
        }
    }
}

/**
 * Computes the true L2 bezier curve error relative to a set of points,
 * using the `tseed` and `ttmp` inputs as initial seed and scratch buffer.
 *
 * N is the number of points being tested. (If possible, call without
 * endpoints, as those yield zero error anyway.)
 */
static double true_bezier_linf_error(int N, const double *xpts,
                                     const double *ypts, const double xctrl[4],
                                     const double yctrl[4], double *tseed,
                                     double *ttmp) {
    double err0 = compute_l2_err(N, xpts, ypts, tseed, xctrl, yctrl);

    int *scratch0 = new int[N], *scratch1 = new int[N + 1];
    double mulby = 1.;
    for (int iteration = 0; iteration < 10; iteration++) {
        newton_t_step(N, xpts, ypts, xctrl, yctrl, tseed, ttmp);
        force_monotonic(N, ttmp, 0., 1., scratch0, scratch1);

        double err1 = compute_l2_err(N, xpts, ypts, ttmp, xctrl, yctrl);
        if (err1 < err0) {
            std::swap(ttmp, tseed);
            err0 = err1;
            mulby = 1.;
        } else {
            mulby *= 0.5;
        }
    }

    delete[] scratch0;
    delete[] scratch1;

    return compute_linf_err(N, xpts, ypts, tseed, xctrl, yctrl);
    ;
}

static void compute_deriv(int N, const double *xpts, const double *ypts,
                          const double *tvs, const double xctrl[4],
                          const double yctrl[4], double dctrl[4],
                          double *dtval) {
    dctrl[0] = 0.;
    dctrl[1] = 0.;
    dctrl[2] = 0.;
    dctrl[3] = 0.;
    for (int i = 0; i < N; i++) {
        double xdt =
            compute_dt(tvs[i], xpts[i], xctrl[0], xctrl[1], xctrl[2], xctrl[3]);
        double ydt =
            compute_dt(tvs[i], ypts[i], yctrl[0], yctrl[1], yctrl[2], yctrl[3]);
        dtval[i] = xdt + ydt;

        dctrl[0] += compute_dx1(tvs[i], xpts[i], xctrl[0], xctrl[1], xctrl[2],
                                xctrl[3]);
        dctrl[1] += compute_dx2(tvs[i], xpts[i], xctrl[0], xctrl[1], xctrl[2],
                                xctrl[3]);
        dctrl[2] += compute_dx1(tvs[i], ypts[i], yctrl[0], yctrl[1], yctrl[2],
                                yctrl[3]);
        dctrl[3] += compute_dx2(tvs[i], ypts[i], yctrl[0], yctrl[1], yctrl[2],
                                yctrl[3]);
    }
}

/** Returns true on success/nonsingularity */
static inline bool matrix_4x4_inverse(const double A[4][4], double Ainv[4][4]) {
    // todo: replace with specialized symmetric version
    double a00 = A[0][0], a01 = A[0][1], a02 = A[0][2], a03 = A[0][3];
    double a10 = A[1][0], a11 = A[1][1], a12 = A[1][2], a13 = A[1][3];
    double a20 = A[2][0], a21 = A[2][1], a22 = A[2][2], a23 = A[2][3];
    double a30 = A[3][0], a31 = A[3][1], a32 = A[3][2], a33 = A[3][3];

    double det =
        a00 * a11 * a22 * a33 - a00 * a11 * a23 * a32 - a00 * a12 * a21 * a33 +
        a00 * a12 * a23 * a31 + a00 * a13 * a21 * a32 - a00 * a13 * a22 * a31 -
        a01 * a10 * a22 * a33 + a01 * a10 * a23 * a32 + a01 * a12 * a20 * a33 -
        a01 * a12 * a23 * a30 - a01 * a13 * a20 * a32 + a01 * a13 * a22 * a30 +
        a02 * a10 * a21 * a33 - a02 * a10 * a23 * a31 - a02 * a11 * a20 * a33 +
        a02 * a11 * a23 * a30 + a02 * a13 * a20 * a31 - a02 * a13 * a21 * a30 -
        a03 * a10 * a21 * a32 + a03 * a10 * a22 * a31 + a03 * a11 * a20 * a32 -
        a03 * a11 * a22 * a30 - a03 * a12 * a20 * a31 + a03 * a12 * a21 * a30;
    if (std::abs(det) < 1e-100) {
        return false;
    }
    double idet = 1 / det;

    double b00 = a11 * a22 * a33 - a11 * a23 * a32 - a12 * a21 * a33 +
                 a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31;
    double b01 = a10 * a22 * a33 - a10 * a23 * a32 - a12 * a20 * a33 +
                 a12 * a23 * a30 + a13 * a20 * a32 - a13 * a22 * a30;
    double b02 = a10 * a21 * a33 - a10 * a23 * a31 - a11 * a20 * a33 +
                 a11 * a23 * a30 + a13 * a20 * a31 - a13 * a21 * a30;
    double b03 = a10 * a21 * a32 - a10 * a22 * a31 - a11 * a20 * a32 +
                 a11 * a22 * a30 + a12 * a20 * a31 - a12 * a21 * a30;
    double b10 = a01 * a22 * a33 - a01 * a23 * a32 - a02 * a21 * a33 +
                 a02 * a23 * a31 + a03 * a21 * a32 - a03 * a22 * a31;
    double b11 = a00 * a22 * a33 - a00 * a23 * a32 - a02 * a20 * a33 +
                 a02 * a23 * a30 + a03 * a20 * a32 - a03 * a22 * a30;
    double b12 = a00 * a21 * a33 - a00 * a23 * a31 - a01 * a20 * a33 +
                 a01 * a23 * a30 + a03 * a20 * a31 - a03 * a21 * a30;
    double b13 = a00 * a21 * a32 - a00 * a22 * a31 - a01 * a20 * a32 +
                 a01 * a22 * a30 + a02 * a20 * a31 - a02 * a21 * a30;
    double b20 = a01 * a12 * a33 - a01 * a13 * a32 - a02 * a11 * a33 +
                 a02 * a13 * a31 + a03 * a11 * a32 - a03 * a12 * a31;
    double b21 = a00 * a12 * a33 - a00 * a13 * a32 - a02 * a10 * a33 +
                 a02 * a13 * a30 + a03 * a10 * a32 - a03 * a12 * a30;
    double b22 = a00 * a11 * a33 - a00 * a13 * a31 - a01 * a10 * a33 +
                 a01 * a13 * a30 + a03 * a10 * a31 - a03 * a11 * a30;
    double b23 = a00 * a11 * a32 - a00 * a12 * a31 - a01 * a10 * a32 +
                 a01 * a12 * a30 + a02 * a10 * a31 - a02 * a11 * a30;
    double b30 = a01 * a12 * a23 - a01 * a13 * a22 - a02 * a11 * a23 +
                 a02 * a13 * a21 + a03 * a11 * a22 - a03 * a12 * a21;
    double b31 = a00 * a12 * a23 - a00 * a13 * a22 - a02 * a10 * a23 +
                 a02 * a13 * a20 + a03 * a10 * a22 - a03 * a12 * a20;
    double b32 = a00 * a11 * a23 - a00 * a13 * a21 - a01 * a10 * a23 +
                 a01 * a13 * a20 + a03 * a10 * a21 - a03 * a11 * a20;
    double b33 = a00 * a11 * a22 - a00 * a12 * a21 - a01 * a10 * a22 +
                 a01 * a12 * a20 + a02 * a10 * a21 - a02 * a11 * a20;

    Ainv[0][0] = b00 * idet;
    Ainv[1][0] = -b01 * idet;
    Ainv[2][0] = b02 * idet;
    Ainv[3][0] = -b03 * idet;

    Ainv[0][1] = -b10 * idet;
    Ainv[1][1] = b11 * idet;
    Ainv[2][1] = -b12 * idet;
    Ainv[3][1] = b13 * idet;

    Ainv[0][2] = b20 * idet;
    Ainv[1][2] = -b21 * idet;
    Ainv[2][2] = b22 * idet;
    Ainv[3][2] = -b23 * idet;

    Ainv[0][3] = -b30 * idet;
    Ainv[1][3] = b31 * idet;
    Ainv[2][3] = -b32 * idet;
    Ainv[3][3] = b33 * idet;
    return true;
}

/**
 * Determine the best direction based on the Newton heuristic.
 * (Note: many opportunities for work reduction, thanks to symmetry, etc.)
 */
static void newton_xt_direction(int N, const double *xpts, const double *ypts,
                                const double *tvs, const double xctrl[4],
                                const double yctrl[4], double dctrl[4],
                                double *dtval) {
    /* compute x,t derivative vector */
    double *treg = new double[N];
    double areg[4];
    compute_deriv(N, xpts, ypts, tvs, xctrl, yctrl, areg, treg);

    /* compute hessian components */
    double a11, a12, a22;
    compute_A(N, tvs, &a11, &a12, &a22);
    /* double, since only for one variable at a time*/
    a11 *= 2, a12 *= 2, a22 *= 2;
    // control order is: [x1,x2,y1,y2]
    double E[4][4] = {
        {a11, a12, 0, 0}, {a12, a22, 0, 0}, {0, 0, a11, a12}, {0, 0, a12, a22}};

    double *bbx1 = new double[N];
    double *bbx2 = new double[N];
    double *bby1 = new double[N];
    double *bby2 = new double[N];
    double *tinv = new double[N];
    for (int i = 0; i < N; i++) {
        double xtder = compute_dt_dt(tvs[i], xpts[i], xctrl[0], xctrl[1],
                                     xctrl[2], xctrl[3]);
        double ytder = compute_dt_dt(tvs[i], ypts[i], yctrl[0], yctrl[1],
                                     yctrl[2], yctrl[3]);
        double tder = xtder + ytder;
        tinv[i] = tder == 0. ? 0. : 1 / tder;

        bbx1[i] = compute_dx1_dt(tvs[i], xpts[i], xctrl[0], xctrl[1], xctrl[2],
                                 xctrl[3]);
        bbx2[i] = compute_dx2_dt(tvs[i], xpts[i], xctrl[0], xctrl[1], xctrl[2],
                                 xctrl[3]);
        bby1[i] = compute_dx1_dt(tvs[i], ypts[i], yctrl[0], yctrl[1], yctrl[2],
                                 yctrl[3]);
        bby2[i] = compute_dx2_dt(tvs[i], ypts[i], yctrl[0], yctrl[1], yctrl[2],
                                 yctrl[3]);
    }

    /* apply specialized hessian inverse to derivative */
    /*
     * E = A - (BB.T * Tinv * BB)
     * BI = -Einv @ [BB*tinv]
     * fixi = tinv * bblock @ Einv @ bblock.T @ tinv * treg
     *
     * dxval = - [ Einv @ dx + BI @ dt ]
     * dxval = - [ BI.T @ dx + tinv @ dt + fixi ]
     */

    if (0) {
        for (int i = 0; i < 4; i++) {
            fprintf(stderr, "E: %f\t%f\t%f\t%f\n", E[i][0], E[i][1], E[i][2],
                    E[i][3]);
        }
    }

    const double *bb[4] = {bbx1, bbx2, bby1, bby2};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            // might be able to save runtime by symmetry.
            // (similarly with the inverse.)
            for (int k = 0; k < N; k++) {
                E[i][j] -= bb[i][k] * tinv[k] * bb[j][k];
            }
        }
    }

    double Einv[4][4];
    if (!matrix_4x4_inverse(E, Einv)) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                E[i][j] = 0;
            }
        }
    }

    // dxval = -Einv @ dx
    for (int i = 0; i < 4; i++) {
        dctrl[i] = 0.;
        for (int j = 0; j < 4; j++) {
            dctrl[i] -= Einv[i][j] * areg[j];
        }
    }
    // -= BI@dt ; Einv @ [BB*tinv] @ dt
    for (int i = 0; i < 4; i++) {
        double s = 0.;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < N; k++) {
                s += Einv[i][j] * bb[j][k] * tinv[k] * treg[k];
            }
        }
        dctrl[i] += s;
    }
    // dtval = -BI.T @ dx ; tinv*BB @ Einv @ dx
    for (int k = 0; k < N; k++) {
        double s = 0.;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                s += tinv[k] * bb[i][k] * Einv[i][j] * areg[j];
            }
        }
        s -= tinv[k] * treg[k];
        dtval[k] = s;
    }
    // dtval -= fixi; fixi = tinv * bblock @ Einv @ bblock.T @ tinv * treg
    // (a.k.a, filtering through a C[4] = Einv @ bblock.T @ tinv * treg)
    double C[4];
    for (int i = 0; i < 4; i++) {
        double s = 0.;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < N; k++) {
                s += Einv[i][j] * bb[j][k] * tinv[k] * treg[k];
            }
        }
        C[i] = s;
    }
    for (int k = 0; k < N; k++) {
        double s = 0.;
        for (int i = 0; i < 4; i++) {
            s += bb[i][k] * C[i];
        }
        dtval[k] -= s * tinv[k];
    }

    if (0) {
        // gradient descent example...
        compute_deriv(N, xpts, ypts, tvs, xctrl, yctrl, dctrl, dtval);
        for (int i = 0; i < 4; i++) {
            dctrl[i] = -dctrl[i];
        }
        for (int i = 0; i < N; i++) {
            dtval[i] = -dtval[i];
        }
    }

    delete[] bbx1;
    delete[] bbx2;
    delete[] bby1;
    delete[] bby2;
    delete[] tinv;
    delete[] treg;
}

/**
 * Compute derivative a given distance along a ray, in the direction of the
 * ray.
 */
static double
query_directional_derivative(int N, const double *xpts, const double *ypts,
                             const double *tvs, const double xctrl[4],
                             const double yctrl[4], const double dxval[4],
                             const double *dtval, double offset) {
    double mxctrl[4] = {xctrl[0], xctrl[1] + dxval[0] * offset,
                        xctrl[2] + dxval[1] * offset, xctrl[3]};
    double myctrl[4] = {yctrl[0], yctrl[1] + dxval[2] * offset,
                        yctrl[2] + dxval[3] * offset, yctrl[3]};

    double direct_deriv = 0.;
    for (int i = 0; i < N; i++) {
        double t = tvs[i] + dtval[i] * offset;
        double xdt =
            compute_dt(t, xpts[i], mxctrl[0], mxctrl[1], mxctrl[2], mxctrl[3]);
        double ydt =
            compute_dt(t, ypts[i], myctrl[0], myctrl[1], myctrl[2], myctrl[3]);
        direct_deriv += (xdt + ydt) * dtval[i];

        double dx1 =
            compute_dx1(t, xpts[i], mxctrl[0], mxctrl[1], xctrl[2], mxctrl[3]);
        double dx2 =
            compute_dx2(t, xpts[i], mxctrl[0], mxctrl[1], mxctrl[2], mxctrl[3]);
        double dy1 =
            compute_dx1(t, ypts[i], myctrl[0], myctrl[1], myctrl[2], myctrl[3]);
        double dy2 =
            compute_dx2(t, ypts[i], myctrl[0], myctrl[1], myctrl[2], myctrl[3]);

        direct_deriv += dx1 * dxval[0];
        direct_deriv += dx2 * dxval[1];
        direct_deriv += dy1 * dxval[2];
        direct_deriv += dy2 * dxval[3];
    }

    return direct_deriv;
}

/**
 * Minimize the objective function along a ray, and return the minimal
 * value.
 */
static double directional_min(int N, const double *xpts, const double *ypts,
                              const double *tvs, const double xctrl[4],
                              const double yctrl[4], const double dxval[4],
                              const double *dtval, double *s,
                              double *t_scratch) {
    for (int i = 0; i < N - 2; i++) {
        t_scratch[i] = tvs[i];
    }

    double lower_limit = 0.;
    double upper_limit = 1.;
    double upper_deriv = -std::numeric_limits<double>::infinity();
    while (upper_limit < 10. && upper_deriv < 0) {
        upper_deriv = query_directional_derivative(
            N, xpts, ypts, tvs, xctrl, yctrl, dxval, dtval, upper_limit);
        upper_limit *= 2;
    }
    upper_limit /= 2;

    for (int i = 0; i < 16; i++) {
        double mid_limit = (lower_limit + upper_limit) / 2;
        double mid_deriv = query_directional_derivative(
            N, xpts, ypts, tvs, xctrl, yctrl, dxval, dtval, mid_limit);
        if (mid_deriv > 0) {
            upper_limit = mid_limit;
        } else {
            lower_limit = mid_limit;
        }
    }

    /* select the average point of the two */
    double target = (lower_limit + upper_limit) / 2;

    double mxctrl[4] = {xctrl[0], xctrl[1] + dxval[0] * target,
                        xctrl[2] + dxval[1] * target, xctrl[3]};
    double myctrl[4] = {yctrl[0], yctrl[1] + dxval[2] * target,
                        yctrl[2] + dxval[3] * target, yctrl[3]};

    *s = target;

    // Q: do we need to indirect?
    for (int i = 0; i < N; i++) {
        t_scratch[i] = tvs[i] + dtval[i] * target;
    }
    return compute_l2_err(N, xpts, ypts, t_scratch, mxctrl, myctrl);
}

/**
 * Fit bezier coordinates to interpolate a point path, and calculate error.
 */
double bezier_segment_fit(const double *points, int N, double xctrl[4],
                          double yctrl[4]) {
    /* Start and end control points fixed */
    double *xpts = new double[N];
    double *ypts = new double[N];
    for (int i = 0; i < N; i++) {
        xpts[i] = points[2 * i];
        ypts[i] = points[2 * i + 1];
    }
    double ss = 0.;
    double *arc_length = new double[N - 2];
    for (int i = 1; i < N - 1; i++) {
        ss += std::sqrt(e2(xpts[i] - xpts[i - 1]) + e2(ypts[i] - ypts[i - 1]));
        arc_length[i - 1] = ss;
    }
    ss += std::sqrt(e2(xpts[N - 1] - xpts[N - 2]) +
                    e2(ypts[N - 1] - ypts[N - 2]));
    for (int i = 0; i < N - 2; i++) {
        arc_length[i] /= ss;
    }

    fast_bezier_1d(N, arc_length, xpts, xctrl);
    fast_bezier_1d(N, arc_length, ypts, yctrl);

    /* Newton refinement. tseed, ctrl are known good; ttmp/dctrl are candidates
     */
    double *tseed = new double[N - 2];
    for (int i = 0; i < N - 2; i++) {
        tseed[i] = arc_length[i];
    }
    double *dt = new double[N - 2];
    double *ttmp = new double[N - 2];
    double err = compute_l2_err(N - 2, xpts + 1, ypts + 1, tseed, xctrl, yctrl);
    for (int iteration = 0; iteration < 100; iteration++) {
        double dctrl[4];
        newton_xt_direction(N - 2, xpts + 1, ypts + 1, tseed, xctrl, yctrl,
                            dctrl, dt);
        double s = 0.;
        double n_err = directional_min(N - 2, xpts + 1, ypts + 1, tseed, xctrl,
                                       yctrl, dctrl, dt, &s, ttmp);
        if (n_err < err) {
            xctrl[1] += dctrl[0] * s;
            xctrl[2] += dctrl[1] * s;
            yctrl[1] += dctrl[2] * s;
            yctrl[2] += dctrl[3] * s;
            for (int i = 0; i < N - 2; i++) {
                tseed[i] += dt[i] * s;
            }
            err = n_err;
        } else {
            /* have reached noise threshold or other steady state */
            break;
        }
    }

    /* Final error estimation; we reuse the latest good T values */
    double error = true_bezier_linf_error(N - 2, xpts + 1, ypts + 1, xctrl,
                                          yctrl, tseed, ttmp);
    delete[] dt;
    delete[] ttmp;
    delete[] tseed;
    delete[] xpts;
    delete[] ypts;
    delete[] arc_length;
    return error;
}
