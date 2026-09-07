// Cable winch stepper kinematics with flex compensation
//
// Copyright (C) 2018-2019  Kevin O'Connor <kevin@koconnor.net>
// Copyright (C) 2024       Contributors
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h>   // sqrt, fabs, fmax, fmin
#include <stddef.h> // offsetof
#include <stdlib.h> // malloc
#include <stdio.h> // malloc
#include <string.h> // memset
#include <float.h> // DBL_MAX
#include "compiler.h" // __visible
#include "pyhelper.h" // errorf
#include "itersolve.h" // struct stepper_kinematics
#include "trapq.h" // move_get_coord

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define WINCH_MAX_ANCHORS 20
#define EPSILON 1e-9
#define G_ACCEL 9.81
// LAMBDA is Tikhonov regularization weight in both the Tikhonov and the QP solver.
#define LAMBDA 1e-3

enum winch_force_algorithm {
    WINCH_FORCE_ALGO_TIKHONOV = 0,
    WINCH_FORCE_ALGO_QP = 1,
};

struct winch_flex {
    struct coord anchors[WINCH_MAX_ANCHORS];
    int num_anchors;
    int enabled;
    int flex_compensation_algorithm;
    double buildup_factor;
    double mover_weight;
    double spring_constant;
    double min_force[WINCH_MAX_ANCHORS];
    double max_force[WINCH_MAX_ANCHORS];
    double guy_wires[WINCH_MAX_ANCHORS];
    int ignore_gravity;
    int ignore_pretension;
    double distances_origin[WINCH_MAX_ANCHORS];
    int mechanical_advantage[WINCH_MAX_ANCHORS];
    double spool_radius[WINCH_MAX_ANCHORS];
    double spool_radius_sq[WINCH_MAX_ANCHORS];
    double k0[WINCH_MAX_ANCHORS];
    double k2[WINCH_MAX_ANCHORS];
    double steps_per_mm[WINCH_MAX_ANCHORS];
    double inv_steps_per_mm[WINCH_MAX_ANCHORS];
    int use_constant_spool_model[WINCH_MAX_ANCHORS];
};

struct winch_stepper {
    struct stepper_kinematics sk;
    struct winch_flex *wf;
    int index;
    struct coord anchor;
};

static inline double
hypot2(double dx, double dy)
{
    return sqrt(dx*dx + dy*dy);
}

static int
invert2x2(const double M[2][2], double Minv[2][2])
{
    double a = M[0][0], b = M[0][1];
    double c = M[1][0], d = M[1][1];

    double det = a * d - b * c;
    if (fabs(det) < EPSILON)
        return 0;

    double invdet = 1. / det;

    Minv[0][0] =  d * invdet;
    Minv[0][1] = -b * invdet;
    Minv[1][0] = -c * invdet;
    Minv[1][1] =  a * invdet;

    return 1;
}

static int
build_direction_matrix(struct winch_flex *wf, const struct coord *pos,
                       double *A)
{
    int N = wf->num_anchors;
    int valid = 0;

    for (int j = 0; j < N; ++j) {
        double dx = wf->anchors[j].x - pos->x;
        double dy = wf->anchors[j].y - pos->y;
        double norm = hypot2(dx, dy);

        if (norm < EPSILON) {
            A[0 * N + j] = 0.;
            A[1 * N + j] = 0.;
            continue;
        }

        double inv = 1.0 / norm;
        A[0 * N + j] = dx * inv;
        A[1 * N + j] = dy * inv;
        valid++;
    }

    return valid >= 2;
}

static void
solve_min_norm_T(const double *A, int N, const double Fext[2], double lambda,
                 double *T)
{
    double S[2][2] = {
        {lambda, 0.},
        {0., lambda}
    };

    for (int j = 0; j < N; ++j) {
        double ax = A[0 * N + j];
        double ay = A[1 * N + j];

        S[0][0] += ax * ax;
        S[0][1] += ax * ay;
        S[1][0] += ay * ax;
        S[1][1] += ay * ay;
    }

    double Sinv[2][2];

    if (!invert2x2(S, Sinv)) {
        S[0][0] += 1e-6;
        S[1][1] += 1e-6;
        invert2x2(S, Sinv);
    }

    double y0 = Sinv[0][0] * Fext[0]
              + Sinv[0][1] * Fext[1];

    double y1 = Sinv[1][0] * Fext[0]
              + Sinv[1][1] * Fext[1];

    for (int j = 0; j < N; ++j) {
        double ax = A[0 * N + j];
        double ay = A[1 * N + j];

        T[j] = ax * y0 + ay * y1;
    }
}

static void
build_null_projector(const double *A, int N, double lambda, double *P)
{
    double S[2][2] = {
        {lambda, 0.},
        {0., lambda}
    };

    for (int j = 0; j < N; ++j) {
        double ax = A[0 * N + j];
        double ay = A[1 * N + j];

        S[0][0] += ax * ax;
        S[0][1] += ax * ay;
        S[1][0] += ay * ax;
        S[1][1] += ay * ay;
    }

    double Sinv[2][2];

    if (!invert2x2(S, Sinv)) {
        S[0][0] += 1e-6;
        S[1][1] += 1e-6;
        invert2x2(S, Sinv);
    }

    /*
     * P = I - A^T * (A * A^T + lambda I)^-1 * A
     *
     * P is an N x N projector onto the null-space
     * of the 2D cable direction matrix A.
     */
    for (int i = 0; i < N; ++i) {
        double aix = A[0 * N + i];
        double aiy = A[1 * N + i];

        for (int j = 0; j < N; ++j) {
            double ajx = A[0 * N + j];
            double ajy = A[1 * N + j];

            double v0 = Sinv[0][0] * ajx
                      + Sinv[0][1] * ajy;

            double v1 = Sinv[1][0] * ajx
                      + Sinv[1][1] * ajy;

            double value = aix * v0 + aiy * v1;

            P[i * N + j] = (i == j ? 1. : 0.) - value;
        }
    }
}

static void
project_nullspace(const double *P, int N, const double *v, double *out)
{
    for (int r = 0; r < N; ++r) {
        double acc = 0.;
        for (int c = 0; c < N; ++c)
            acc += P[r * N + c] * v[c];
        out[r] = acc;
    }
}
