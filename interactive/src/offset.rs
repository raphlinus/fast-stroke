//! An implementation of curve offset.

use kurbo::{
    common::{solve_itp, solve_quadratic},
    BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, QuadBez, Vec2,
};

use crate::cusp::CuspAnalysis;

// Copyright 2022 the Kurbo Authors
// SPDX-License-Identifier: Apache-2.0 OR MIT

/// Info that's constant for the curve.
struct CubicOffset {
    /// Source curve.
    c: CubicBez,
    /// Derivative of source curve.
    q: QuadBez,
    /// Offset.
    d: f64,
    // c0 + c1 t + c2 t^2 is the cross product of second and first
    // derivatives of the underlying cubic, multiplied by offset (for
    // computing cusp).
    c0: f64,
    c1: f64,
    c2: f64,
    tolerance: f64,
    // Maybe make this optional, or make creation cheap in smooth case
    cusp: CuspAnalysis,
    // move &mut result here?
}

/// Info specific to a subdivision.
struct OffsetRec {
    t0: f64,
    t1: f64,
    // unit tangent at t0
    utan0: Vec2,
    // unit tangent at t1
    utan1: Vec2,
    cusp0: f64,
    cusp1: f64,
    depth: usize,
}

// We never let cusp values haven an absolute value smaller than
// this. When a cusp is found, determine its sign and use this value.
const CUSP_EPSILON: f64 = 1e-12;

/// Maximum recursion depth
///
/// Perhaps should be configurable.
const MAX_DEPTH: usize = 8;

// t values for minmax and error estimation
const T1_FOR_OFFSET: f64 = 1. / 6.;
const OFFSET_TS: [f64; 3] = [T1_FOR_OFFSET, 0.5, 1. - T1_FOR_OFFSET];

pub fn offset_cubic(c: CubicBez, d: f64, tolerance: f64) -> BezPath {
    let mut result = BezPath::new();
    let co = CubicOffset::new(c, d, tolerance);
    // TODO: cusp analysis to avoid divide by 0 in following math.
    let (cusp0, utan0) = co.cusp_and_utan(co.q.p0, co.c0);
    let (cusp1, utan1) = co.cusp_and_utan(co.q.p2, co.c0 + co.c1 + co.c2);
    result.move_to(c.p0 + d * turn(utan0));
    let rec = OffsetRec {
        t0: 0.0,
        t1: 1.0,
        utan0,
        utan1,
        cusp0,
        cusp1,
        depth: 0,
    };
    co.offset_rec(&rec, &mut result);
    result
}

impl CubicOffset {
    /// Create a new curve from Bézier segment and offset.
    fn new(c: CubicBez, d: f64, tolerance: f64) -> Self {
        let q = c.deriv();
        let d2 = 2.0 * d;
        let p1xp0 = q.p1.to_vec2().cross(q.p0.to_vec2());
        let p2xp0 = q.p2.to_vec2().cross(q.p0.to_vec2());
        let p2xp1 = q.p2.to_vec2().cross(q.p1.to_vec2());
        let cusp = CuspAnalysis::new(q);
        CubicOffset {
            c,
            q,
            d,
            c0: d2 * p1xp0,
            c1: d2 * (p2xp0 - 2.0 * p1xp0),
            c2: d2 * (p2xp1 - p2xp0 + p1xp0),
            tolerance,
            cusp,
        }
    }

    // Compute a function which has a zero-crossing at cusps, and is
    // positive at low curvatures on the source curve.
    fn cusp_sign(&self, t: f64) -> f64 {
        let ds2 = self.q.eval(t).to_vec2().hypot2();
        ((self.c2 * t + self.c1) * t + self.c0) / (ds2 * ds2.sqrt()) + 1.0
    }

    /// Compute cusp and unit tangent of endpoint.
    ///
    /// Equivalent to computing both separately, but more optimized.
    ///
    /// The y parameter should be c0 for the start point, and c0 + c1 + c2
    /// for the end point; it's just evaluating the polynomial at t=0 and
    /// t=1.
    fn cusp_and_utan(&self, tan: Point, y: f64) -> (f64, Vec2) {
        let rsqrt = 1.0 / tan.to_vec2().hypot();
        let cusp = y * (rsqrt * rsqrt * rsqrt) + 1.0;
        let utan = rsqrt * tan.to_vec2();
        (cusp, utan)
    }

    fn offset_rec(&self, rec: &OffsetRec, result: &mut BezPath) {
        if rec.cusp0 * rec.cusp1 < 0.0 {
            let a = rec.t0;
            let b = rec.t1;
            let s = rec.cusp1.signum();
            let f = |t| s * self.cusp_sign(t);
            let k1 = 0.2 / (b - a);
            const ITP_EPS: f64 = 1e-12;
            let t = solve_itp(f, a, b, ITP_EPS, 1, k1, s * rec.cusp0, s * rec.cusp1);
            // TODO(robustness): If we're unlucky, there will be 3 cusps between t0
            // and t1, and the solver will land on the middle one. In that case, the
            // derivative on cusp will be the opposite sign as expected.
            //
            // If we're even more unlucky, there is a second-order cusp, both zero
            // cusp value and zero derivative.
            let utan_t = self.q.eval(t).to_vec2().normalize();
            let cusp_t_minus = CUSP_EPSILON.copysign(rec.cusp0);
            let cusp_t_plus = CUSP_EPSILON.copysign(rec.cusp1);
            self.subdivide(rec, result, t, utan_t, cusp_t_minus, cusp_t_plus);
            return;
        }
        let utans = self.utans(rec);
        let (a, b) = self.one_point(rec, &utans);
        let mut ts = OFFSET_TS;
        let mut c_approx = self.apply(rec, a, b);
        let err_one_point = self.est_error(rec, &utans, c_approx, &mut ts);
        let mut err = err_one_point;
        // early out if err is in tolerance?
        let (a_minmax, b_minmax) = self.linear_minmax(rec, &utans);
        let c_minmax = self.apply(rec, a_minmax, b_minmax);
        ts = OFFSET_TS;
        let err_minmax = self.est_error(rec, &utans, c_minmax, &mut ts);
        if err_minmax < err {
            err = err_minmax;
            c_approx = c_minmax;
        }
        let (mut a_refined, mut b_refined) = (a_minmax, b_minmax);
        (a_refined, b_refined) =
            self.refine_minmax(rec, &utans, c_approx, a_refined, b_refined, ts);
        let c_refined = self.apply(rec, a_refined, b_refined);
        let err_refined = self.est_error(rec, &utans, c_refined, &mut ts);
        if err_refined < err {
            err = err_refined;
            c_approx = c_refined;
        }
        // web_sys::console::log_1(
        //     &format!(
        //         "{}{:.3}..{:.3} 1p {err_one_point:.6} mm {err_minmax:.6} r {err_refined:.6}",
        //         " ".repeat(rec.depth),
        //         rec.t0,
        //         rec.t1
        //     )
        //     .into(),
        // );

        if rec.depth < MAX_DEPTH && err > self.tolerance {
            let t = self.find_subdivision_point(rec);
            let utan_t = self.q.eval(t).to_vec2().normalize();
            // TODO(robustness): deal with derivative near-zero
            let cusp = self.cusp_sign(t);
            self.subdivide(rec, result, t, utan_t, cusp, cusp);
        } else {
            result.curve_to(c_approx.p1, c_approx.p2, c_approx.p3);
        }
    }

    fn subdivide(
        &self,
        rec: &OffsetRec,
        result: &mut BezPath,
        t: f64,
        utan_t: Vec2,
        cusp_t_minus: f64,
        cusp_t_plus: f64,
    ) {
        let rec0 = OffsetRec {
            t0: rec.t0,
            t1: t,
            utan0: rec.utan0,
            utan1: utan_t,
            cusp0: rec.cusp0,
            cusp1: cusp_t_minus,
            depth: rec.depth + 1,
        };
        self.offset_rec(&rec0, result);
        let rec1 = OffsetRec {
            t0: t,
            t1: rec.t1,
            utan0: utan_t,
            utan1: rec.utan1,
            cusp0: cusp_t_plus,
            cusp1: rec.cusp1,
            depth: rec.depth + 1,
        };
        self.offset_rec(&rec1, result);
    }

    fn apply(&self, rec: &OffsetRec, a: f64, b: f64) -> CubicBez {
        // Discussion question: should this clamp the derivatives to
        // point in the right direction?

        // wondering if p0 and p3 should be in rec
        // Scale factor from derivatives to displacements
        let s = (1. / 3.) * (rec.t1 - rec.t0);
        let p0 = self.c.eval(rec.t0) + self.d * turn(rec.utan0);
        let l0 = s * self.q.eval(rec.t0).to_vec2().length() + a * self.d;
        let mut p1 = p0;
        if l0 * rec.cusp0 > 0.0 {
            p1 += l0 * rec.utan0;
        }
        let p3 = self.c.eval(rec.t1) + self.d * turn(rec.utan1);
        let mut p2 = p3;
        let l1 = s * self.q.eval(rec.t1).to_vec2().length() - b * self.d;
        if l1 * rec.cusp1 > 0.0 {
            p2 -= l1 * rec.utan1;
        }
        //web_sys::console::log_1(&format!("l0={l0} l1={l1}").into());
        CubicBez::new(p0, p1, p2, p3)
    }

    fn utans(&self, rec: &OffsetRec) -> [Vec2; 3] {
        let q = self.q.subsegment(rec.t0..rec.t1);
        OFFSET_TS.map(|t| q.eval(t).to_vec2().normalize())
    }

    /// Approximate curve using one-point shape control.
    ///
    /// This solves (a, b) parameters (in terms of unit tangents) to place
    /// the t = 0.5 point on the approximation at the offset normal to
    /// t = 0.5 on the generatrix.
    fn one_point(&self, rec: &OffsetRec, utans: &[Vec2; 3]) -> (f64, f64) {
        let ca = rec.utan0;
        let cb = rec.utan1;
        let z = utans[1] - 0.5 * (ca + cb);
        let idet = (8.0 / 3.0) / ca.cross(cb);
        let a = -z.dot(cb) * idet;
        let b = z.dot(ca) * idet;
        (a, b)
    }

    /// Estimate the error of the curve approximation.
    ///
    /// The error estimate is good but not guaranteed to be conservative in
    /// all cases. Experiment suggests underestimate can happen when there
    /// is large curvature variation; perhaps a fudge factor based on that
    /// would help.
    ///
    /// Also do a Newton step to refine the `ts` values, to place them in
    /// the normal ray of the generatrix. The `ts` values are the parameter
    /// values on the approximation corresponding to `OFFSET_TS` on the
    /// generatrix.
    fn est_error(
        &self,
        rec: &OffsetRec,
        utans: &[Vec2; 3],
        c_approx: CubicBez,
        ts: &mut [f64; 3],
    ) -> f64 {
        let qa = c_approx.deriv();
        let mut max_err = 0.0;
        for i in 0..3 {
            let t = OFFSET_TS[i];
            let utan = utans[i];
            let t_orig = rec.t0 + t * (rec.t1 - rec.t0);
            let p = self.c.eval(t_orig) + self.d * turn(utan);
            let mut ta = ts[i];
            // Newton step to place c_approx(ta) in normal ray
            let pa = c_approx.eval(ta);
            let tana = qa.eval(ta).to_vec2();
            ta -= utan.dot(pa - p) / utan.dot(tana);
            let dist_err = p.distance(c_approx.eval(ta));
            let angle_err = qa.eval(ta).to_vec2().cross(utan);
            let err = dist_err + 0.15 * angle_err.abs();
            max_err = err.max(max_err);
            ts[i] = ta;
        }
        max_err
    }

    /// Approximate offset curve using linear minmax.
    fn linear_minmax(&self, rec: &OffsetRec, utans: &[Vec2; 3]) -> (f64, f64) {
        // Bezier weights; 2 and 3 are by symmetry
        const fn w0(t: f64) -> f64 {
            (1.0 - t) * (1.0 - t) * (1.0 - t)
        }
        const fn w1(t: f64) -> f64 {
            3.0 * (1.0 - t) * (1.0 - t) * t
        }
        const fn w01(t: f64) -> f64 {
            w0(t) + w1(t)
        }
        let t0 = T1_FOR_OFFSET;
        let t1 = 0.5;
        let t2 = 1.0 - t0;
        let ca = w1(t1) * rec.utan0.cross(utans[1]);
        let ca1 = ca + w1(t0) * rec.utan0.cross(utans[0]);
        let ca2 = ca + w1(t2) * rec.utan0.cross(utans[2]);
        let cb = w1(t1) * rec.utan1.cross(utans[1]);
        let cb1 = cb + w1(t2) * rec.utan1.cross(utans[0]);
        let cb2 = cb + w1(t0) * rec.utan1.cross(utans[2]);
        let cc = w01(t1) * (rec.utan0.dot(utans[1]) + rec.utan1.dot(utans[1])) - 2.0;
        let cc1 = cc + w01(t0) * rec.utan0.dot(utans[0]) + w01(t2) * rec.utan1.dot(utans[0]);
        let cc2 = cc + w01(t2) * rec.utan0.dot(utans[2]) + w01(t0) * rec.utan1.dot(utans[2]);
        let idet = 1.0 / (ca1 * cb2 - ca2 * cb1);
        let a = (cc1 * cb2 - cc2 * cb1) * idet;
        let b = (ca1 * cc2 - ca2 * cc1) * idet;
        (a, b)
    }

    /// Refine the minmax approximation, taking offset into account.
    ///
    /// Returns new (a, b) values
    fn refine_minmax(
        &self,
        rec: &OffsetRec,
        utans: &[Vec2; 3],
        c_approx: CubicBez,
        a: f64,
        b: f64,
        ts: [f64; 3],
    ) -> (f64, f64) {
        let mut ca = [0.0; 3];
        let mut cb = [0.0; 3];
        let mut cc = [0.0; 3];
        for i in 0..3 {
            let t = OFFSET_TS[i];
            let utan = utans[i];
            let t_orig = rec.t0 + t * (rec.t1 - rec.t0);
            let p = self.c.eval(t_orig);
            let n = turn(utan);
            let ta = ts[i];
            let mta = 1. - ta;
            let pa = c_approx.eval(ta);
            ca[i] = 3. * mta * ta * mta * rec.utan0.dot(n);
            cb[i] = 3. * mta * ta * ta * rec.utan1.dot(n);
            cc[i] = (pa - p).dot(n) - self.d;
        }
        let ca01 = ca[0] + ca[1];
        let cb01 = cb[0] + cb[1];
        let cc01 = cc[0] + cc[1];
        let ca12 = ca[1] + ca[2];
        let cb12 = cb[1] + cb[2];
        let cc12 = cc[1] + cc[2];
        let det = self.d * (ca01 * cb12 - ca12 * cb01);
        let da = (cc01 * cb12 - cc12 * cb01) / det;
        let db = (ca01 * cc12 - ca12 * cc01) / det;
        (a - da, b - db)
    }

    // Note: probably want to return unit tangent for robustness, rather
    // than re-computing it from the cubic.
    fn find_subdivision_point(&self, rec: &OffsetRec) -> f64 {
        // outline of work:
        // if segment contains inflection point, then t = 0.5 (can defer this, and might not be ideal)
        // else, average
        // alternative idea wrt inflection point: if solution count in t range is not 1, then
        // fall back to t = 0.5
        let mut t = 0.0;
        let mut n_soln = 0;
        // Note: do we want to track p0 & p3 in rec, to avoid repeated eval?
        let chord = self.c.eval(rec.t1) - self.c.eval(rec.t0);
        if chord.cross(rec.utan0) * chord.cross(rec.utan1) < 0.0 {
            let tan = rec.utan0 + rec.utan1;
            // set up quadratic equation for matching tangents
            let z0 = tan.cross(self.q.p0.to_vec2());
            let z1 = tan.cross(self.q.p1.to_vec2());
            let z2 = tan.cross(self.q.p2.to_vec2());
            let c0 = z0;
            let c1 = 2.0 * (z1 - z0);
            let c2 = (z2 - z1) - (z1 - z0);
            for root in solve_quadratic(c0, c1, c2) {
                if root > rec.t0 && root < rec.t1 {
                    t = root;
                    n_soln += 1;
                }
            }
        }
        if n_soln == 1 {
            //web_sys::console::log_1(&format!("{}..{} -> {t}", rec.t0, rec.t1).into());
            t
        } else {
            //web_sys::console::log_1(&format!("{}..{} -> midpoint", rec.t0, rec.t1).into());
            0.5 * (rec.t0 + rec.t1)
        }
    }
}

/// Rotate 90 degrees (clockwise in y-down).
///
/// a.dot(b) = a.cross(turn(b))
/// a.cross(y) = turn(a).dot(b)
fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}
