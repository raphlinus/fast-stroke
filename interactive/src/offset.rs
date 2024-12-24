//! An implementation of curve offset.

use kurbo::{
    common::solve_itp, BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, QuadBez, Vec2,
};

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
        CubicOffset {
            c,
            q,
            d,
            c0: d2 * p1xp0,
            c1: d2 * (p2xp0 - 2.0 * p1xp0),
            c2: d2 * (p2xp1 - p2xp0 + p1xp0),
            tolerance,
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
        let (a_minmax, b_minmax) = self.linear_minmax(rec);
        let c_minmax = self.apply(rec, a_minmax, b_minmax);
        ts = OFFSET_TS;
        let err_minmax = self.est_error(rec, &utans, c_minmax, &mut ts);
        if err_minmax < err {
            err = err_minmax;
            c_approx = c_minmax;
        }
        let (a_refined, b_refined) =
            self.refine_minmax(rec, &utans, c_approx, a_minmax, b_minmax, ts);
        let c_refined = self.apply(rec, a_refined, b_refined);
        let err_refined = self.est_error(rec, &utans, c_refined, &mut ts);
        if err_refined < err {
            err = err_refined;
            c_approx = c_refined;
        }

        web_sys::console::log_1(
            &format!(
                "{}{:.3}..{:.3} 1p {err_one_point:.6} mm {err_minmax:.6} r {err_refined:.6}",
                " ".repeat(rec.depth),
                rec.t0,
                rec.t1
            )
            .into(),
        );
        if rec.depth < MAX_DEPTH && err > self.tolerance {
            let t = rec.t0 + 0.5 * (rec.t1 - rec.t0);
            let utan_t = utans[1];
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
        let p1 = p0 + s * self.q.eval(rec.t0).to_vec2() + a * self.d * rec.utan0;
        let p3 = self.c.eval(rec.t1) + self.d * turn(rec.utan1);
        let p2 = p3 - s * self.q.eval(rec.t1).to_vec2() + b * self.d * rec.utan1;
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
    fn linear_minmax(&self, rec: &OffsetRec) -> (f64, f64) {
        // TODO: maybe take q as arg?
        let q = self.q.subsegment(rec.t0..rec.t1);
        let dt = rec.t1 - rec.t0;
        let b01 = dt * q.p0.to_vec2();
        let b12 = dt * q.p1.to_vec2();
        let b23 = dt * q.p2.to_vec2();
        let u0xb12 = rec.utan0.cross(b12);
        let u0xb23 = rec.utan0.cross(b23);
        let b01xu1 = b01.cross(rec.utan1);
        let b12xu1 = b12.cross(rec.utan1);
        let n0xb01 = rec.utan0.dot(b01);
        let n1xb01 = rec.utan1.dot(b01);
        let n0xb12 = rec.utan0.dot(b12);
        let n1xb12 = rec.utan1.dot(b12);
        let n0xb23 = rec.utan0.dot(b23);
        let n1xb23 = rec.utan1.dot(b23);
        // TODO: oppportunity to simplify algebra. But it's possible
        // we rework in favor of utans[].
        let coefs_at = |t: f64| {
            let mt = 1.0 - t;
            let ca = 6.0 * mt.powi(3) * t * t * u0xb12 + 3.0 * mt * mt * t.powi(3) * u0xb23;
            let cb = -3.0 * mt.powi(3) * t * t * b01xu1 + -6.0 * mt * mt * t.powi(3) * b12xu1;
            let cc = mt.powi(4) * (mt + 3.0 * t) * n0xb01
                + mt * mt * t * t * (3.0 * mt + t) * n1xb01
                + mt.powi(3) * t * (2.0 * mt + 6.0 * t) * n0xb12
                + mt * t.powi(3) * (6.0 * mt + 2.0 * t) * n1xb12
                + mt * mt * t * t * (mt + 3.0 * t) * n0xb23
                + t.powi(4) * (3.0 * mt + t) * n1xb23;
            (ca, cb, cc)
        };
        let [t0, t1, t2] = OFFSET_TS;
        let (ca0, cb0, cc0) = coefs_at(t0);
        let (ca1, cb1, cc1) = coefs_at(t1);
        let (ca2, cb2, cc2) = coefs_at(t2);
        let z0 = cc0 - dt * q.eval(t0).to_vec2().hypot();
        let z1 = cc1 - dt * q.eval(t1).to_vec2().hypot();
        let z2 = cc2 - dt * q.eval(t2).to_vec2().hypot();
        let ca01 = ca0 + ca1;
        let cb01 = cb0 + cb1;
        let z01 = z0 + z1;
        let ca12 = ca1 + ca2;
        let cb12 = cb1 + cb2;
        let z12 = z1 + z2;
        let det = ca01 * cb12 - ca12 * cb01;
        let a = (z01 * cb12 - z12 * cb01) / det;
        let b = (ca01 * z12 - ca12 * z01) / det;
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
}

/// Rotate 90 degrees (clockwise in y-down).
///
/// a.dot(b) = a.cross(turn(b))
/// a.cross(y) = turn(a).dot(b)
fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}
