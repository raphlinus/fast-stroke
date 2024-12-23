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

struct OffsetMath {
    // This is the derivative of the subsegment. We might
    // get rid of it in favor of sampling the original `q`.
    q: QuadBez,
    // Unit tangents for `OFFSET_TS`
    utans: [Vec2; 3],
    b01xb12: f64,
    b01xb23: f64,
    b12xb23: f64,
    n0xb01: f64,
    n1xb01: f64,
    n0xb12: f64,
    n1xb12: f64,
    n0xb23: f64,
    n1xb23: f64,
}

// We never let cusp values haven an absolute value smaller than
// this. When a cusp is found, determine its sign and use this value.
const CUSP_EPSILON: f64 = 1e-12;

/// Maximum recursion depth
///
/// Perhaps should be configurable.
const MAX_DEPTH: usize = 4;

// t values for minmax and error estimation
const T1_FOR_OFFSET: f64 = 0.2;
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
        let d0 = q.p0.to_vec2();
        let d1 = 2.0 * (q.p1 - q.p0);
        let d2 = (q.p2 - q.p1) - (q.p1 - q.p0);
        CubicOffset {
            c,
            q,
            d,
            c0: d * d1.cross(d0),
            c1: d * 2.0 * d2.cross(d0),
            c2: d * d2.cross(d1),
            tolerance,
        }
        // possible optimizations:
        // c0 = 2 p1 x p0
        // c1 = 2 (p2 x p0 - 2 p1 x p0)
        // c2 = 2 (p2 - p1) x (p1 - p0) = 2 (p2 x p1 - p2 x p0 + p1 x p0)
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
        let rsqrt = 1.0 / tan.to_vec2().hypot2().sqrt();
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
        let math = OffsetMath::new(self, rec);
        let (a, b) = math.one_point(rec);
        let mut ts = OFFSET_TS;
        let c_approx = self.apply(rec, a, b);
        let err = math.est_error(self, rec, c_approx, &mut ts);
        web_sys::console::log_1(
            &format!("{}..{} ab {a:.6} {b:.6}, err {err}", rec.t0, rec.t1).into(),
        );
        if rec.depth < MAX_DEPTH && err > self.tolerance {
            let t = rec.t0 + 0.5 * (rec.t1 - rec.t0);
            let utan_t = math.utans[1];
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
}

impl OffsetMath {
    fn new(co: &CubicOffset, rec: &OffsetRec) -> Self {
        let q = co.q.subsegment(rec.t0..rec.t1);
        let utans = OFFSET_TS.map(|t| q.eval(t).to_vec2().normalize());
        let b01 = q.p0.to_vec2();
        let b12 = q.p1.to_vec2();
        let b23 = q.p2.to_vec2();
        // TODO: these should be reformulated in terms
        // of unit tangent vectors.
        let b01xb12 = b01.cross(b12);
        let b01xb23 = b01.cross(b23);
        let b12xb23 = b12.cross(b23);
        let n0 = turn(rec.utan0);
        let n1 = turn(rec.utan1);
        // These are actually dot products, which might be
        // faster with SIMD. TODO: get sign right.
        let n0xb01 = n0.cross(b01);
        let n1xb01 = n1.cross(b01);
        let n0xb12 = n0.cross(b12);
        let n1xb12 = n1.cross(b12);
        let n0xb23 = n0.cross(b23);
        let n1xb23 = n1.cross(b23);
        Self {
            q,
            utans,
            b01xb12,
            b01xb23,
            b12xb23,
            n0xb01,
            n1xb01,
            n0xb12,
            n1xb12,
            n0xb23,
            n1xb23,
        }
    }

    // Return (a, b) parameters in terms of unit tangents
    fn one_point(&self, rec: &OffsetRec) -> (f64, f64) {
        let ca = rec.utan0;
        let cb = rec.utan1;
        let z = self.utans[1] - 0.5 * (ca + cb);
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
        co: &CubicOffset,
        rec: &OffsetRec,
        c_approx: CubicBez,
        ts: &mut [f64; 3],
    ) -> f64 {
        let qa = c_approx.deriv();
        let mut max_err = 0.0;
        for i in 0..3 {
            let t = OFFSET_TS[i];
            let utan = self.utans[i];
            let t_orig = rec.t0 + t * (rec.t1 - rec.t0);
            let p = co.c.eval(t_orig) + co.d * turn(utan);
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
}

/// Rotate 90 degrees (clockwise in y-down).
///
/// a.dot(b) = a.cross(turn(b))
/// a.cross(y) = turn(a).dot(b)
fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}
