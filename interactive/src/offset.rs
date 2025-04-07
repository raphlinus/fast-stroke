// Copyright 2022 the Kurbo Authors
// SPDX-License-Identifier: Apache-2.0 OR MIT

//! An implementation of curve offset.

use kurbo::{
    common::{solve_itp, solve_quadratic},
    BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, QuadBez, Shape, Vec2,
};

use crate::cusp::CuspAnalysis;

const N_LSE: usize = 8;

const BLEND: f64 = 1e-2;

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
    utans: [Vec2; N_LSE],
    p_offset: [Point; N_LSE],
}

// We never let cusp values haven an absolute value smaller than
// this. When a cusp is found, determine its sign and use this value.
const CUSP_EPSILON: f64 = 1e-12;

/// Maximum recursion depth
///
/// Perhaps should be configurable.
const MAX_DEPTH: usize = 8;

pub fn offset_cubic(c: CubicBez, d: f64, tolerance: f64) -> BezPath {
    let mut result = BezPath::new();
    let co = CubicOffset::new(c, d, tolerance);
    // TODO: cusp analysis to avoid divide by 0 in following math.
    let (cusp0, utan0) = co.cusp_and_utan(co.q.p0, co.c0);
    let (cusp1, utan1) = co.cusp_and_utan(co.q.p2, co.c0 + co.c1 + co.c2);
    result.move_to(c.p0 + d * turn(utan0));
    let rec = OffsetRec::new(&co, 0., 1., utan0, utan1, cusp0, cusp1, 0);
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
        let (a, b) = self.least_squares(rec);
        let mut ts = [0.0; N_LSE];
        for i in 0..N_LSE {
            ts[i] = (i + 1) as f64 * (1.0 / (N_LSE + 1) as f64);
        }
        let mut c_approx = self.apply(rec, a, b);
        let err_init = self.eval_err(rec, c_approx, &mut ts);
        let mut err = err_init;
        // speed/quality tradeoff: skip refinement if in tolerance
        let (mut a2, mut b2) = self.refine_least_squares(rec, c_approx, ts, a, b);
        for _ in 0..1 {
            self.eval_err(rec, c_approx, &mut ts);
            (a2, b2) = self.refine_least_squares(rec, c_approx, ts, a2, b2);
        }
        let c_approx2 = self.apply(rec, a2, b2);
        let err2 = self.eval_err(rec, c_approx2, &mut ts);
        if err2 < err {
            c_approx = c_approx2;
            err = err2;
        }
        let c_subseg = self.c.subsegment(rec.t0..rec.t1);
        web_sys::console::log_1(&c_subseg.into_path(0.0).to_svg().into());
        web_sys::console::log_1(
            &format!(
                "{}{:.3}..{:.3} init {err_init:.6} refined {err2:.6}",
                " ".repeat(rec.depth),
                rec.t0,
                rec.t1
            )
            .into(),
        );

        if rec.depth < MAX_DEPTH && err > self.tolerance.powi(2) {
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
        let rec0 = OffsetRec::new(
            self,
            rec.t0,
            t,
            rec.utan0,
            utan_t,
            rec.cusp0,
            cusp_t_minus,
            rec.depth + 1,
        );
        self.offset_rec(&rec0, result);
        let rec1 = OffsetRec::new(
            self,
            t,
            rec.t1,
            utan_t,
            rec.utan1,
            cusp_t_plus,
            rec.cusp1,
            rec.depth + 1,
        );
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

    // Compute least squares error approximation, returning (a, b)
    fn least_squares(&self, rec: &OffsetRec) -> (f64, f64) {
        let c = CubicBez::new(
            turn(rec.utan0).to_point(),
            turn(rec.utan0).to_point(),
            turn(rec.utan1).to_point(),
            turn(rec.utan1).to_point(),
        );
        let mut aa = 0.0;
        let mut ab = 0.0;
        let mut ac = 0.0;
        let mut bb = 0.0;
        let mut bc = 0.0;
        for i in 0..N_LSE {
            let t = (i + 1) as f64 * (1.0 / (N_LSE + 1) as f64);
            let n = turn(rec.utans[i]);
            let p = c.eval(t).to_vec2();
            let c_n = p.dot(n) - 1.0;
            let c_t = p.cross(n);
            let mt = 1.0 - t;
            let a_n = 3.0 * mt * t * mt * rec.utan0.dot(n);
            let a_t = 3.0 * mt * t * mt * rec.utan0.cross(n);
            let b_n = 3.0 * mt * t * t * rec.utan1.dot(n);
            let b_t = 3.0 * mt * t * t * rec.utan1.cross(n);
            aa += a_n * a_n + BLEND * a_t * a_t;
            ab += a_n * b_n + BLEND * a_t * b_t;
            ac += a_n * c_n + BLEND * a_t * c_t;
            bb += b_n * b_n + BLEND * b_t * b_t;
            bc += b_n * c_n + BLEND * b_t * c_t;
        }
        let idet = 1.0 / (aa * bb - ab * ab);
        let a_lse = -idet * (ac * bb - ab * bc);
        let b_lse = -idet * (aa * bc - ac * ab);

        // Cusp / very thick case
        let th = rec.utan1.cross(rec.utan0).atan2(rec.utan1.dot(rec.utan0));
        let arc_weight = self.d * th.abs();
        let lse_weight = self.c.eval(rec.t0).distance(self.c.eval(rec.t1));
        let blend = arc_weight / (arc_weight + lse_weight);
        let a_arc = (2. / 3.) / (1.0 + (0.5 * th).cos()) * 2.0 * (0.5 * th).sin();
        let b_arc = -a_arc;
        let a = a_arc * blend + a_lse * (1.0 - blend);
        let b = b_arc * blend + b_lse * (1.0 - blend);

        //let a = a_lse;
        //let b = b_lse;
        (a, b)
    }

    /// Evaluate error and also refine t values
    ///
    /// Returns squared absolute distance error.
    fn eval_err(&self, rec: &OffsetRec, c_approx: CubicBez, ts: &mut [f64; N_LSE]) -> f64 {
        let qa = c_approx.deriv();
        let mut max_err = 0.0;
        for i in 0..N_LSE {
            let mut ta = ts[i];
            // TODO: probably should also store in rec, we always need this
            let utan = rec.utans[i];
            let p = rec.p_offset[i];
            // Newton step to refine ta value
            let pa = c_approx.eval(ta);
            let tana = qa.eval(ta).to_vec2();
            ta -= utan.dot(pa - p) / utan.dot(tana);
            ts[i] = ta;
            let dist_err_squared = p.distance_squared(c_approx.eval(ta));
            // Note: would be very cheap to also include angle error
            // let angle_err = qa.eval(ta).to_vec2().cross(utan);
            // 0.15 is based off n=3, should decrease
            // let err = dist_err + 0.15 * angle_err.abs();
            max_err = dist_err_squared.max(max_err);
        }
        max_err
    }

    fn refine_least_squares(
        &self,
        rec: &OffsetRec,
        c_approx: CubicBez,
        ts: [f64; N_LSE],
        a: f64,
        b: f64,
    ) -> (f64, f64) {
        let mut aa = 0.0;
        let mut ab = 0.0;
        let mut ac = 0.0;
        let mut bb = 0.0;
        let mut bc = 0.0;
        for i in 0..N_LSE {
            let n = turn(rec.utans[i]);
            let t = ts[i];
            let err_vec = c_approx.eval(t) - rec.p_offset[i];
            let c_n = err_vec.dot(n);
            let c_t = err_vec.cross(n);
            let mt = 1.0 - t;
            let a_n = 3.0 * mt * t * mt * rec.utan0.dot(n);
            let a_t = 3.0 * mt * t * mt * rec.utan0.cross(n);
            let b_n = 3.0 * mt * t * t * rec.utan1.dot(n);
            let b_t = 3.0 * mt * t * t * rec.utan1.cross(n);
            aa += a_n * a_n + BLEND * a_t * a_t;
            ab += a_n * b_n + BLEND * a_t * b_t;
            ac += a_n * c_n + BLEND * a_t * c_t;
            bb += b_n * b_n + BLEND * b_t * b_t;
            bc += b_n * c_n + BLEND * b_t * c_t;
        }
        let idet = 1.0 / (self.d * (aa * bb - ab * ab));
        let delta_a = idet * (ac * bb - ab * bc);
        let delta_b = idet * (aa * bc - ac * ab);
        (a - delta_a, b - delta_b)
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

impl OffsetRec {
    fn new(
        co: &CubicOffset,
        t0: f64,
        t1: f64,
        utan0: Vec2,
        utan1: Vec2,
        cusp0: f64,
        cusp1: f64,
        depth: usize,
    ) -> Self {
        let mut utans = [Vec2::ZERO; N_LSE];
        let mut p_offset = [Point::ZERO; N_LSE];
        let dt = (t1 - t0) * (1.0 / (N_LSE + 1) as f64);
        for i in 0..N_LSE {
            let t = t0 + (i + 1) as f64 * dt;
            // TODO: deal with zero derivative
            let utan = co.q.eval(t).to_vec2().normalize();
            utans[i] = utan;
            p_offset[i] = co.c.eval(t) + co.d * turn(utan);
        }
        OffsetRec {
            t0,
            t1,
            utan0,
            utan1,
            cusp0,
            cusp1,
            depth,
            utans,
            p_offset,
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
