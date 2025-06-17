// Copyright 2025 the Kurbo Authors
// SPDX-License-Identifier: Apache-2.0 OR MIT

use kurbo::{
    common::solve_itp, BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, QuadBez, Vec2,
};

const CUSP_EPSILON: f64 = 1e-12;

const N_LSE: usize = 8;

const BLEND: f64 = 1e-2;

struct CubicEvolute {
    c: CubicBez,
    q: QuadBez,
    c0: f64,
    c1: f64,
    c2: f64,
}

impl CubicEvolute {
    fn new(c: CubicBez) -> Self {
        let q = c.deriv();
        let p1xp0 = q.p1.to_vec2().cross(q.p0.to_vec2());
        let p2xp0 = q.p2.to_vec2().cross(q.p0.to_vec2());
        let p2xp1 = q.p2.to_vec2().cross(q.p1.to_vec2());
        Self {
            c,
            q,
            c0: 2.0 * p1xp0,
            c1: 2.0 * p2xp0 - 4.0 * p1xp0,
            c2: 2.0 * p2xp1 - 2.0 * p2xp0 + 2.0 * p1xp0,
        }
    }

    #[allow(unused)]
    fn curvature(&self, t: f64) -> f64 {
        let num = self.c0 + t * self.c1 + t * t * self.c2;
        let denom = self.q.eval(t).to_vec2().length_squared().powf(1.5);
        num / denom
    }

    // reciprocal of curvature
    fn radius(&self, t: f64) -> f64 {
        let num = self.c0 + t * self.c1 + t * t * self.c2;
        let denom = self.q.eval(t).to_vec2().length_squared().powf(1.5);
        denom / num
    }

    fn radius_deriv(&self, t: f64) -> f64 {
        // TODO: consider analytical method
        const DELTA: f64 = 1e-6;
        (0.5 / DELTA) * (self.radius(t + DELTA) - self.radius(t - DELTA))
    }

    fn eval(&self, t: f64) -> Point {
        let n = turn(self.q.eval(t).to_vec2().normalize());
        self.c.eval(t) - n * self.radius(t)
    }

    fn apply(&self, rec: &EvoluteRec, a: f64, b: f64) -> CubicBez {
        let p1 = rec.p0 + a * rec.utan0;
        let p2 = rec.p1 + b * rec.utan1;
        CubicBez::new(rec.p0, p1, p2, rec.p1)
    }

    fn evolute_rec(&self, rec: &EvoluteRec, path: &mut BezPath) {
        if rec.dr0 * rec.dr1 < 0.0 {
            let a = rec.t0;
            let b = rec.t1;
            let s = rec.dr1.signum();
            let f = |t| s * self.radius_deriv(t);
            let k1 = 0.2 / (b - a);
            const ITP_EPS: f64 = 1e-12;
            let t = solve_itp(f, a, b, ITP_EPS, 1, k1, s * rec.dr0, s * rec.dr1);
            let p = self.eval(t);
            let dr_minus = CUSP_EPSILON.copysign(rec.dr0);
            let dr_plus = CUSP_EPSILON.copysign(rec.dr1);
            self.subdivide(rec, path, t, p, dr_minus, dr_plus);
            return;
        }
        let chord = rec.p0.distance(rec.p1);
        let scale = (1. / 3.) * chord;
        let (mut a, mut b) = (-scale, scale);
        let dt = (rec.t1 - rec.t0) * (1.0 / (N_LSE + 1) as f64);
        let mut ts: [f64; N_LSE] = core::array::from_fn(|i| rec.t0 + (i + 1) as f64 * dt);
        //web_sys::console::log_1(&format!("pre: {ts:?}").into());
        let mut ca = self.apply(rec, a, b);
        let mut max_err = self.eval_err(rec, ca, &mut ts);
        for i in 0..5 {
            let (a2, b2) = self.refine_least_squares(rec, ca, &ts, a, b);
            let ca2 = self.apply(rec, a2, b2);
            let err = self.eval_err(rec, ca2, &mut ts);
            web_sys::console::log_1(&format!("{i}: {err}?{max_err}").into());
            if err < max_err {
                (a, b) = (a2, b2);
                ca = ca2;
                max_err = err;
            } else {
                break;
            }
        }
        path.curve_to(ca.p1, ca.p2, ca.p3);
    }

    fn subdivide(
        &self,
        rec: &EvoluteRec,
        path: &mut BezPath,
        t: f64,
        p: Point,
        dr_minus: f64,
        dr_plus: f64,
    ) {
        let utan = turn(self.q.eval(t).to_vec2()).normalize();
        let rec0 = EvoluteRec {
            t0: rec.t0,
            t1: t,
            p0: rec.p0,
            p1: p,
            utan0: rec.utan0,
            utan1: utan * dr_minus.signum(),
            dr0: rec.dr0,
            dr1: dr_minus,
        };
        self.evolute_rec(&rec0, path);
        let rec1 = EvoluteRec {
            t0: t,
            t1: rec.t1,
            p0: p,
            p1: rec.p1,
            utan0: utan * dr_plus.signum(),
            utan1: rec.utan1,
            dr0: dr_plus,
            dr1: rec.dr1,
        };
        self.evolute_rec(&rec1, path);
    }

    /// Evaluate error and refine t values.
    ///
    /// Returns squared absolute difference error.
    ///
    /// Adjust t values on the source curve.
    fn eval_err(&self, _rec: &EvoluteRec, c_approx: CubicBez, ts: &mut [f64; N_LSE]) -> f64 {
        let qa = c_approx.deriv();
        let mut max_err = 0.0;
        for i in 0..N_LSE {
            let t_approx = (i + 1) as f64 * (1.0 / (N_LSE + 1) as f64);
            let t_orig = ts[i];
            let p = self.eval(t_orig);
            let pa = c_approx.eval(t_approx);
            //web_sys::console::log_1(&format!("p={p:?}@{t_orig} pa={pa:?}@{t_approx}").into());
            let tana = qa.eval(t_approx).to_vec2();
            const DT: f64 = 1e-6;
            // Numerical differentiation of evolute; might be a good
            // idea to replace with analytical
            let p_plus = self.eval(t_orig + DT);
            let dp_dt = (p_plus - p) * (1.0 / DT);
            let error = tana.dot(pa - p);
            let t_orig_new = t_orig + error / tana.dot(dp_dt);
            ts[i] = t_orig_new;
            let p_new = self.eval(t_orig_new);
            let dist_err_squared = p_new.distance_squared(pa);
            max_err = dist_err_squared.max(max_err);
        }
        max_err
    }

    fn refine_least_squares(
        &self,
        rec: &EvoluteRec,
        c_approx: CubicBez,
        ts: &[f64; N_LSE],
        a: f64,
        b: f64,
    ) -> (f64, f64) {
        let mut aa = 0.0;
        let mut ab = 0.0;
        let mut ac = 0.0;
        let mut bb = 0.0;
        let mut bc = 0.0;
        for i in 0..N_LSE {
            let t = (i + 1) as f64 * (1.0 / (N_LSE + 1) as f64);
            // n is tangent to source curve, so normal to evolute
            let t_orig = ts[i];
            let n = self.q.eval(t_orig).to_vec2().normalize();
            let p_orig = self.eval(t_orig);
            let err_vec = c_approx.eval(t) - p_orig;
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
        let idet = 1.0 / (aa * bb - ab * ab);
        let delta_a = idet * (ac * bb - ab * bc);
        let delta_b = idet * (aa * bc - ac * ab);
        (a - delta_a, b - delta_b)
    }
}

// Approximate an evolute by sampling.
pub fn evolute_hacky_approx(c: CubicBez) -> BezPath {
    let ev = CubicEvolute::new(c);
    let mut path = BezPath::new();
    const N: usize = 50;
    for i in 0..=N {
        let t = i as f64 / N as f64;
        let p = ev.eval(t);
        if i == 0 {
            path.move_to(p);
        } else {
            path.line_to(p);
        }
    }
    path
}

// Approximate with a single cubic Bezier
//
// This is only going to work reasonably well in the monotonic curvature case.
// Likely we should skip this and go to subdivision. Subdivide at curvature extrema
// in the source curve (these are cusps in the evolute).
pub fn evolute_approx_cubic(c: CubicBez) -> CubicBez {
    let ev = CubicEvolute::new(c);
    let p0 = ev.eval(0.0);
    let p3 = ev.eval(1.0);
    let chord = p0.distance(p3);
    let scale = (1. / 3.) * chord;
    let n0 = turn(ev.q.p0.to_vec2().normalize()) * scale * ev.radius_deriv(0.0).signum();
    let p1 = p0 - n0;
    let n1 = turn(ev.q.p2.to_vec2().normalize()) * scale * ev.radius_deriv(1.0).signum();
    let p2 = p3 + n1;
    CubicBez::new(p0, p1, p2, p3)
}

struct EvoluteRec {
    t0: f64,
    t1: f64,
    p0: Point,
    p1: Point,
    utan0: Vec2,
    utan1: Vec2,
    // derivative of radius. This is in the rec because it vanishes at subdivision points
    dr0: f64,
    dr1: f64,
}

pub fn evolute_approx(c: CubicBez) -> BezPath {
    let ev = CubicEvolute::new(c);
    let mut path = BezPath::new();
    let t0 = 0.0;
    let t1 = 1.0;
    let p0 = ev.eval(t0);
    let p1 = ev.eval(t1);
    path.move_to(p0);
    let dr0 = ev.radius_deriv(t0);
    let dr1 = ev.radius_deriv(t1);
    let q = c.deriv();
    let utan0 = turn(q.p0.to_vec2()).normalize() * dr0.signum();
    let utan1 = turn(q.p2.to_vec2()).normalize() * dr1.signum();
    let rec = EvoluteRec {
        t0,
        t1,
        p0,
        p1,
        utan0,
        utan1,
        dr0,
        dr1,
    };
    ev.evolute_rec(&rec, &mut path);
    path
}

fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}
