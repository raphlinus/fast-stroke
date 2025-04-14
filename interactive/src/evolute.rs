// Copyright 2025 the Kurbo Authors
// SPDX-License-Identifier: Apache-2.0 OR MIT

use kurbo::{
    common::solve_itp, BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, QuadBez, Vec2,
};

const CUSP_EPSILON: f64 = 1e-12;

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
        let n0 = turn(self.q.eval(rec.t0).to_vec2().normalize()) * scale * rec.dr0.signum();
        let n1 = turn(self.q.eval(rec.t1).to_vec2().normalize()) * scale * rec.dr1.signum();
        let p1 = rec.p0 - n0;
        let p2 = rec.p1 + n1;
        path.curve_to(p1, p2, rec.p1);
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
        let rec0 = EvoluteRec {
            t0: rec.t0,
            t1: t,
            p0: rec.p0,
            p1: p,
            dr0: rec.dr0,
            dr1: dr_minus,
        };
        self.evolute_rec(&rec0, path);
        let rec1 = EvoluteRec {
            t0: t,
            t1: rec.t1,
            p0: p,
            p1: rec.p1,
            dr0: dr_plus,
            dr1: rec.dr1,
        };
        self.evolute_rec(&rec1, path);
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
    let rec = EvoluteRec {
        t0,
        t1,
        p0,
        p1,
        dr0,
        dr1,
    };
    ev.evolute_rec(&rec, &mut path);
    path
}

fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}
