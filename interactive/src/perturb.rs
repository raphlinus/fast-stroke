use kurbo::QuadBez;
// Doing lots of experiments, will try things and then move on.
#[allow(unused)]
use kurbo::{
    common::solve_cubic, Affine, BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, Vec2,
};

fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}

const ERROR_SCALE: f64 = 20.0;

pub struct CurveOffset {
    c: CubicBez,
    q: QuadBez,
    n0: Vec2,
    n1: Vec2,
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

const T1_FOR_OFFSET: f64 = 0.2;
const OFFSET_TS: [f64; 3] = [T1_FOR_OFFSET, 0.5, 1. - T1_FOR_OFFSET];

pub struct OffsetSolution {
    // offset
    d: f64,
    // parameters
    a: f64,
    b: f64,
    // t values on approximation corresponding to T1_FOR_OFFSET, 1/2 and 1 - T1_FOR_OFFSET
    ts: [f64; 3],
}

impl CurveOffset {
    pub fn new(c: CubicBez) -> Self {
        let q = c.deriv();
        let b01 = q.p0.to_vec2();
        let b12 = q.p1.to_vec2();
        let b23 = q.p2.to_vec2();
        let n0 = turn(b01).normalize();
        let n1 = turn(b23).normalize();
        let b01xb12 = b01.cross(b12);
        let b01xb23 = b01.cross(b23);
        let b12xb23 = b12.cross(b23);
        let n0xb01 = n0.cross(b01);
        let n1xb01 = n1.cross(b01);
        let n0xb12 = n0.cross(b12);
        let n1xb12 = n1.cross(b12);
        let n0xb23 = n0.cross(b23);
        let n1xb23 = n1.cross(b23);
        CurveOffset {
            c,
            q,
            n0,
            n1,
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

    fn coefs_at(&self, t: f64) -> (f64, f64, f64) {
        let mt = 1.0 - t;
        let ca = 6.0 * mt.powi(3) * t * t * self.b01xb12 + 3.0 * mt * mt * t.powi(3) * self.b01xb23;
        let cb =
            -3.0 * mt.powi(3) * t * t * self.b01xb23 + -6.0 * mt * mt * t.powi(3) * self.b12xb23;
        let cc = mt.powi(4) * (mt + 3.0 * t) * self.n0xb01
            + mt * mt * t * t * (3.0 * mt + t) * self.n1xb01
            + mt.powi(3) * t * (2.0 * mt + 6.0 * t) * self.n0xb12
            + mt * t.powi(3) * (6.0 * mt + 2.0 * t) * self.n1xb12
            + mt * mt * t * t * (mt + 3.0 * t) * self.n0xb23
            + t.powi(4) * (3.0 * mt + t) * self.n1xb23;
        (ca, cb, cc)
    }

    fn error_at(&self, t: f64, a: f64, b: f64) -> f64 {
        let (ca, cb, cc) = self.coefs_at(t);
        -(ca * a + cb * b + cc) / self.q.eval(t).to_vec2().hypot() - 1.0
    }

    fn deriv_error_at(&self, t: f64, a: f64, b: f64) -> f64 {
        // Sure we could do this analytically, but just experimenting
        const DT: f64 = 1e-6;
        (self.error_at(t + DT, a, b) - self.error_at(t - DT, a, b)) * (0.5 / DT)
    }

    pub fn make_delta(&self, a: f64, b: f64) -> CubicBez {
        let b01 = self.q.p0.to_vec2();
        let b23 = self.q.p2.to_vec2();
        CubicBez::new(
            self.n0.to_point(),
            (self.n0 + a * b01).to_point(),
            (self.n1 + b * b23).to_point(),
            self.n1.to_point(),
        )
    }

    pub fn apply(&self, a: f64, b: f64, d: f64) -> CubicBez {
        let delta = self.make_delta(a, b);
        CubicBez::new(
            self.c.p0 + d * delta.p0.to_vec2(),
            self.c.p1 + d * delta.p1.to_vec2(),
            self.c.p2 + d * delta.p2.to_vec2(),
            self.c.p3 + d * delta.p3.to_vec2(),
        )
    }

    /// Do a Newton step to refine a t value on the approximation.
    ///
    /// Given curve parameters, a t value on the generatrix, and an
    /// approximate t value on the approximation, refine the latter.
    fn newton_step_t(&self, a: f64, b: f64, d: f64, t: f64, ta: f64) -> f64 {
        let ca = self.apply(a, b, d);
        let p = self.c.eval(t);
        let tan = self.q.eval(t).to_vec2();
        let pa = ca.eval(ta);
        let error = tan.dot(pa - p);
        let qa = ca.deriv();
        ta - error / tan.dot(qa.eval(ta).to_vec2())
    }
}

// Produce the (a, b) for the given curve, making error zero at
// t = 1/3 and 2/3 in the linear approximation.
pub fn two_point_approx(c: CubicBez) -> (f64, f64) {
    let co = CurveOffset::new(c);
    let t0 = 1.0 / 3.0;
    let t1 = 2.0 / 3.0;
    let (ca0, cb0, cc0) = co.coefs_at(t0);
    let (ca1, cb1, cc1) = co.coefs_at(t1);
    let z0 = -co.q.eval(t0).to_vec2().hypot() - cc0;
    let z1 = -co.q.eval(t1).to_vec2().hypot() - cc1;
    let det = ca0 * cb1 - ca1 * cb0;
    let a = (z0 * cb1 - z1 * cb0) / det;
    let b = (ca0 * z1 - ca1 * z0) / det;
    (a, b)
}

const PLOT_TRANSLATE: Vec2 = Vec2::new(100., 200.);
const PLOT_SCALE: Affine = Affine::scale(500.0);

pub fn plot_error(c: CubicBez, a: f64, b: f64) -> BezPath {
    let mut result = BezPath::new();
    let co = CurveOffset::new(c);
    const N: usize = 50;
    for i in 0..=N {
        let t = i as f64 / N as f64;
        let y = ERROR_SCALE * co.error_at(t, a, b);
        let p = PLOT_SCALE * Point::new(t, y) + PLOT_TRANSLATE;
        if i == 0 {
            result.move_to(p);
        } else {
            result.line_to(p);
        }
    }
    result
}

/// Approximate error using cubic spline Hermite interpolation.
///
/// This works well for smooth errors, not so much when velocity varies.
pub fn spline_error(c: CubicBez, a: f64, b: f64) -> BezPath {
    let mut result = BezPath::new();
    let co = CurveOffset::new(c);
    result.move_to(Point::ZERO);
    let mut hermite = |t0, y0, d0, t1, y1, d1| {
        let d_scale = (1. / 3.) * (t1 - t0);
        let p1 = (t0 + d_scale, ERROR_SCALE * (y0 + d_scale * d0));
        let p2 = (t1 - d_scale, ERROR_SCALE * (y1 - d_scale * d1));
        let p3 = (t1, ERROR_SCALE * y1);
        result.curve_to(p1, p2, p3);
    };
    let t1 = 1. / 6.;
    let y1 = co.error_at(t1, a, b);
    let d1 = co.deriv_error_at(t1, a, b);
    hermite(0.0, 0.0, 0.0, t1, y1, d1);
    let t2 = 0.5;
    let y2 = co.error_at(t2, a, b);
    let d2 = co.deriv_error_at(t2, a, b);
    hermite(t1, y1, d1, t2, y2, d2);
    let t3 = 1.0 - t1;
    let y3 = co.error_at(t3, a, b);
    let d3 = co.deriv_error_at(t3, a, b);
    hermite(t2, y2, d2, t3, y3, d3);
    hermite(t3, y3, d3, 1.0, 0.0, 0.0);
    Affine::translate(PLOT_TRANSLATE) * PLOT_SCALE * result
}

pub fn est_err_bounds(c: CubicBez, a: f64, b: f64) -> [f64; 3] {
    let co = CurveOffset::new(c);
    let mut max_err = 0.0;
    let mut max_d_err = 0.0;
    let t1 = 1. / 6.;
    [t1, 0.5, 1.0 - t1].map(|t| {
        let y = co.error_at(t, a, b);
        let d = co.deriv_error_at(t, a, b);
        max_err = y.abs();
        max_d_err = d.abs();
        let d_scale = 0.3;
        ERROR_SCALE * 500.0 * max_err.hypot(d_scale * max_d_err)
    })
}

pub fn plot(xys: &[(f64, f64)]) -> BezPath {
    let mut result = BezPath::new();
    for (x, y) in xys {
        let y = ERROR_SCALE * y;
        let p = Affine::translate((100., 200.)) * Affine::scale(500.) * Point::new(*x, y);
        if result.elements().is_empty() {
            result.move_to(p);
        } else {
            result.line_to(p);
        }
    }
    result
}

pub fn error_by_rays(c: CubicBez, offset: f64, approx: CubicBez) -> Vec<(f64, f64)> {
    let q = c.deriv();
    const N: usize = 50;
    (1..N)
        .map(|i| {
            let t = i as f64 / N as f64;
            let p = c.eval(t);
            let d = q.eval(t).to_vec2().normalize();
            let dnorm = d.normalize();
            let p1 = 3.0 * (approx.p1 - approx.p0);
            let p2 =
                3.0 * approx.p2.to_vec2() - 6.0 * approx.p1.to_vec2() + 3.0 * approx.p0.to_vec2();
            let p3 = (approx.p3 - approx.p0) - 3.0 * (approx.p2 - approx.p1);
            let c0 = (approx.p0 - p).dot(d);
            let c1 = p1.dot(d);
            let c2 = p2.dot(d);
            let c3 = p3.dot(d);
            let mut best = 10f64;
            for t in solve_cubic(c0, c1, c2, c3) {
                if (0.0..=1.0).contains(&t) {
                    let approx_p = approx.eval(t);
                    let err = dnorm.cross(approx_p - p) - offset;
                    //web_sys::console::log_1(&format!("{:?}", err).into());
                    if err.abs() < best.abs() {
                        best = err;
                    }
                }
            }
            (t, best / offset)
        })
        .collect()
}

// This is a hacky thing to get quick empirical estimates of scaling exponents.
// It should probably become a systematic thing.
pub fn scaling_test(c: CubicBez, d: f64) {
    let mut errs = vec![];
    for i in 0..5 {
        let scale = 0.5f64.powf(i as f64 / 5.0);
        let c_s = c.subsegment(0.5 - 0.5 * scale..0.5 + 0.5 * scale);
        let (a, b) = two_point_approx(c_s);
        let co = CurveOffset::new(c_s);
        let c_offset = co.apply(a, b, d);
        let mut max_err = 0.0;
        for (_, e) in error_by_rays(c_s, d, c_offset) {
            max_err = e.max(max_err);
        }
        errs.push(max_err);
    }
    web_sys::console::log_1(&format!("{errs:?}").into());
}

/// Linear approximation to error minmax
///
/// The basic idea is to make e(1/6) + e(1/2) and e(1/2) + e(5/6) both
/// zero. This roughly minimizes the maximum error for the first Chebyshev
/// polynomial that can't be zeroed exactly.
pub fn linear_minmax(c: CubicBez) -> (f64, f64) {
    let co = CurveOffset::new(c);
    // probably want to renumber to 1, 2, 3, as 0 is not start point
    let t0 = 1. / 6.;
    let t1 = 0.5;
    let t2 = 1.0 - t0;
    let (ca0, cb0, cc0) = co.coefs_at(t0);
    let (ca1, cb1, cc1) = co.coefs_at(t1);
    let (ca2, cb2, cc2) = co.coefs_at(t2);
    let z0 = -co.q.eval(t0).to_vec2().hypot() - cc0;
    let z1 = -co.q.eval(t1).to_vec2().hypot() - cc1;
    let z2 = -co.q.eval(t2).to_vec2().hypot() - cc2;
    let ca01 = ca0 + ca1;
    let cb01 = cb0 + cb1;
    let z01 = z0 + z1;
    let ca12 = ca1 + ca2;
    let cb12 = cb1 + cb2;
    let z12 = z1 + z2;
    let det = ca01 * cb12 - ca12 * cb01;
    let dot = ca01 * cb01 + ca12 * cb12;
    web_sys::console::log_1(&format!("{:.6} {:.6}", det * 1e-9, dot * 1e-9).into());
    let a = (z01 * cb12 - z12 * cb01) / det;
    let b = (ca01 * z12 - ca12 * z01) / det;
    (a, b)
}

/// One point shape control, straight up
pub fn one_point(c: CubicBez) -> (f64, f64) {
    let co = CurveOffset::new(c);
    let b01 = co.q.p0.to_vec2();
    let b23 = co.q.p2.to_vec2();
    let t1 = 0.5;
    let (w0, w1, w2, w3) = (0.125, 0.375, 0.375, 0.125);
    let n1 = turn(co.q.eval(t1).to_vec2().normalize());
    let ca = w1 * b01;
    let cb = w2 * b23;
    let cc = (w0 + w1) * co.n0 + (w2 + w3) * co.n1;
    let z = n1 - cc;
    let det = ca.cross(cb);
    let a = z.cross(cb) / det;
    let b = ca.cross(z) / det;
    (a, b)
}

/// One point shape control, adjustable ta
pub fn one_point_at(c: CubicBez, d: f64, ta: f64) -> (f64, f64) {
    let co = CurveOffset::new(c);
    let b01 = co.q.p0.to_vec2();
    let b23 = co.q.p2.to_vec2();
    let n1 = turn(co.q.eval(0.5).to_vec2().normalize());
    let dp = c.eval(0.5) + d * n1 - c.eval(ta);
    let mt = 1.0 - ta;
    let w0 = mt * mt * mt;
    let w1 = 3.0 * mt * mt * ta;
    let w2 = 3.0 * mt * ta * ta;
    let w3 = ta * ta * ta;
    let ca = w1 * b01;
    let cb = w2 * b23;
    let cc = (w0 + w1) * co.n0 + (w2 + w3) * co.n1;
    let z = dp - d * cc;
    let det = d * ca.cross(cb);
    let a = z.cross(cb) / det;
    let b = ca.cross(z) / det;
    (a, b)
}

// Given a t value, do a Newton step to minimize angle error
// at generatrix t = 0.5
pub fn refine_one_point(c: CubicBez, d: f64, t: f64) -> f64 {
    let co = CurveOffset::new(c);
    const DT: f64 = 1e-6;
    let errs = [-DT, DT].map(|dt| {
        let (a, b) = one_point_at(c, d, t + dt);
        let approx = co.apply(a, b, d);
        let tana = approx.deriv().eval(t + dt).to_vec2();
        let tan = co.q.eval(t).to_vec2();
        let angle_err = tana.cross(tan).abs() / tan.hypot();
        angle_err
    });
    let new_t = t - (errs[0] + errs[1]) / (errs[1] - errs[0]) * DT;
    web_sys::console::log_1(&format!("{errs:?} {new_t}").into());
    new_t
}

pub fn brute_one_point(c: CubicBez, d: f64, t: f64) -> f64 {
    let co = CurveOffset::new(c);
    let mut best_err = 1e9;
    let mut best_t = 0.5;
    const DT: f64 = 1e-4;
    for i in -1000..=1000 {
        let dt = DT * i as f64;
        let (a, b) = one_point_at(c, d, t + dt);
        let mut soln = OffsetSolution::from_a_b(a, b, d);
        let [y0, y1, y2] = soln.refine_ts(&co);
        let y = y0.max(y1).max(y2);
        if y < best_err {
            best_err = y;
            best_t = t + dt;
        }
    }
    web_sys::console::log_1(&format!("{best_t}").into());
    best_t
}

pub fn linear_center(c: CubicBez) -> (f64, f64) {
    let co = CurveOffset::new(c);
    // probably want to renumber to 1, 2, 3, as 0 is not start point
    let t0 = 1. / 6.;
    let t1 = 0.5;
    let t2 = 1.0 - t0;
    let (ca0, cb0, cc0) = co.coefs_at(t0);
    let (ca1, cb1, cc1) = co.coefs_at(t1);
    let (ca2, cb2, cc2) = co.coefs_at(t2);
    let z0 = -co.q.eval(t0).to_vec2().hypot() - cc0;
    let z1 = -co.q.eval(t1).to_vec2().hypot() - cc1;
    let z2 = -co.q.eval(t2).to_vec2().hypot() - cc2;
    let denom = ca0 + 2.0 * ca1 + ca2 - cb0 - 2. * cb1 - cb2;
    let a = (z0 + 2.0 * z1 + z2) / denom;
    let b = -a;
    (a, b)
}

/// Logs refinement of t on normal ray. Just for exploring the convergence of that.
pub fn refine(c: CubicBez, a: f64, b: f64, d: f64) {
    let co = CurveOffset::new(c);
    let t = 1. / 6.;
    let mut ta = t;
    for i in 0..5 {
        ta = co.newton_step_t(a, b, d, t, ta);
        web_sys::console::log_1(&format!("{i}: {:.6}", ta).into());
    }
}

// Estimate error by refining t approximation, Euclidean distance to
// approximation, and tangent discrepancy.
//
// This is very good but not a rigorous bound that could be enforced by
// property testing. I would consider that future work.
pub fn est_err_refined(c: CubicBez, a: f64, b: f64, d: f64) -> [f64; 3] {
    let co = CurveOffset::new(c);
    let ca = co.apply(a, b, d);
    let t1 = 0.2;
    [t1, 0.5, 1.0 - t1].map(|t| {
        let tan = co.q.eval(t).to_vec2();
        let p = c.eval(t) + d * turn(tan.normalize());
        let mut ta = t;
        for _ in 0..1 {
            ta = co.newton_step_t(a, b, d, t, ta);
        }
        let pa = ca.eval(ta);
        let dist_err = p.distance(pa);
        let tana = ca.deriv().eval(ta).to_vec2();
        let angle_err = tana.cross(tan).abs() / tan.hypot();
        let err = dist_err + 0.15 * angle_err;
        ERROR_SCALE * 500.0 * err / d
    })
}

impl OffsetSolution {
    pub fn from_a_b(a: f64, b: f64, d: f64) -> Self {
        let ts = OFFSET_TS;
        Self { d, a, b, ts }
    }

    pub fn apply(&self, co: &CurveOffset) -> CubicBez {
        co.apply(self.a, self.b, self.d)
    }

    /// Refine t values, using a Newton step to place them closer to the
    /// normal ray from the generatrix.
    ///
    /// Return an error estimate. This is 3 values for the purpose of
    /// visualization, but for production we probably only need 1. (That
    /// said, it's not out of the question these can guide subdivision,
    /// though my current thinking is to use the velocity denominator)
    pub fn refine_ts(&mut self, co: &CurveOffset) -> [f64; 3] {
        let ca = self.apply(co);
        let mut errors = [0.0; 3];
        for i in 0..3 {
            let t = OFFSET_TS[i];
            let tan = co.q.eval(t).to_vec2();
            let p = co.c.eval(t) + self.d * turn(tan.normalize());
            let ta = co.newton_step_t(self.a, self.b, self.d, t, self.ts[i]);
            let pa = ca.eval(ta);
            let dist_err = p.distance(pa);
            let tana = ca.deriv().eval(ta).to_vec2();
            let angle_err = tana.cross(tan).abs() / tan.hypot();
            let err = dist_err + 0.15 * angle_err;
            errors[i] = ERROR_SCALE * 500.0 * err / self.d;
            self.ts[i] = ta;
        }
        errors
    }

    /// Refine a and b values, using a Newton step towards minmax.
    pub fn refine_minmax(&mut self, co: &CurveOffset) {
        let b01 = co.q.p0.to_vec2();
        let b23 = co.q.p2.to_vec2();
        let mut ca = [0.0; 3];
        let mut cb = [0.0; 3];
        let mut cc = [0.0; 3];
        let c_approx = self.apply(co);
        for i in 0..3 {
            let t = OFFSET_TS[i];
            let p = co.c.eval(t);
            let tan = co.q.eval(t).to_vec2();
            let n = turn(tan).normalize();
            let ta = self.ts[i];
            let mta = 1. - ta;
            let pa = c_approx.eval(ta);
            ca[i] = 3. * mta * mta * ta * b01.dot(n);
            cb[i] = 3. * mta * ta * ta * b23.dot(n);
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
        self.a -= da;
        self.b -= db;
    }
}

#[test]
fn foo() {
    use kurbo::Vec2;
    let x = Vec2::new(1.0, 0.0);
    let y = Vec2::new(0.0, 1.0);
    println!("{}", x.cross(y));
}
