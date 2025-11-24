use kurbo::{ParamCurveCurvature, QuadBez};
// Doing lots of experiments, will try things and then move on.
#[allow(unused)]
use kurbo::{
    common::solve_cubic, Affine, BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, Vec2,
};

pub fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}

const ERROR_SCALE: f64 = 1.0e2;
const BLEND: f64 = 1e-3;

pub struct CurveOffset {
    c: CubicBez,
    pub(crate) q: QuadBez,
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

pub struct OffsetSolutionLse {
    // offset
    d: f64,
    // parameters
    a: f64,
    b: f64,
    // t values on approximation corresponding to T1_FOR_OFFSET, 1/2 and 1 - T1_FOR_OFFSET
    ts: [f64; LSE_N],
    // t values on generatrix
    ts_rev: [f64; LSE_N],
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

    /// Do a Newton step to refine a t value on the generatrix.
    ///
    /// Given curve parameters, a t value on the generatrix, and an
    /// approximate t value on the approximation, refine the former.
    fn newton_step_rev_t(&self, a: f64, b: f64, d: f64, t: f64, ta: f64) -> f64 {
        let ca = self.apply(a, b, d);
        let p = self.c.eval(t);
        let qa = ca.deriv();
        let tan = qa.eval(ta).to_vec2();
        let pa = ca.eval(ta);
        let error = tan.dot(pa - p);
        t + error / tan.dot(self.q.eval(t).to_vec2())
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
    let dot = ca0 * cb0 + ca1 * cb1;
    let ratio = det / dot;
    web_sys::console::log_1(&format!("two_point ratio = {ratio:.6}").into());
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
    plot_scaled(xys, ERROR_SCALE)
}

pub fn plot_scaled(xys: &[(f64, f64)], scale: f64) -> BezPath {
    let mut result = BezPath::new();
    for (x, y) in xys {
        let p = Affine::translate((100., 200.)) * Affine::scale(500.) * Point::new(*x, scale * y);
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
            let mut best = 1e3f64;
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
    //let dot = ca01 * cb01 + ca12 * cb12;
    //web_sys::console::log_1(&format!("{:.6} {:.6}", det * 1e-9, dot * 1e-9).into());
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

/// One point shape control, adjustable t on source curve
///
/// Note: this function originally adjusted the t on the approximation,
/// keeping t = 0.5 on source curve fixed.
pub fn one_point_at(c: CubicBez, d: f64, t: f64) -> (f64, f64) {
    let co = CurveOffset::new(c);
    let b01 = co.q.p0.to_vec2();
    let b23 = co.q.p2.to_vec2();
    let n1 = turn(co.q.eval(t).to_vec2().normalize());
    let dp = c.eval(t) + d * n1 - c.eval(0.5);
    let (w0, w1, w2, w3) = (0.125, 0.375, 0.375, 0.125);
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
        let tana = approx.deriv().eval(0.5).to_vec2();
        let tan = co.q.eval(t + dt).to_vec2().normalize();
        let angle_err = tana.cross(tan);
        angle_err
    });
    let new_t = t - (errs[0] + errs[1]) / (errs[1] - errs[0]) * DT;
    web_sys::console::log_1(&format!("refine one point {errs:?} {new_t}").into());
    new_t
}

pub fn compute_angle_err(co: &CurveOffset, d: f64, t: f64) -> f64 {
    let (a, b) = one_point_at(co.c, d, t);
    let approx = co.apply(a, b, d);
    let tana = approx.deriv().eval(0.5).to_vec2();
    let tan = co.q.eval(t).to_vec2().normalize();
    tana.cross(tan)
}

/// Compute the derivative of angle error wrt t
///
/// Work in progress: compute analytically
pub fn angle_err_deriv(co: &CurveOffset, d: f64, t: f64) -> f64 {
    // first compute numerically
    //const DT: f64 = 1e-6;
    //let errs = [-DT, DT].map(|dt| compute_angle_err(&co, d, t + dt));
    //let num_deriv = (errs[1] - errs[0]) / (2.0 * DT);

    let k = co.c.curvature(t);

    // Note: these should probably be computed from unit tangents
    let b01 = co.q.p0.to_vec2();
    let b23 = co.q.p2.to_vec2();
    // Note: this det is missing a factor of d, also scaled 3/8 compared with other code
    //let det = (3. / 8.) * b01.cross(b23);
    let idet = (8. / 3.) / b01.cross(b23);

    let tan = co.q.eval(t).to_vec2();
    let utan = tan.normalize();

    let dp = co.c.eval(t) + d * turn(utan) - co.c.eval(0.5);
    let cc = co.n0 + co.n1;
    let z = dp - d * 0.5 * cc;
    // Note: these (a, b) have a factor of d, compared with one_point_at
    let a = z.cross(b23) * idet;
    let b = b01.cross(z) * idet;

    // tana_base is invariant to t
    let tana_base = co.q.eval(0.5).to_vec2() + 1.5 * d * (co.n1 - co.n0);
    let tana = tana_base + 0.75 * (b * b23 - a * b01);

    // Multiply these by 3/8(1 + k d) / det to get true da/dt and db/dt
    let dadt = tan.cross(b23);
    let dbdt = b01.cross(tan);

    //web_sys::console::log_1(&format!("d ab {dadt} {dbdt} {da_dt_a} {db_dt_a}").into());

    let dutan_dt = -k * turn(tan);

    /*
    let tanp = co.q.eval(t + DT).to_vec2().normalize();
    let tanm = co.q.eval(t - DT).to_vec2().normalize();
    let dtan_dt_numeric = (tanp - tanm) / (2.0 * DT);
    web_sys::console::log_1(&format!("num = {dtan_dt_numeric} deriv = {dtan_dt}").into());
    */

    let factor = 0.75 * (1.0 + k * d) * idet;
    let dtana_dt = factor * (dbdt * b23 - dadt * b01);

    /*
    let approxp = co.apply(a + DT * dadt, b + DT * dbdt, d);
    let approxm = co.apply(a - DT * dadt, b - DT * dbdt, d);
    let tanap = approxp.deriv().eval(0.5).to_vec2();
    let tanam = approxm.deriv().eval(0.5).to_vec2();
    let dtana_dt = (tanap - tanam) / (2.0 * DT);
    */

    let deriv = tana.cross(dutan_dt) + dtana_dt.cross(utan);

    //web_sys::console::log_1(&format!("num = {num_deriv} deriv = {deriv}").into());

    deriv
}

const TRY_ITP: bool = true;

/// Solve for t such that t=0.5 on approximation matches position and tangent.
pub fn solve_midpoint(co: &CurveOffset, d: f64) -> Option<f64> {
    if let Some(t) = try_midpoint_newton(co, d, 0.5) {
        return Some(t);
    }
    let roots = midpoint_cubic_approx(co, d);
    for root in &roots {
        if let Some(t) = try_midpoint_newton(co, d, *root) {
            return Some(t);
        }
    }
    try_midpoint_itp(co, d)
}

fn try_midpoint_newton(co: &CurveOffset, d: f64, t0: f64) -> Option<f64> {
    const THRESH: f64 = 1e-12; // TODO: make configurable
    const MAX_ITERS: usize = 10;
    let mut t = t0;
    let mid_err = compute_angle_err(co, d, t);
    web_sys::console::log_1(&format!("init err = {mid_err}").into());
    let mut err = mid_err;
    for _ in 0..MAX_ITERS {
        if err.abs() < THRESH {
            let (a, b) = one_point_at(co.c, d, t);
            if a * d <= -1. / 3. || b * d >= 1. / 3. {
                break;
            }
            return Some(t);
        }
        //let errs = [-DT, DT].map(|dt| compute_angle_err(&co, d, t + dt));
        //t -= (errs[0] + errs[1]) / (errs[1] - errs[0]) * DT;
        t -= err / angle_err_deriv(co, d, t);
        if t < 0.0 || t > 1.0 {
            break;
        }
        let new_err = compute_angle_err(co, d, t);
        //web_sys::console::log_1(&format!("t = {t}, err = {new_err}").into());
        if new_err.abs() > err.abs() {
            break;
        }
        err = new_err;
    }
    //web_sys::console::log_1(&format!("newton failure t0={t0}").into());
    None
}

fn try_midpoint_itp(co: &CurveOffset, d: f64) -> Option<f64> {
    const THRESH: f64 = 1e-12; // TODO: make configurable
                               // First, find brackets containing solution
    const INIT_DELTA: f64 = 1. / 256.;
    let mut delta = INIT_DELTA;
    let mut t = 0.5;
    let mut t_lo = t;
    let mut t_hi = t;
    let mid_err = compute_angle_err(co, d, t);
    let mut err = mid_err;
    let mut err_lo = mid_err;
    let mut err_hi = mid_err;
    let mid_sign = mid_err.signum();
    loop {
        t = 0.5 + delta;
        let (a, b) = one_point_at(co.c, d, t);
        web_sys::console::log_1(&format!("t {t} a {} b {}", a * d, b * d).into());
        if a * d > -1. / 3. && b * d < 1. / 3. {
            err = compute_angle_err(&co, d, t);
            if err * mid_sign < 0.0 {
                break;
            }
            if delta > 0. {
                err_hi = err;
                t_hi = t;
            } else {
                err_lo = err;
                t_lo = t;
            }
        } else {
        }
        if delta > 0. {
            delta = -delta;
        } else {
            delta = -2.0 * delta;
            if delta == 1.0 {
                // failure
                return None;
            }
        }
    }
    let (t0, t1, e0, e1, sign) = if delta > 0.0 {
        (t_hi, t, err_hi, err, -mid_sign)
    } else {
        (t, t_lo, err, err_lo, mid_sign)
    };
    let f = |t| sign * compute_angle_err(&co, d, t);
    let ya = sign * e0;
    let yb = sign * e1;
    let k1 = 0.2 / (t1 - t0);
    let soln = kurbo::common::solve_itp(f, t0, t1, THRESH, 1, k1, ya, yb);
    web_sys::console::log_1(&format!("itp {t0}..{t1} ya = {ya} yb = {yb} -> {soln}").into());
    Some(soln)
}

/// Find cubic approximation of midpoint G1 function and solve.
///
/// Results are filtered to have positive control point distances
fn midpoint_cubic_approx(co: &CurveOffset, d: f64) -> Vec<f64> {
    const T0: f64 = 0.5;
    const DT: f64 = 1e-3;
    let e0 = compute_angle_err(co, d, T0);
    let e1 = compute_angle_err(co, d, T0 + DT);
    let d0 = DT * angle_err_deriv(co, d, T0);
    let d1 = DT * angle_err_deriv(co, d, T0 + DT);
    let c0 = e0;
    let c1 = d0;
    let c2 = 3. * (e1 - e0) - 2. * d0 - d1;
    let c3 = -2. * (e1 - e0) + d0 + d1;
    let roots = solve_cubic(c0, c1, c2, c3);
    let mut filtered = roots.iter().map(|x| T0 + x * DT).collect::<Vec<_>>();
    filtered.sort_by_key(|x| (x - (T0 + 0.5 * DT)).abs().to_bits());
    web_sys::console::log_1(&format!("{filtered:?}").into());
    filtered
}

/// Brute-force search for root of angle error.
///
/// Selects the second one (useful for middle of cubic-like soln) but this
/// is not fully robust.
///
/// Useful for experimentation to decide if this middle soln is good.
pub fn brute_midpoint(c: CubicBez, d: f64) -> f64 {
    let co = CurveOffset::new(c);
    let mut last = 0.0;
    let mut crossings = 0;
    let mut t_best = 0.0;
    const N: usize = 200;
    for i in 0..=N {
        let t = i as f64 / N as f64;
        let err = compute_angle_err(&co, d, t);
        if i > 0 && err * last < 0.0 {
            t_best = t - (0.5 / N as f64);
            crossings += 1;
            if crossings == 2 {
                break;
            }
        }
        last = err;
    }
    t_best
}

pub fn onept_err_plot(c: CubicBez, d: f64) -> Vec<(f64, f64)> {
    let co = CurveOffset::new(c);
    const N: usize = 50;
    (0..=N)
        .map(|i| {
            let t = i as f64 / N as f64;
            let (a, b) = one_point_at(c, d, t);
            let approx = co.apply(a, b, d);
            let tana = approx.deriv().eval(0.5).to_vec2();
            let tan = co.q.eval(t).to_vec2().normalize();
            let angle_err = tana.cross(tan);
            (t, 2e-4 * angle_err)
        })
        .collect()
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

/// Number of subdivisions for least squares integration
/// We do dumb integration, but maybe something like Gauss-Lagrange would be better.
const LSE_N: usize = 8;

/// Compute solution based on minimizing least-squares distance error
pub fn least_squares(c: CubicBez) -> (f64, f64) {
    let co = CurveOffset::new(c);
    let mut aa_n = 0.0;
    let mut aa_t = 0.0;
    let mut ab_n = 0.0;
    let mut ab_t = 0.0;
    let mut ac_n = 0.0;
    let mut ac_t = 0.0;
    let mut bb_n = 0.0;
    let mut bb_t = 0.0;
    let mut bc_n = 0.0;
    let mut bc_t = 0.0;
    let mut cc_n = 0.0;
    let mut cc_t = 0.0;
    let utan0 = co.q.p0.to_vec2().normalize();
    let utan1 = co.q.p2.to_vec2().normalize();
    let c = CubicBez::new(
        turn(utan0).to_point(),
        turn(utan0).to_point(),
        turn(utan1).to_point(),
        turn(utan1).to_point(),
    );
    for i in 0..LSE_N {
        let t = (i + 1) as f64 / (LSE_N + 1) as f64;
        let utan = co.q.eval(t).to_vec2().normalize();
        let n = turn(utan);
        let c_n = c.eval(t).to_vec2().dot(n) - 1.0;
        let c_t = c.eval(t).to_vec2().cross(n);
        let mt = 1.0 - t;
        let a_n = 3.0 * mt * t * mt * utan0.dot(n);
        let a_t = 3.0 * mt * t * mt * utan0.cross(n);
        let b_n = 3.0 * mt * t * t * utan1.dot(n);
        let b_t = 3.0 * mt * t * t * utan1.cross(n);
        aa_n += a_n * a_n;
        aa_t += a_t * a_t;
        ab_n += a_n * b_n;
        ab_t += a_t * b_t;
        ac_n += a_n * c_n;
        ac_t += a_t * c_t;
        bb_n += b_n * b_n;
        bb_t += b_t * b_t;
        bc_n += b_n * c_n;
        bc_t += b_t * c_t;
        cc_n += c_n * c_n;
        cc_t += c_t * c_t;
    }
    let blend = BLEND;
    let aa = aa_n + blend * aa_t;
    let ab = ab_n + blend * ab_t;
    let ac = ac_n + blend * ac_t;
    let bb = bb_n + blend * bb_t;
    let bc = bc_n + blend * bc_t;
    let idet = 1.0 / (aa * bb - ab * ab);
    // web_sys::console::log_1(
    //     &format!("det = {:.3e}, aa * bb = {:.3e}", aa * bb - ab * ab, aa * bb).into(),
    // );
    let a = -idet * (ac * bb - ab * bc);
    let b = -idet * (aa * bc - ac * ab);
    // Don't do the scaling by length when moving to offset.rs
    let a_scaled = a / co.q.p0.to_vec2().length();
    let b_scaled = b / co.q.p2.to_vec2().length();
    let err_min =
        a * a * aa_n + 2.0 * a * b * ab_n + 2.0 * a * ac_n + b * b * bb_n + 2.0 * b * bc_n + cc_n;
    let err_max = a * a * (aa_n + aa_t)
        + 2.0 * a * b * (ab_n + ab_t)
        + 2.0 * a * (ac_n + ac_t)
        + b * b * (bb_n + bb_t)
        + 2.0 * b * (bc_n + bc_t)
        + (cc_n + cc_t);
    // web_sys::console::log_1(&format!("0: err_min = {err_min:.8}, err_max = {err_max:.8}").into());
    (a_scaled, b_scaled)
}

/// Compute solution based on simply drawing an arc approximation
pub fn draw_arc(c: CubicBez) -> (f64, f64) {
    let co = CurveOffset::new(c);
    let utan0 = co.q.p0.to_vec2().normalize();
    let utan1 = co.q.p2.to_vec2().normalize();
    let co = CurveOffset::new(c);
    let th = utan1.cross(utan0).atan2(utan1.dot(utan0));
    let a_arc = (2. / 3.) / (1.0 + (0.5 * th).cos()) * 2.0 * (0.5 * th).sin();
    let b_arc = -a_arc;
    let a_scaled = a_arc / co.q.p0.to_vec2().length();
    let b_scaled = b_arc / co.q.p2.to_vec2().length();
    (a_scaled, b_scaled)
}

/// Refine arc drawing solution to have zero linear error at midpoint
pub fn arc_onept_linear(c: CubicBez) -> (f64, f64) {
    let co = CurveOffset::new(c);
    let utan0 = co.q.p0.to_vec2().normalize();
    let utan1 = co.q.p2.to_vec2().normalize();
    let th = utan1.cross(utan0).atan2(utan1.dot(utan0));
    let a_arc = (2. / 3.) / (1.0 + (0.5 * th).cos()) * 2.0 * (0.5 * th).sin();
    let b_arc = -a_arc;
    let mid = 0.5 * turn(utan0 + utan1) + (3. / 8.) * (a_arc * utan0 + b_arc * utan1);
    let nm = turn(co.q.eval(0.5).to_vec2().normalize());
    let delta = 1.0 - mid.dot(nm);
    let dot0 = utan0.dot(nm);
    let dot1 = utan1.dot(nm);
    let ratio = (8. / 3.) * delta / (dot0 * dot0 + dot1 * dot1);
    let a_refined = a_arc + dot0 * ratio;
    let b_refined = b_arc + dot1 * ratio;
    // web_sys::console::log_1(&format!("mid = {mid:.4} nm = {nm:.4} delta = {delta:.4}").into());
    // let mid_refined = 0.5 * turn(utan0 + utan1) + (3. / 8.) * (a_refined * utan0 + b_refined * utan1);
    // let dot_refined = mid_refined.dot(nm);
    // web_sys::console::log_1(&format!("dot_refined = {dot_refined:.4}").into());

    let a_scaled = a_refined / co.q.p0.to_vec2().length();
    let b_scaled = b_refined / co.q.p2.to_vec2().length();
    (a_scaled, b_scaled)
}

/// Variant of arc onepoint with second order error approximation
pub fn arc_onept_quadratic(c: CubicBez, d: f64) -> (f64, f64) {
    let co = CurveOffset::new(c);
    let utan0 = co.q.p0.to_vec2().normalize();
    let utan1 = co.q.p2.to_vec2().normalize();
    let th = utan1.cross(utan0).atan2(utan1.dot(utan0));
    let a_arc = (2. / 3.) / (1.0 + (0.5 * th).cos()) * 2.0 * (0.5 * th).sin();
    let b_arc = -a_arc;

    let r = 1.0 / c.curvature(0.5);
    let r_offset = r + d;
    // center of osculating circle
    let center = c.eval(0.5) - turn(co.q.eval(0.5).to_vec2().normalize()) * r;
    web_sys::console::log_1(&format!("r = {r:.4}, center {center:.4}").into());
    let a_scaled = a_arc / co.q.p0.to_vec2().length();
    let b_scaled = b_arc / co.q.p2.to_vec2().length();
    // Note: we could compute this much more directly, but this is clear & concise
    let p_approx = co.apply(a_scaled, b_scaled, d).eval(0.5);
    let goal_vec = (p_approx - center).normalize();
    let goal = center + r_offset.abs() * goal_vec;
    web_sys::console::log_1(&format!("p_approx = {p_approx:.4}, goal = {goal:.4}").into());
    let delta = (goal - p_approx).dot(goal_vec) / d;
    web_sys::console::log_1(&format!("delta = {delta:.4}").into());
    let dot0 = utan0.dot(goal_vec);
    let dot1 = utan1.dot(goal_vec);
    let ratio = (8. / 3.) * delta / (dot0 * dot0 + dot1 * dot1);

    let a_refined = a_arc + dot0 * ratio;
    let b_refined = b_arc + dot1 * ratio;

    let ar_scaled = a_refined / co.q.p0.to_vec2().length();
    let br_scaled = b_refined / co.q.p2.to_vec2().length();
    (ar_scaled, br_scaled)
}

/// New version of one-point that places t = 0.5 on true curve
///
/// Return value is (a, b, t) where t is point on true curve.
pub fn arc_onept(c: CubicBez, d: f64) -> (f64, f64, f64) {
    let co = CurveOffset::new(c);
    let utan0 = co.q.p0.to_vec2().normalize();
    let utan1 = co.q.p2.to_vec2().normalize();
    let th = utan1.cross(utan0).atan2(utan1.dot(utan0));
    let a_arc = (2. / 3.) / (1.0 + (0.5 * th).cos()) * 2.0 * (0.5 * th).sin();
    let b_arc = -a_arc;

    let mid = 0.5 * turn(utan0 + utan1) + (3. / 8.) * (a_arc * utan0 + b_arc * utan1);
    let p_approx = c.eval(0.5) + d * mid;

    let nm = turn(co.q.eval(0.5).to_vec2().normalize());
    let mut n = nm;
    let dot0 = utan0.dot(nm);
    let dot1 = utan1.dot(nm);
    let dir = utan0 * dot0 + utan1 * dot1;
    // now solve for t so that true offset(t) x dir = nm x dir
    let mut t = 0.5;
    // Point on true offset
    let mut p = co.c.eval(t) + d * nm;
    // TODO: detect convergence and break
    // hypothesis: one Newton step is enough to establish scaling
    for i in 0..2 {
        // Newton step
        let dpdt = co.q.eval(t).to_vec2() * (1.0 + d * c.curvature(t));
        let err = (p - p_approx).cross(dir);
        web_sys::console::log_1(&format!("{i}: err = {err}").into());
        t -= err / dpdt.cross(dir);
        n = turn(co.q.eval(t).to_vec2().normalize());
        p = co.c.eval(t) + d * n;
    }
    let delta = (p - p_approx).dot(nm);
    let ratio = (8. / 3.) * delta / (d * (dot0 * dot0 + dot1 * dot1));

    let a_refined = a_arc + dot0 * ratio;
    let b_refined = b_arc + dot1 * ratio;

    let ar_scaled = a_refined / co.q.p0.to_vec2().length();
    let br_scaled = b_refined / co.q.p2.to_vec2().length();

    let c_soln = co.apply(ar_scaled, br_scaled, d);
    let p2 = c_soln.eval(0.5);
    web_sys::console::log_1(&format!("true {p:.4?} approx {p2:.4?}").into());
    (ar_scaled, br_scaled, t)
}

/// Two point along with estimate of transverse component of error.
///
/// Multiply third tuple value by d^2 to get estimate.
pub fn two_pt_w_err(c: CubicBez) -> (f64, f64, f64) {
    let co = CurveOffset::new(c);
    let mut aa = 0.0;
    let mut ab = 0.0;
    let mut ac = 0.0;
    let mut bb = 0.0;
    let mut bc = 0.0;
    let utan0 = co.q.p0.to_vec2().normalize();
    let utan1 = co.q.p2.to_vec2().normalize();
    let c = CubicBez::new(
        turn(utan0).to_point(),
        turn(utan0).to_point(),
        turn(utan1).to_point(),
        turn(utan1).to_point(),
    );
    for t in [1. / 3., 2. / 3.] {
        let utan = co.q.eval(t).to_vec2().normalize();
        let n = turn(utan);
        let c_n = c.eval(t).to_vec2().dot(n) - 1.0;
        let mt = 1.0 - t;
        let a_n = 3.0 * mt * t * mt * utan0.dot(n);
        let b_n = 3.0 * mt * t * t * utan1.dot(n);
        aa += a_n * a_n;
        ab += a_n * b_n;
        ac += a_n * c_n;
        bb += b_n * b_n;
        bc += b_n * c_n;
    }
    let idet = 1.0 / (aa * bb - ab * ab);
    let a = -idet * (ac * bb - ab * bc);
    let b = -idet * (aa * bc - ac * ab);
    let mut max_err = 0.0;
    // Compute transverse error
    for t in [1. / 3., 2. / 3.] {
        let mt = 1.0 - t;
        let utan = co.q.eval(t).to_vec2().normalize();
        let n = turn(utan);
        let c_t = c.eval(t).to_vec2().cross(n);
        let a_t = 3.0 * mt * t * mt * utan0.cross(n);
        let b_t = 3.0 * mt * t * t * utan1.cross(n);
        let transverse = a_t * a + b_t * b + c_t;
        let k = co.c.curvature(t);
        let err = 0.5 * k.abs() * transverse * transverse;
        web_sys::console::log_1(&format!("transverse {transverse:.4} err {err:.4}").into());
        max_err = err.max(max_err)
    }
    // Don't do the scaling by length when moving to offset.rs
    let a_scaled = a / co.q.p0.to_vec2().length();
    let b_scaled = b / co.q.p2.to_vec2().length();
    (a_scaled, b_scaled, max_err)
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

impl OffsetSolutionLse {
    pub fn from_a_b(a: f64, b: f64, d: f64) -> Self {
        let mut ts = [0.0; LSE_N];
        for i in 0..LSE_N {
            ts[i] = (i + 1) as f64 / (LSE_N + 1) as f64;
        }
        let ts_rev = ts;
        Self {
            d,
            a,
            b,
            ts,
            ts_rev,
        }
    }

    pub fn apply(&self, co: &CurveOffset) -> CubicBez {
        co.apply(self.a, self.b, self.d)
    }

    pub fn refine_lse(&mut self, co: &CurveOffset) {
        let mut aa_n = 0.0;
        let mut aa_t = 0.0;
        let mut ab_n = 0.0;
        let mut ab_t = 0.0;
        let mut ac_n = 0.0;
        let mut ac_t = 0.0;
        let mut bb_n = 0.0;
        let mut bb_t = 0.0;
        let mut bc_n = 0.0;
        let mut bc_t = 0.0;
        let mut cc_n = 0.0;
        let mut cc_t = 0.0;
        let utan0 = co.q.p0.to_vec2().normalize();
        let utan1 = co.q.p2.to_vec2().normalize();
        let c_approx = self.apply(co);
        let q_approx = c_approx.deriv();
        for i in 0..LSE_N {
            let t_approx = (i + 1) as f64 / (LSE_N + 1) as f64;
            // Newton step is done in eval_err
            //let t = co.newton_step_t(self.a, self.b, self.d, t_gen, self.ts[i]);
            //self.ts[i] = t;
            let t = self.ts[i];
            // May need to flip sign inside cusp
            let utan = q_approx.eval(t_approx).to_vec2().normalize();
            let p = co.c.eval(t);
            let n = turn(utan);
            let p_off = p + self.d * n;
            let err_vec = c_approx.eval(t_approx) - p_off;
            let c_n = err_vec.dot(n);
            let c_t = err_vec.cross(n);
            let mt = 1.0 - t;
            let a_n = 3.0 * mt * t * mt * utan0.dot(n);
            let a_t = 3.0 * mt * t * mt * utan0.cross(n);
            let b_n = 3.0 * mt * t * t * utan1.dot(n);
            let b_t = 3.0 * mt * t * t * utan1.cross(n);
            aa_n += a_n * a_n;
            aa_t += a_t * a_t;
            ab_n += a_n * b_n;
            ab_t += a_t * b_t;
            ac_n += a_n * c_n;
            ac_t += a_t * c_t;
            bb_n += b_n * b_n;
            bb_t += b_t * b_t;
            bc_n += b_n * c_n;
            bc_t += b_t * c_t;
            cc_n += c_n * c_n;
            cc_t += c_t * c_t;
        }
        let blend = BLEND;
        let aa = aa_n + blend * aa_t;
        let ab = ab_n + blend * ab_t;
        let ac = ac_n + blend * ac_t;
        let bb = bb_n + blend * bb_t;
        let bc = bc_n + blend * bc_t;
        let idet = 1.0 / (aa * bb - ab * ab);
        // Don't do the scaling by length when moving to offset.rs
        let a = -idet * (ac * bb - ab * bc);
        let b = -idet * (aa * bc - ac * ab);
        let a_scaled = a / (self.d * co.q.p0.to_vec2().length());
        let b_scaled = b / (self.d * co.q.p2.to_vec2().length());
        self.a += a_scaled;
        self.b += b_scaled;
        let err_min = a * a * aa_n
            + 2.0 * a * b * ab_n
            + 2.0 * a * ac_n
            + b * b * bb_n
            + 2.0 * b * bc_n
            + cc_n;
        let err_min = err_min / self.d.powi(2);
        let err_max = a * a * (aa_n + aa_t)
            + 2.0 * a * b * (ab_n + ab_t)
            + 2.0 * a * (ac_n + ac_t)
            + b * b * (bb_n + bb_t)
            + 2.0 * b * (bc_n + bc_t)
            + (cc_n + cc_t);
        let err_max = err_max / self.d.powi(2);
        // web_sys::console::log_1(
        //     &format!("   err_min = {err_min:.8}, err_max = {err_max:.8}").into(),
        // );
    }

    /// Evaluate error and also do Newton step on `ts`.
    pub fn eval_err(&mut self, co: &CurveOffset) -> f64 {
        let c_approx = self.apply(co);
        let mut err = 0.0;
        for i in 0..LSE_N {
            let t_gen = (i + 1) as f64 / (LSE_N + 1) as f64;
            let t = co.newton_step_t(self.a, self.b, self.d, t_gen, self.ts[i]);
            self.ts[i] = t;
            let utan = co.q.eval(t_gen).to_vec2().normalize();
            let p = co.c.eval(t_gen);
            let n = turn(utan);
            let p_off = p + self.d * n;
            let err_vec = c_approx.eval(t) - p_off;
            err += err_vec.length_squared();
        }
        err /= self.d.powi(2);
        // web_sys::console::log_1(&format!("err = {err:.8}").into());
        err
    }

    // Do the Newton step to refine generatrix t's
    pub fn newton_step_rev(&mut self, co: &CurveOffset) {
        for i in 0..LSE_N {
            let ta = (i + 1) as f64 / (LSE_N + 1) as f64;
            let t = co.newton_step_rev_t(self.a, self.b, self.d, self.ts_rev[i], ta);
            self.ts_rev[i] = t;
        }
    }

    // Osculating parabola error measurement
    pub fn est_arc_error_dot(&self, co: &CurveOffset) -> f64 {
        let delta = co.make_delta(self.a, self.b);
        const N: usize = 8;
        let mut err = 0.0;
        for i in 0..N {
            let t = (i as f64 + 0.5) * (1. / N as f64);
            let unorm = turn(co.q.eval(t).to_vec2().normalize());
            let approx = delta.eval(t).to_vec2();
            let dot = unorm.dot(approx);
            let cross = unorm.cross(approx);
            let k = co.c.curvature(t);
            let k_off = k / (1. + k * self.d);
            let adj_err = dot - 1. + 0.5 * k_off * self.d * cross * cross;
            err = adj_err.abs().max(err);
            web_sys::console::log_1(&format!("{i}: {:.3} {adj_err:.3} {k:.4}", dot - 1.).into());
        }
        1.2 * err
    }

    pub fn est_arc_error_dot_no_k(&self, co: &CurveOffset) -> f64 {
        let delta = co.make_delta(self.a, self.b);
        const N: usize = 8;
        let mut err = 0.0;
        for i in 0..N {
            let t = (i as f64 + 0.5) * (1. / N as f64);
            let unorm = turn(co.q.eval(t).to_vec2().normalize());
            let approx = delta.eval(t).to_vec2();
            let dot = unorm.dot(approx);
            let adj_err = dot - 1.;
            err = adj_err.abs().max(err);
        }
        1.2 * err
    }
}

#[test]
fn foo() {
    use kurbo::Vec2;
    let x = Vec2::new(1.0, 0.0);
    let y = Vec2::new(0.0, 1.0);
    println!("{}", x.cross(y));
}
