//! New module for midpoint G1 solutions

use kurbo::{
    common::solve_itp, BezPath, CubicBez, ParamCurve, ParamCurveCurvature, ParamCurveDeriv, Point,
    QuadBez, Vec2,
};

// We can use the Vec2 method when we upgrade kurbo
use crate::perturb::turn;

/// Info that's constant for the curve.
///
/// This is basically the same as offset.
pub struct CubicOffset {
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
}

/// Info specific to a subdivision.
pub struct OffsetRec {
    t0: f64,
    t1: f64,
    // unit tangent at t0
    utan0: Vec2,
    // unit tangent at t1
    utan1: Vec2,
    cusp0: f64,
    cusp1: f64,
    depth: usize,
    idet: f64,
    c_base: Point,
    tana_base: Vec2,
    extreme: bool,
}

/// Intermediate results from sampling the curve at some t
struct Sample {
    t: f64,
    a: f64,
    b: f64,
    tan: Vec2,
    utan: Vec2,
    tana: Vec2,
    angle_err: f64,
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
    let (cusp0, utan0) = co.cusp_and_utan_inner(co.q.p0, co.c0);
    let (cusp1, utan1) = co.cusp_and_utan_inner(co.q.p2, co.c0 + co.c1 + co.c2);
    result.move_to(c.p0 + d * turn(utan0));
    let rec = OffsetRec::new(&co, 0., 1., utan0, utan1, cusp0, cusp1, 0, false);
    co.offset_rec(&rec, &mut result);
    result
}

impl CubicOffset {
    /// Create a new curve from Bézier segment and offset.
    pub fn new(c: CubicBez, d: f64, tolerance: f64) -> Self {
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

    pub fn approximate(&self) -> CubicBez {
        let (c, _err, _t) = self.compute_rec(&self.init_rec());
        web_sys::console::log_1(&format!("err = {}", _err.sqrt()).into());
        c
    }

    /// Compute one recursion.
    ///
    /// Return value: cubic, error squared, t
    fn compute_rec(&self, rec: &OffsetRec) -> (CubicBez, f64, f64) {
        let (a, b, t) = if let Some(sample) = self.solve_midpoint(&rec) {
            (sample.a, sample.b, sample.t)
        } else {
            let (a, b) = self.arc_draw(&rec);
            (a, b, 0.5 * (rec.t0 + rec.t1))
        };
        let c = self.apply(&rec, a, b);
        let err = self.eval_err(&rec, t, c);
        web_sys::console::log_1(&format!("err = {}", err.sqrt()).into());
        (c, err, t)
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
            self.subdivide(rec, result, t, utan_t, cusp_t_minus, cusp_t_plus, false);
            return;
        }
        let (c, err, mut t) = self.compute_rec(rec);
        let dt = rec.t1 - rec.t0;
        let mut child_extreme = t < rec.t0 + 0.1 * dt || t > rec.t0 + 0.9 * dt;
        if rec.extreme && child_extreme {
            t = 0.5 * (rec.t0 + rec.t1);
            child_extreme = false;
        }

        if rec.depth < MAX_DEPTH && err > self.tolerance.powi(2) {
            let (cusp, utan) = self.cusp_and_utan(t);
            self.subdivide(rec, result, t, utan, cusp, cusp, child_extreme);
        } else {
            result.curve_to(c.p1, c.p2, c.p3);
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
        extreme: bool,
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
            extreme,
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
            extreme,
        );
        self.offset_rec(&rec1, result);
    }

    fn init_rec(&self) -> OffsetRec {
        let (cusp0, utan0) = self.cusp_and_utan_inner(self.q.p0, self.c0);
        let (cusp1, utan1) = self.cusp_and_utan_inner(self.q.p2, self.c0 + self.c1 + self.c2);
        OffsetRec::new(self, 0., 1.0, utan0, utan1, cusp0, cusp1, 0, false)
    }

    fn cusp_sign(&self, t: f64) -> f64 {
        let ds2 = self.q.eval(t).to_vec2().hypot2();
        ((self.c2 * t + self.c1) * t + self.c0) / (ds2 * ds2.sqrt()) + 1.0
    }

    fn cusp_and_utan(&self, t: f64) -> (f64, Vec2) {
        self.cusp_and_utan_inner(self.q.eval(t), (self.c2 * t + self.c1) * t + self.c0)
    }

    fn cusp_and_utan_inner(&self, tan: Point, y: f64) -> (f64, Vec2) {
        let rsqrt = 1.0 / tan.to_vec2().hypot();
        let cusp = y * (rsqrt * rsqrt * rsqrt) + 1.0;
        let utan = rsqrt * tan.to_vec2();
        (cusp, utan)
    }

    fn apply(&self, rec: &OffsetRec, a: f64, b: f64) -> CubicBez {
        // wondering if p0 and p3 should be in rec
        let p0 = self.c.eval(rec.t0) + self.d * turn(rec.utan0);
        let l0 = a;
        let mut p1 = p0;
        if l0 * rec.cusp0 > 0.0 {
            p1 += l0 * rec.utan0;
        }
        let p3 = self.c.eval(rec.t1) + self.d * turn(rec.utan1);
        let mut p2 = p3;
        let l1 = -b;
        if l1 * rec.cusp1 > 0.0 {
            p2 -= l1 * rec.utan1;
        }
        CubicBez::new(p0, p1, p2, p3)
    }

    fn one_point_at(&self, rec: &OffsetRec, t: f64) -> Sample {
        let tan = self.q.eval(t).to_vec2();
        let utan = tan.normalize();
        let z = self.c.eval(t) + self.d * turn(utan) - rec.c_base;
        let a = z.cross(rec.utan1) * rec.idet;
        let b = rec.utan0.cross(z) * rec.idet;
        let tana = rec.tana_base + 0.75 * (b * rec.utan1 - a * rec.utan0);
        let angle_err = tana.cross(utan);
        Sample {
            t,
            a,
            b,
            tan,
            utan,
            tana,
            angle_err,
        }
    }

    fn angle_err_deriv(&self, rec: &OffsetRec, sample: &Sample) -> f64 {
        // Multiply these by 3/8 (1 + k d) / det to get true da/dt and db/dt
        let dadt = sample.tan.cross(rec.utan1);
        let dbdt = rec.utan0.cross(sample.tan);
        // TODO: redo this in terms of cusp etc
        let k = self.c.curvature(sample.t);
        //let dutan_dt = -k * turn(sample.tan);
        let factor = 0.75 * (1.0 + k * self.d) * rec.idet;
        let dtana_dt = dbdt * rec.utan1 - dadt * rec.utan0;
        //sample.tana.cross(dutan_dt) + factor * dtana_dt.cross(sample.utan)
        -k * sample.tan.dot(sample.tana) + factor * dtana_dt.cross(sample.utan)
    }

    fn solve_midpoint(&self, rec: &OffsetRec) -> Option<Sample> {
        if let Some(sample) = self.try_midpoint_newton(rec, 0.5 * (rec.t0 + rec.t1)) {
            return Some(sample);
        }
        let roots = self.midpoint_cubic_approx(rec);
        for root in &roots {
            if let Some(sample) = self.try_midpoint_newton(rec, *root) {
                return Some(sample);
            }
        }
        None
    }

    fn try_midpoint_newton(&self, rec: &OffsetRec, t0: f64) -> Option<Sample> {
        // This is a somewhat delicate setting. Observed error tends to be
        // about 1/3 of angle_err (not rigorously validated). If the threshold
        // is set too high, then out-of-tolerance approximations will result,
        // which might have been in tolerance had there been another iteration.
        // If it is set too low, then in addition to wasted work, there is the
        // risk that the initial iteration was good enough but Newton iteration
        // would not converge.
        let threshold = 3. * self.tolerance;
        const MAX_ITERS: usize = 10;
        let mut t = t0;
        let mut sample = self.one_point_at(rec, t);
        for _i in 0..MAX_ITERS {
            web_sys::console::log_1(&format!("newton {_i}: {}", sample.angle_err).into());
            let err = sample.angle_err.abs();
            if err < threshold {
                if sample.a * rec.cusp0 < 0. || sample.b * rec.cusp1 > 0. {
                    web_sys::console::log_1(&format!("rejected {}", sample.t).into());
                    break;
                }
                return Some(sample);
            }
            t -= sample.angle_err / self.angle_err_deriv(rec, &sample);
            if !(t >= rec.t0 && t <= rec.t1) {
                break;
            }
            sample = self.one_point_at(rec, t);
            if sample.angle_err.abs() > err {
                break;
            }
        }
        web_sys::console::log_1(&format!("newton failure t0={t0}").into());
        None
    }

    /// Find cubic approximation of midpoint G1 function and solve.
    ///
    /// Results are filtered to have positive control point distances
    fn midpoint_cubic_approx(&self, rec: &OffsetRec) -> Vec<f64> {
        let t0 = 0.5 * (rec.t0 + rec.t1);
        let dt = 1e-3 * (rec.t1 - rec.t0);
        // TODO: harvest s0 and d0 from initial try
        let s0 = self.one_point_at(rec, t0);
        let s1 = self.one_point_at(rec, t0 + dt);
        let d0 = dt * self.angle_err_deriv(rec, &s0);
        let d1 = dt * self.angle_err_deriv(rec, &s1);
        let c0 = s0.angle_err;
        let c1 = d0;
        let derr = s1.angle_err - s0.angle_err;
        let c2 = 3. * derr - 2. * d0 - d1;
        let c3 = -2. * derr + d0 + d1;
        let roots = kurbo::common::solve_cubic(c0, c1, c2, c3);
        let mut sorted = roots.iter().map(|x| t0 + x * dt).collect::<Vec<_>>();
        sorted.sort_by_key(|x| (x - (t0 + 0.5 * dt)).abs().to_bits());
        web_sys::console::log_1(&format!("{sorted:?}").into());
        sorted
    }

    // Evaluate the error of the approximation
    fn eval_err(&self, rec: &OffsetRec, tm: f64, c: CubicBez) -> f64 {
        let mut max_err = 0.0;
        for t in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] {
            let p_approx = c.eval(t);
            let q_approx = c.deriv().eval(t).to_vec2();
            let mut u = if t < 0.5 {
                rec.t0 + 2.0 * t * (tm - rec.t0)
            } else {
                rec.t1 + 2.0 * (1.0 - t) * (tm - rec.t1)
            };
            // bring p into normal ray of approx
            // TODO: make iteration count adaptive
            const N: usize = 3;
            for i in 0..N {
                let tan = self.q.eval(u).to_vec2();
                let p = self.c.eval(u) + self.d * turn(tan.normalize());
                let dpdu = tan * self.cusp_sign(u);
                u -= (p - p_approx).dot(q_approx) / dpdu.dot(q_approx);
                let err2 = p_approx.distance_squared(p);
                if i == N - 1 {
                    let z = q_approx.cross(dpdu);
                    web_sys::console::log_1(&format!("{t}: {} {:.3}", err2.sqrt(), z).into());
                    max_err = err2.max(max_err);
                }
            }
        }
        max_err
    }

    fn arc_draw(&self, rec: &OffsetRec) -> (f64, f64) {
        let sum = rec.utan0 + rec.utan1;
        // TODO: tune threshold to optimize numerical precision
        let d = if rec.idet.abs() < 1e9 {
            0.5 * (sum.length() - sum.dot(rec.utan0)) * rec.idet
        } else {
            (1. / 3.) * rec.utan0.cross(rec.utan1)
        };
        let scale = (1. / 3.) * (rec.t1 - rec.t0);
        let a = scale * self.q.eval(rec.t0).to_vec2().length() - d * self.d;
        let b = -scale * self.q.eval(rec.t1).to_vec2().length() + d * self.d;
        (a, b)
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
        extreme: bool,
    ) -> Self {
        let idet = (8. / 3.) / utan0.cross(utan1);
        // Approximation and deriv at approx t=0.5 when (a, b) are zero
        let p0 = co.c.eval(t0);
        let p1 = co.c.eval(t1);
        let c_base = p0.midpoint(p1) + co.d * 0.5 * turn(utan0 + utan1);
        let tana_base = 1.5 * ((p1 - p0) + co.d * turn(utan1 - utan0));
        OffsetRec {
            t0,
            t1,
            utan0,
            utan1,
            cusp0,
            cusp1,
            depth,
            idet,
            c_base,
            tana_base,
            extreme,
        }
    }
}
