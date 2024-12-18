// Doing lots of experiments, will try things and then move on.
#[allow(unused)]
use kurbo::{
    common::solve_cubic, Affine, BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, Vec2,
};

fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}

const ERROR_SCALE: f64 = 5.0;

// Produce the delta for the given curve
pub fn linear_approx(c: CubicBez) -> CubicBez {
    let q = c.deriv();
    let b01 = q.p0.to_vec2();
    let b12 = q.p1.to_vec2();
    let b23 = q.p2.to_vec2();
    let n0 = turn(b01).normalize();
    let n1 = turn(b23).normalize();
    let coefs = |t: f64| {
        let mt = 1.0 - t;
        let ca =
            6.0 * mt.powi(3) * t * t * b01.cross(b12) + 3.0 * mt * mt * t.powi(3) * b01.cross(b23);
        let cb =
            3.0 * mt.powi(3) * t * t * b23.cross(b01) + 6.0 * mt * mt * t.powi(3) * b23.cross(b12);
        let cc = mt.powi(4) * (mt + 3.0 * t) * n0.cross(b01)
            + mt * mt * t * t * (3.0 * mt + t) * n1.cross(b01)
            + mt.powi(3) * t * (2.0 * mt + 6.0 * t) * n0.cross(b12)
            + mt * t.powi(3) * (6.0 * mt + 2.0 * t) * n1.cross(b12)
            + mt * mt * t * t * (mt + 3.0 * t) * n0.cross(b23)
            + t.powi(4) * (3.0 * mt + t) * n1.cross(b23);
        (ca, cb, cc)
    };
    let t0 = 1.0 / 3.0;
    let t1 = 2.0 / 3.0;
    let (ca0, cb0, cc0) = coefs(t0);
    let (ca1, cb1, cc1) = coefs(t1);
    let z0 = -q.eval(t0).to_vec2().hypot() - cc0;
    let z1 = -q.eval(t1).to_vec2().hypot() - cc1;
    let det = ca0 * cb1 - ca1 * cb0;
    let a = (z0 * cb1 - z1 * cb0) / det;
    let b = (ca0 * z1 - ca1 * z0) / det;
    CubicBez::new(
        n0.to_point(),
        (n0 + a * b01).to_point(),
        (n1 + b * b23).to_point(),
        n1.to_point(),
    )
}

pub fn plot_error(c: CubicBez, delta: CubicBez) -> BezPath {
    let mut result = BezPath::new();
    let q = c.deriv();
    let b01 = q.p0.to_vec2();
    let b12 = q.p1.to_vec2();
    let b23 = q.p2.to_vec2();
    let n0 = turn(b01).normalize();
    let n1 = turn(b23).normalize();
    let a = delta.p1 - delta.p0;
    let b = delta.p2 - delta.p3;
    let error = |t: f64| {
        let mt = 1.0 - t;
        let ca = 6.0 * mt.powi(3) * t * t * a.cross(b12) + 3.0 * mt * mt * t.powi(3) * a.cross(b23);
        let cb = 3.0 * mt.powi(3) * t * t * b.cross(b01) + 6.0 * mt * mt * t.powi(3) * b.cross(b12);
        let cc = mt.powi(4) * (mt + 3.0 * t) * n0.cross(b01)
            + mt * mt * t * t * (3.0 * mt + t) * n1.cross(b01)
            + mt.powi(3) * t * (2.0 * mt + 6.0 * t) * n0.cross(b12)
            + mt * t.powi(3) * (6.0 * mt + 2.0 * t) * n1.cross(b12)
            + mt * mt * t * t * (mt + 3.0 * t) * n0.cross(b23)
            + t.powi(4) * (3.0 * mt + t) * n1.cross(b23);
        -(ca + cb + cc) / q.eval(t).to_vec2().hypot() - 1.0
    };
    const N: usize = 50;
    for i in 0..=N {
        let t = i as f64 / N as f64;
        let y = ERROR_SCALE * error(t);
        let p = Affine::translate((100., 200.)) * Affine::scale(500.) * Point::new(t, y);
        if i == 0 {
            result.move_to(p);
        } else {
            result.line_to(p);
        }
    }
    result
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
        let delta = linear_approx(c_s);
        let c_offset = CubicBez::new(
            c_s.p0 + d * delta.p0.to_vec2(),
            c_s.p1 + d * delta.p1.to_vec2(),
            c_s.p2 + d * delta.p2.to_vec2(),
            c_s.p3 + d * delta.p3.to_vec2(),
        );
        let mut max_err = 0.0;
        for (_, e) in error_by_rays(c_s, d, c_offset) {
            max_err = e.max(max_err);
        }
        errs.push(max_err);
    }
    web_sys::console::log_1(&format!("{:?}", errs).into());
}

/// Linear approximation to error minmax
///
/// The basic idea is to make e(1/6) + e(1/2) and e(1/2) + e(5/6) both
/// zero. This roughly minimizes the maximum error for the first Chebyshev
/// polynomial that can't be zeroed exactly.
pub fn linear_minmax(c: CubicBez) -> CubicBez {
    let q = c.deriv();
    let b01 = q.p0.to_vec2();
    let b12 = q.p1.to_vec2();
    let b23 = q.p2.to_vec2();
    let n0 = turn(b01).normalize();
    let n1 = turn(b23).normalize();
    let coefs = |t: f64| {
        let mt = 1.0 - t;
        let ca =
            6.0 * mt.powi(3) * t * t * b01.cross(b12) + 3.0 * mt * mt * t.powi(3) * b01.cross(b23);
        let cb =
            3.0 * mt.powi(3) * t * t * b23.cross(b01) + 6.0 * mt * mt * t.powi(3) * b23.cross(b12);
        let cc = mt.powi(4) * (mt + 3.0 * t) * n0.cross(b01)
            + mt * mt * t * t * (3.0 * mt + t) * n1.cross(b01)
            + mt.powi(3) * t * (2.0 * mt + 6.0 * t) * n0.cross(b12)
            + mt * t.powi(3) * (6.0 * mt + 2.0 * t) * n1.cross(b12)
            + mt * mt * t * t * (mt + 3.0 * t) * n0.cross(b23)
            + t.powi(4) * (3.0 * mt + t) * n1.cross(b23);
        (ca, cb, cc)
    };
    let t0 = 1.0 / 6.0;
    let t1 = 1.0 / 2.0;
    let t2 = 5.0 / 6.0;
    let (ca0, cb0, cc0) = coefs(t0);
    let (ca1, cb1, cc1) = coefs(t1);
    let (ca2, cb2, cc2) = coefs(t2);
    let z0 = -q.eval(t0).to_vec2().hypot() - cc0;
    let z1 = -q.eval(t1).to_vec2().hypot() - cc1;
    let z2 = -q.eval(t2).to_vec2().hypot() - cc2;
    let ca01 = ca0 + ca1;
    let cb01 = cb0 + cb1;
    let z01 = z0 + z1;
    let ca12 = ca1 + ca2;
    let cb12 = cb1 + cb2;
    let z12 = z1 + z2;
    let det = ca01 * cb12 - ca12 * cb01;
    let a = (z01 * cb12 - z12 * cb01) / det;
    let b = (ca01 * z12 - ca12 * z01) / det;
    CubicBez::new(
        n0.to_point(),
        (n0 + a * b01).to_point(),
        (n1 + b * b23).to_point(),
        n1.to_point(),
    )
}

#[test]
fn foo() {
    use kurbo::Vec2;
    let x = Vec2::new(1.0, 0.0);
    let y = Vec2::new(0.0, 1.0);
    println!("{}", x.cross(y));
}
