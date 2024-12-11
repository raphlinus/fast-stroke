use kurbo::{Affine, BezPath, CubicBez, ParamCurve, ParamCurveDeriv, Point, Vec2};

fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}

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
        (ca + cb + cc) / q.eval(t).to_vec2().hypot()
    };
    const N: usize = 20;
    for i in 0..=N {
        let t = i as f64 / N as f64;
        let y = 5e-2 * error(t);
        let p = Affine::translate((100., 200.)) * Affine::scale(500.) * Point::new(t, y);
        if i == 0 {
            result.move_to(p);
        } else {
            result.line_to(p)
        }
    }
    result
}

#[test]
fn foo() {
    use kurbo::Vec2;
    let x = Vec2::new(1.0, 0.0);
    let y = Vec2::new(0.0, 1.0);
    println!("{}", x.cross(y));
}
