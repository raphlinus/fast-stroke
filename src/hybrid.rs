//! This is the hybrid approach to parallel curves, intended for very high performance.

use kurbo::{BezPath, CubicBez, ParamCurveDeriv, QuadBez, Vec2};

/// Derivatives of source curve.
struct Derivs {
    q: QuadBez,
    // The cross product of second and first derivatives is the polynomial
    // c0 + c1 * t + c2 * t^2
    c0: f64,
    c1: f64,
    c2: f64,
}

fn parallel_cubic(c: CubicBez, tolerance: f64, out: &mut BezPath) {
    // First step: determine whether the Bézier has a single inflection point; try one-point only
    // if not.
    let derivs = Derivs::from_cubic(c);
    let d2xd1_0 = derivs.c0;
    // Note: we can micro-optimize this, it's 2 * p2xp1
    let d2xd1_1 = derivs.c0 + derivs.c1 + derivs.c2;
    if d2xd1_0 * d2xd1_1 >= 0.0 {
        // not a single inflection point, ok to try one-point
    }
}

impl Derivs {
    fn from_cubic(c: CubicBez) -> Self {
        let q = c.deriv();
        let p1xp0 = q.p1.to_vec2().cross(q.p0.to_vec2());
        let p2xp0 = q.p2.to_vec2().cross(q.p0.to_vec2());
        let p2xp1 = q.p2.to_vec2().cross(q.p1.to_vec2());
        // Possible micro-optimization: don't do multiplication by 2.0 here,
        // fold that into downstream calculations. Not needed for detecting
        // the presence of inflection points.
        Derivs {
            q,
            c0: 2.0 * p1xp0,
            c1: 2.0 * (p2xp0 - 2.0 * p1xp0),
            c2: 2.0 * (p2xp1 - p2xp0 + p1xp0),
        }
    }
}

// TODO: use turn_90 method when upgrade kurbo
pub fn turn(v: Vec2) -> Vec2 {
    Vec2::new(-v.y, v.x)
}
