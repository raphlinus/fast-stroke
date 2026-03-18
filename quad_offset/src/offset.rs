use xilem_web::svg::kurbo::{CubicBez, ParamCurve, QuadBez, Vec2};

/// Unit normal (rotated 90° ccw from tangent) of a QuadBez at parameter t.
fn quad_normal(q: &QuadBez, t: f64) -> Vec2 {
    let d0 = q.p1 - q.p0;
    let d1 = q.p2 - q.p1;
    let tangent = (1.0 - t) * d0 + t * d1;
    tangent.normalize().turn_90()
}

/// Fit a cubic Bezier to the offset curve using one-point shape control.
///
/// Endpoints and tangent directions match the true offset, and the midpoint
/// t=0.5 is interpolated exactly.
pub fn offset_cubic(q: &QuadBez, d: f64) -> CubicBez {
    let utan0 = (q.p1 - q.p0).normalize();
    let utan1 = (q.p2 - q.p1).normalize();

    let p0 = q.p0 + d * utan0.turn_90();
    let p3 = q.p2 + d * utan1.turn_90();
    let m = q.eval(0.5) + d * (q.p2 - q.p0).normalize().turn_90();

    let rhs = m - p0.midpoint(p3);
    let det = utan1.cross(utan0);
    let (p1, p2) = if det.abs() > 1e-10 {
        let idet = (8. / 3.) / det;
        let a = (utan1.cross(rhs) * idet).max(0.0);
        let b = (utan0.cross(rhs) * idet).max(0.0);
        (p0 + a * utan0, p3 - b * utan1)
    } else {
        let chord = p3 - p0;
        (p0 + chord * (1. / 3.), p3 - chord * (1. / 3.))
    };
    CubicBez::new(p0, p1, p2, p3)
}

/// Approximate the offset of a QuadBez by offsetting the endpoints and setting
/// the control point to the intersection of the endpoint tangents.
pub fn offset_quad(q: &QuadBez, d: f64) -> QuadBez {
    let p0 = q.p0 + d * quad_normal(q, 0.0);
    let p2 = q.p2 + d * quad_normal(q, 1.0);
    // Tangent directions at t=0 and t=1 (same as the original curve).
    let t0 = q.p1 - q.p0;
    let t1 = q.p2 - q.p1;
    // Intersect: p0 + s*t0 = p2 + u*t1  =>  s*t0 - u*t1 = p2 - p0
    let rhs = p2 - p0;
    let det = t0.cross(t1);
    let p1 = if det.abs() > 1e-10 {
        let s = rhs.cross(t1) / det;
        (p0.to_vec2() + s * t0).to_point()
    } else {
        // Parallel tangents: fall back to midpoint of the chord.
        (0.5 * (p0.to_vec2() + p2.to_vec2())).to_point()
    };
    QuadBez::new(p0, p1, p2)
}
