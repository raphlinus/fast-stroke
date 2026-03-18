use xilem_web::svg::kurbo::{CubicBez, ParamCurve, QuadBez, Vec2};

/// Unit normal (rotated 90° ccw from tangent) of a QuadBez at parameter t.
fn quad_normal(q: &QuadBez, t: f64) -> Vec2 {
    let d0 = q.p1 - q.p0;
    let d1 = q.p2 - q.p1;
    let tangent = 2.0 * ((1.0 - t) * d0 + t * d1);
    tangent.normalize().turn_90()
}

/// Fit a cubic Bezier to the offset curve using one-point shape control:
/// endpoints and tangent directions match the true offset, and the midpoint
/// t=0.5 is interpolated exactly.
///
/// With B = A + α·n0 and C = D − β·n1, the midpoint condition
/// C(0.5) = M gives α·n0 − β·n1 = (8M − 4A − 4D)/3, a 2×2 linear system.
pub fn offset_cubic(q: &QuadBez, d: f64) -> CubicBez {
    let a = q.p0 + d * quad_normal(q, 0.0);
    let dd = q.p2 + d * quad_normal(q, 1.0);
    let m = q.eval(0.5) + d * quad_normal(q, 0.5);

    let n0 = (q.p1 - q.p0).normalize();
    let n1 = (q.p2 - q.p1).normalize();

    let rhs = (8.0 * m.to_vec2() - 4.0 * a.to_vec2() - 4.0 * dd.to_vec2()) / 3.0;
    let det = n1.cross(n0);
    let (b, c) = if det.abs() > 1e-10 {
        let alpha = n1.cross(rhs) / det;
        let beta = n0.cross(rhs) / det;
        (a + alpha * n0, dd - beta * n1)
    } else {
        // Parallel tangents (straight-ish curve): space control points uniformly.
        let chord = dd - a;
        (a + chord / 3.0, a + 2.0 * chord / 3.0)
    };
    CubicBez::new(a, b, c, dd)
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
