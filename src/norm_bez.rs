use kurbo::{CubicBez, Point};

/// Normalize a cubic Bézier to the unit chord.
///
/// Returns the length of the chord.
pub fn normalize_bez(c: CubicBez) -> (CubicBez, f64) {
    let chord = c.p3 - c.p0;
    let scale = 1.0 / chord.hypot2();
    let u = chord.x * scale;
    let v = chord.y * scale;
    let d01 = c.p1 - c.p0;
    let p1 = Point::new(d01.x * u + d01.y * v, d01.y * u - d01.x * v);
    let d02 = c.p2 - c.p0;
    let p2 = Point::new(d02.x * u + d02.y * v, d02.y * u - d02.x * v);
    let chord_len = chord.hypot();
    (
        CubicBez::new(Point::new(0., 0.), p1, p2, Point::new(1., 0.)),
        chord_len,
    )
}

/// This function is the inverse of `params_to_cubic`.
pub fn normed_bez_to_params(c: CubicBez) -> [f64; 4] {
    let th0 = c.p1.y.atan2(c.p1.x);
    let th1 = c.p2.y.atan2(1.0 - c.p2.x);
    let d0 = c.p1.y.hypot(c.p1.x);
    let d1 = c.p2.y.hypot(1.0 - c.p2.x);
    let e0_inv = 1.5 * (1.0 + th0.cos());
    let e1_inv = 1.5 * (1.0 + th1.cos());
    let a = 0.5 * (th0 + th1);
    let b = 0.5 * (th0 - th1);
    let d0_normed = d0 * e0_inv - 1.0;
    let d1_normed = d1 * e1_inv - 1.0;
    let c = 0.5 * (d0_normed + d1_normed);
    let d = 0.5 * (d0_normed - d1_normed);
    [a, b, c, d]
}
