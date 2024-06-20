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
