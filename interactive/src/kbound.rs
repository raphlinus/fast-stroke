use kurbo::{CubicBez, ParamCurveDeriv};

/// Compute bounds on curvature
pub fn kbound(c: CubicBez) -> (f64, f64) {
    let q = c.deriv();
    let p1xp0 = q.p1.to_vec2().cross(q.p0.to_vec2());
    let p2xp0 = q.p2.to_vec2().cross(q.p0.to_vec2());
    let p2xp1 = q.p2.to_vec2().cross(q.p1.to_vec2());
    let c0 = 2. * p1xp0;
    let c1 = 2. * (p2xp0 - 2.0 * p1xp0);
    let c2 = 2. * (p2xp1 - p2xp0 + p1xp0);
    let (num_min, num_max) = quadratic_min_max(c0, c1, c2);
    let norm_chord = (c.p3 - c.p0).normalize();
    let r0 = q.p0.to_vec2().dot(norm_chord);
    let r1 = q.p1.to_vec2().dot(norm_chord);
    let r2 = q.p2.to_vec2().dot(norm_chord);
    let d0 = r0;
    let d1 = 2. * (r1 - r0);
    let d2 = (r2 - r1) - (r1 - r0);
    let (deriv_min, _deriv_max) = quadratic_min_max(d0, d1, d2);
    if deriv_min <= 0.0 {
        return (-f64::INFINITY, f64::INFINITY);
    }
    let scale = 1.0 / (deriv_min * deriv_min * deriv_min);
    (num_min.min(0.) * scale, num_max.max(0.) * scale)
}

fn quadratic_min_max(c0: f64, c1: f64, c2: f64) -> (f64, f64) {
    let y1 = c0 + c1 + c2;
    let mut min = c0.min(y1);
    let mut max = c0.max(y1);
    if c2 * c1 < 0. && c2 * (c1 + 2.0 * c2) > 0. {
        let x = -c1 / (2.0 * c2);
        let y = c0 + x * (c1 + x * c2);
        min = min.min(y);
        max = max.max(y);
    }
    (min, max)
}

const N: usize = 1;

/// Compute bounds on curvature
pub fn kbound_refined(c: CubicBez) -> (f64, f64) {
    let q = c.deriv();
    let p1xp0 = q.p1.to_vec2().cross(q.p0.to_vec2());
    let p2xp0 = q.p2.to_vec2().cross(q.p0.to_vec2());
    let p2xp1 = q.p2.to_vec2().cross(q.p1.to_vec2());
    let c0 = 2. * p1xp0;
    let c1 = 2. * (p2xp0 - 2.0 * p1xp0);
    let c2 = 2. * (p2xp1 - p2xp0 + p1xp0);
    let (num_min, num_max) = q_min_max_refined(c0, c1, c2);
    let norm_chord = (c.p3 - c.p0).normalize();
    let r0 = q.p0.to_vec2().dot(norm_chord);
    let r1 = q.p1.to_vec2().dot(norm_chord);
    let r2 = q.p2.to_vec2().dot(norm_chord);
    let d0 = r0;
    let d1 = 2. * (r1 - r0);
    let d2 = (r2 - r1) - (r1 - r0);
    let (deriv_min, _deriv_max) = q_min_max_refined(d0, d1, d2);
    let mut kmin = 0.0;
    let mut kmax = 0.0;
    web_sys::console::log_1(&format!("dmin {deriv_min:?}").into());
    web_sys::console::log_1(&format!("nmin {num_min:?}").into());
    for i in 0..N {
        let dmin = deriv_min[i];
        if dmin <= 0.0 {
            return (-f64::INFINITY, f64::INFINITY);
        }
        let scale = 1.0 / (dmin * dmin * dmin);
        kmin = (num_min[i] * scale).min(kmin);
        kmax = (num_max[i] * scale).max(kmax);
    }
    (kmin, kmax)
}

fn q_min_max_refined(c0: f64, c1: f64, c2: f64) -> ([f64; N], [f64; N]) {
    let mut min = [0.0; N];
    let mut max = [0.0; N];
    let mut last = c0;
    for i in 0..N {
        let x = (i + 1) as f64 * (1.0 / N as f64);
        let y = c0 + x * (c1 + x * c2);
        min[i] = last.min(y);
        max[i] = last.max(y);
        last = y;
    }
    let x = -c1 / (2.0 * c2);
    if x > 0.0 && x < 1.0 {
        let y = c0 + x * (c1 + x * c2);
        let i = (x * N as f64).floor() as usize;
        min[i] = min[i].min(y);
        max[i] = max[i].max(y);
    }

    (min, max)
}
