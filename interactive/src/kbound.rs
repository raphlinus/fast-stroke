use kurbo::{Affine, CubicBez, ParamCurve, ParamCurveCurvature, ParamCurveDeriv, Shape};

use crate::perturb::turn;

/// Compute bounds on curvature
pub fn kbound(c: CubicBez) -> (f64, f64) {
    let q = c.deriv();
    let dd0 = q.p1 - q.p0;
    let dd1 = q.p2 - q.p1;
    let dd2 = dd0.hypot2().max(dd1.hypot2());
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
    let dmin_inv = 1. / deriv_min;
    let dmin_m2 = dmin_inv * dmin_inv;
    let alt = 2.0 * dd2.sqrt();
    let kmin = (num_min.min(0.) * dmin_inv).max(-alt) * dmin_m2;
    let kmax = (num_max.max(0.) * dmin_inv).min(alt) * dmin_m2;
    (kmin, kmax)
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

pub fn kbound_accurate(c: CubicBez) -> (f64, f64) {
    const N: usize = 8;
    let mut min = 0.0;
    let mut max = 0.0;
    for i in 0..N {
        let t0 = (i as f64) * (1.0 / N as f64);
        let t1 = t0 + 1.0 / N as f64;
        let (kmin, kmax) = kbound(c.subsegment(t0..t1));
        max = kmax.max(max);
        min = kmin.min(min);
    }
    (min, max)
}

fn kbound_accurate_abs(c: CubicBez) -> f64 {
    const N: usize = 16;
    let mut max = 0.0;
    for i in 0..N {
        let t0 = (i as f64) * (1.0 / N as f64);
        let t1 = t0 + 1.0 / N as f64;
        let (kmin, kmax) = kbound(c.subsegment(t0..t1));
        max = kmax.max(-kmin).max(max);
    }
    max
}

pub fn est_arc_error_k4(c: CubicBez, d: f64) -> f64 {
    let k = kbound_accurate_abs(c);
    let chord = (c.p3 - c.p0).hypot();
    let ratio = k * chord;
    d * ratio.powi(4)
}

// Other estimation stuff goes here

/// Estimate error based on cubic parameters
pub fn est_arc_error(c: CubicBez, d: f64) -> f64 {
    let c_trans = Affine::translate(-c.p0.to_vec2()) * c;
    let scale = 1.0 / c_trans.p3.to_vec2().hypot2();
    let u = c_trans.p3.x * scale;
    let v = c_trans.p3.y * scale;
    let c_norm = Affine::new([u, -v, v, u, 0., 0.]) * c_trans;
    let th0 = c_norm.p1.y.atan2(c_norm.p1.x);
    let th1 = c_norm.p2.y.atan2(1. - c_norm.p2.x);
    let y_minus = c_norm.p2.y - c_norm.p1.y;
    let y_plus = c_norm.p2.y + c_norm.p1.y;
    let x1 = c_norm.p1.x + 0.75 * c_norm.p1.y.powi(2) - 1. / 3.;
    let x2 = (1.0 - c_norm.p2.x) + 0.75 * c_norm.p2.y.powi(2) - 1. / 3.;
    web_sys::console::log_1(&format!("xs {x1} {x2}").into());
    let dk_err = (th1 - th0) * (th1 - th0);
    let k2 = (th1 + th0) * (th1 + th0);
    let x_plus_err = k2 * (x1 + x2) * (x1 + x2);
    let x_minus_err = k2 * (x1 + x2) * (x1 + x2);
    let k6_err = k2.powi(3);
    let xs = x1 * x1 + x2 * x2;
    let x1_err = k2 * x1 * x1;
    let x2_err = k2 * x2 * x2;
    dk_err * 1. + 1e-4 * k6_err + 3. * k2 * xs
}

pub fn est_arc_error_dot(c: CubicBez, d: f64) -> f64 {
    let q = c.deriv();
    let unorm0 = q.p0.to_vec2().normalize();
    let unorm1 = q.p2.to_vec2().normalize();
    let cross = unorm0.cross(unorm1);
    let sum = unorm0 + unorm1;
    // TODO: numerical robustness when cross is small, see perturb arc_draw
    let d = (4. / 3.) * (sum.length() - sum.dot(unorm0)) / unorm0.cross(unorm1);
    let p0 = turn(unorm0).to_point();
    let p3 = turn(unorm1).to_point();
    let p1 = p0 - d * unorm0;
    let p2 = p3 + d * unorm1;
    let delta = CubicBez::new(p0, p1, p2, p3);
    const N: usize = 4;
    let mut err = 0.0;
    for i in 0..N {
        let t = (i as f64 + 0.5) * (1. / N as f64);
        let dot = turn(q.eval(t).to_vec2().normalize()).dot(delta.eval(t).to_vec2());
        let base_err = dot - 1.;
        let k = c.curvature(t);
        let adj_err = base_err / (1. - k * d);
        err = adj_err.abs().max(err);
        web_sys::console::log_1(&format!("{i}: {dot}").into());
    }
    5. * err
}
