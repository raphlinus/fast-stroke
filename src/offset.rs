use kurbo::{CubicBez, Point, Vec2};

use crate::{
    approx::{self, approx_partials},
    cheb::Cheb,
    norm_bez, params_to_cubic,
};

/// Approximate a parallel curve.
pub fn approx_offset(c: CubicBez, offset: f64) -> CubicBez {
    let params = norm_bez::normed_bez_to_params(c);
    let th0 = params[0] + params[1];
    let th1 = params[0] - params[1];
    let (s0, c0) = th0.sin_cos();
    let (s1, c1) = th1.sin_cos();
    let p0 = Point::new(-s0 * offset, c0 * offset);
    let p1 = Point::new(1.0 + s1 * offset, c1 * offset);
    let ch = p1 - p0;
    let dth = ch.atan2();

    let offset_cheb = approx::approx_offset_chebs(params);
    let [dx0, dx1, dx2, dx3] = approx_partials(params);
    let mut target = [0.0; 6];
    for i in 0..6 {
        target[i] = offset_cheb[i] * offset + dx1[i] * dth;
    }
    // compute inner products
    let mut int_target = Cheb::new(target.into()).integrate();
    let mut int_dx2 = Cheb::new(dx2.into()).integrate();
    let mut int_dx3 = Cheb::new(dx3.into()).integrate();
    int_target.zero_endpoints();
    int_dx2.zero_endpoints();
    int_dx3.zero_endpoints();
    let dx2_dx2 = (&int_dx2 * &int_dx2).def_integral();
    let dx2_dx3 = (&int_dx2 * &int_dx3).def_integral();
    let dx3_dx3 = (&int_dx3 * &int_dx3).def_integral();
    let dx2_target = (&int_dx2 * &int_target).def_integral();
    let dx3_target = (&int_dx3 * &int_target).def_integral();
    let det = dx2_dx2 * dx3_dx3 - dx2_dx3 * dx2_dx3;
    let x2 = (dx2_target * dx3_dx3 - dx2_dx3 * dx3_target) / det;
    let x3 = (dx3_target * dx2_dx2 - dx2_dx3 * dx2_target) / det;

    let offset_params = [params[0], params[1] - dth, params[2] + x2, params[3] + x3];
    println!("params: {params:?}");
    println!("offset: {offset_params:?}");
    let oc = crate::params_to_cubic(offset_params);
    CubicBez::new(
        p0,
        p0 + Vec2::new(oc.p1.x * ch.x - oc.p1.y * ch.y, oc.p1.y * ch.x + oc.p1.x * ch.y),
        p0 + Vec2::new(oc.p2.x * ch.x - oc.p2.y * ch.y, oc.p2.y * ch.x + oc.p2.x * ch.y),
        p1
    )
}
