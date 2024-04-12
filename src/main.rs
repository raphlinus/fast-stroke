use kurbo::{offset::CubicOffset, Affine, CubicBez, ParamCurve, ParamCurveArclen, ParamCurveDeriv, ParamCurveFit, QuadBez};

/// A cubic offset normalized to unit chord.
struct NormOffset {
    co_normed: CubicOffset,
    c_normed: CubicBez,
    // This is also normalized
    d: f64,
    th0: f64,
    th1: f64,
    area: f64,
    /// Arc length of underlying cubic
    c_arclen: f64,
    total_arclen: f64,
}
 
const ARCLEN_EPS: f64 = 1e-12;

impl NormOffset {
    fn new(c: CubicBez, d: f64) -> Self {
        let co = CubicOffset::new(c, d);
        let p0 = co.sample_pt_deriv(0.0).0;
        let p1 = co.sample_pt_deriv(1.0).0;
        let chord = p1 - p0;
        let u = chord.x / chord.hypot2();
        let v = chord.y / chord.hypot2();
        let rot = Affine::new([u, -v, v, u, 0., 0.]);
        let a = Affine::translate(-1. * (rot * p0).to_vec2()) * rot;
        println!("{:?} {:?}", a * p0, a * p1);
        let c_normed = a * c;
        let d = d * u.hypot(v);
        let co_normed = CubicOffset::new(c_normed, d);
        println!("{:?} {:?}", co_normed.sample_pt_deriv(0.0), co_normed.sample_pt_deriv(1.0));
        // TODO: we can get this analytically; verify against numeric solution
        let area = co_normed.moment_integrals(0.0 .. 1.0).0;
        let th0 = co_normed.sample_pt_deriv(0.0).1.angle();
        let th1 = -co_normed.sample_pt_deriv(1.0).1.angle();
        let c_arclen = c_normed.arclen(ARCLEN_EPS);
        let total_arclen = c_arclen + (th0 + th1) * d;
        println!("arclen = {}", total_arclen);
        NormOffset { co_normed, c_normed, d, th0, th1, area, c_arclen, total_arclen }
    }

    fn arclen(&self, t: f64) -> f64 {
        let th_t = self.co_normed.sample_pt_deriv(t).1.angle();
        self.c_normed.subsegment(0.0 .. t).arclen(ARCLEN_EPS) + (self.th0 - th_t) * self.d
    }
}

fn main() {
    let c = CubicBez::new((0.0, 0.0), (0.3, 0.1), (0.7, 0.2), (1.0, 0.0));
    let no = NormOffset::new(c, 0.1);
    println!("{} {} {}", no.area, no.th0, no.th1);
    println!("c_arclen = {}", no.c_arclen);
    println!("{:?} {:?}", no.arclen(0.5001) - no.arclen(0.5), no.co_normed.sample_pt_deriv(0.5).1.hypot());
}
