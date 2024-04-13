use kurbo::{
    common::solve_quadratic, offset::CubicOffset, Affine, CubicBez, ParamCurve, ParamCurveArclen,
    ParamCurveDeriv, ParamCurveFit, QuadBez,
};

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
        let c_normed = a * c;
        let d = d * u.hypot(v);
        let co_normed = CubicOffset::new(c_normed, d);
        // TODO: we can get this analytically; verify against numeric solution
        let area = co_normed.moment_integrals(0.0..1.0).0;
        let th0 = co_normed.sample_pt_deriv(0.0).1.angle();
        let th1 = -co_normed.sample_pt_deriv(1.0).1.angle();
        let c_arclen = c_normed.arclen(ARCLEN_EPS);
        let total_arclen = c_arclen + (th0 + th1) * d;
        NormOffset {
            co_normed,
            c_normed,
            d,
            th0,
            th1,
            area,
            c_arclen,
            total_arclen,
        }
    }

    fn arclen(&self, t: f64) -> f64 {
        let th_t = self.co_normed.sample_pt_deriv(t).1.angle();
        self.c_normed.subsegment(0.0..t).arclen(ARCLEN_EPS) + (self.th0 - th_t) * self.d
    }

    /// Estimate of Fréchet distance.
    ///
    /// This is computed as L-infty norm of distance parametrized by arc length. It might also
    /// be worth considering L2 norms.
    fn cubic_err(&self, c: CubicBez) -> f64 {
        let c_arclen = c.arclen(ARCLEN_EPS);
        const N: usize = 20;
        let mut max_err2 = 0.0;
        for i in 1..N {
            let t = i as f64 * (1.0 / N as f64);
            let p_co = self.co_normed.sample_pt_deriv(t).0;
            let s = self.arclen(t) / self.total_arclen;
            let cubic_t = c.inv_arclen(c_arclen * s, ARCLEN_EPS);
            let p_c = c.eval(cubic_t);
            max_err2 = (p_co - p_c).hypot2().max(max_err2);
        }
        max_err2.sqrt()
    }

    /// Create a cubic approximation given a balance parameter.
    ///
    /// In this case, balance is d0 / (d0 + d1), so meaningfully ranges from 0
    /// to 1. The returned Bézier is C1 continuous.
    fn approx_from_balance(&self, balance: f64) -> CubicBez {
        let s0 = self.th0.sin();
        let s1 = self.th1.sin();
        let s01 = (self.th0 + self.th1).sin();
        let c2 = -0.15 * balance * (1.0 - balance) * s01;
        let c1 = 0.3 * (balance * s0 + (1.0 - balance) * s1);
        let c0 = -self.area;
        let roots = solve_quadratic(c0, c1, c2);
        // TODO: probably will need to pick best root
        let d = roots.get(0).unwrap_or(&0.0);
        let d0 = d * balance;
        let d1 = d * (1.0 - balance);
        CubicBez::new(
            (0., 0.),
            (d0 * self.th0.cos(), d0 * s0),
            (1. - d1 * self.th1.cos(), d1 * s1),
            (1., 0.),
        )
    }
}

fn err_plot_main() {
    let c = CubicBez::new((0.0, 0.0), (0.3, 0.2), (0.6, 0.1), (1.0, 0.0));
    const HEIGHT: usize = 256;
    const WIDTH: usize = 256;
    println!("P3");
    println!("{WIDTH} {HEIGHT}");
    println!("255");
    for y in 0..HEIGHT {
        let d = y as f64 / (HEIGHT - 1) as f64 - 0.5;
        let no = NormOffset::new(c, d);
        for x in 0..WIDTH {
            let balance = x as f64 / (WIDTH - 1) as f64;
            let c = no.approx_from_balance(balance);
            let err = no.cubic_err(c);
            let z = -30. * err.ln() - 100.0;
            let g = z.clamp(0.0, 255.0) as u8;
            let g2 = if y == (HEIGHT - 1) / 2 { z.clamp(0.0, 128.0) as u8 + 127} else { g };
            println!("{g} {g2} {g}");
        }
    }
}

#[allow(unused)]
fn junk_main() {
    let c = CubicBez::new((0.0, 0.0), (0.3, 0.1), (0.7, 0.2), (1.0, 0.0));
    let no = NormOffset::new(c, -0.2);
    println!("{} {} {}", no.area, no.th0, no.th1);
    println!("c_arclen = {}", no.c_arclen);
    println!(
        "{:?} {:?}",
        no.arclen(0.5001) - no.arclen(0.5),
        no.co_normed.sample_pt_deriv(0.5).1.hypot()
    );
    println!("err = {}", no.cubic_err(c));
    println!("{:?}", no.approx_from_balance(0.46725));
}

fn main() {
    err_plot_main();
}