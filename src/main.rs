// This project is full of experiments.
#![allow(unused)]

mod approx;
mod balanced;
mod cheb;
mod norm_bez;
mod offset;

use balanced::BalancedBez;
use cheb::{cubic_to_chebs, Cheb, N_CHEB};
use kurbo::{
    common::solve_quadratic, offset::CubicOffset, Affine, CubicBez, ParamCurve, ParamCurveArclen,
    ParamCurveDeriv, ParamCurveFit, Point, Shape,
};
use rand::{thread_rng, Rng};

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

const ARCLEN_EPS: f64 = 1e-15;

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

    fn inv_arclen(&self, arclen: f64, total_arclen: f64) -> f64 {
        if arclen <= 0.0 {
            0.0
        } else if arclen >= total_arclen {
            1.0
        } else {
            let f = |t| self.arclen(t) - arclen;
            kurbo::common::solve_itp(
                f,
                0.0,
                1.0,
                ARCLEN_EPS,
                1,
                0.2,
                -arclen,
                total_arclen - arclen,
            )
        }
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
    fn approx_from_balance_ratio(&self, balance: f64) -> CubicBez {
        let s0 = self.th0.sin();
        let s1 = self.th1.sin();
        let s01 = (self.th0 + self.th1).sin();
        let c2 = -0.15 * balance * (1.0 - balance) * s01;
        let c1 = 0.3 * (balance * s0 + (1.0 - balance) * s1);
        let c0 = -self.area;
        let roots = solve_quadratic(c0, c1, c2);
        let mut d = 0.0;
        //println!("{roots:?}");
        for candidate in roots {
            if 0.0 <= candidate && candidate < 1.0 && (d == 0.0 || candidate < d) {
                d = candidate;
            }
        }
        let d0 = d * balance;
        let d1 = d * (1.0 - balance);
        //println!("d0 = {d0} d1 = {d1}");
        CubicBez::new(
            (0., 0.),
            (d0 * self.th0.cos(), d0 * s0),
            (1. - d1 * self.th1.cos(), d1 * s1),
            (1., 0.),
        )
    }

    fn cubic_fit(&self) -> CubicBez {
        let cubic_default: CubicBez =
            CubicBez::new((0., 0.), (1. / 3., 0.), (2. / 3., 0.), (1., 0.));
        match kurbo::fit_to_cubic(&self.co_normed, 0.0..1.0, 0.1) {
            Some((c, _err)) => c,
            _ => cubic_default,
        }
    }

    /// Approximate using sweet Bézier.
    fn approx_sweet(&self) -> CubicBez {
        let (s0, c0) = self.th0.sin_cos();
        let (s1, c1) = self.th1.sin_cos();
        let e0 = 2. / (3. + 3. * c0);
        let e1 = 2. / (3. + 3. * c1);
        let s01 = s0 * c1 + s1 * c0;
        let b0 = 2. * s0 - e1 * s01;
        let b1 = 2. * s1 - e0 * s01;
        // Maybe normalize here, or not necessary?
        let a0 = 0.15 * (2. * e0 * s0 + 2. * e1 * s1 - e0 * e1 * s01) - self.area;
        let a1 = 0.15 * (2. * b0 * s0 + 2. * b1 * s1 - (e0 * b1 + e1 * b0) * s01);
        let a2 = -0.15 * (b0 * b1 * s01);
        let roots = solve_quadratic(a0, a1, a2);
        //println!("{a2} u^2 + {a1} u + {a0}");
        //println!("roots: {roots:?}");
        let mut root = roots[0];
        if roots.len() > 1 && roots[1].abs() < root.abs() {
            root = roots[1];
        }
        let d0 = e0 + root * b0;
        let d1 = e1 + root * b1;
        //println!("d0 = {d0} d1 = {d1}");
        CubicBez::new(
            (0., 0.),
            (d0 * c0, d0 * s0),
            (1. - d1 * c1, d1 * s1),
            (1., 0.),
        )
    }

    /// Approximate finding area fit given starter d params.
    fn approx_from_ds(&self, e0: f64, e1: f64) -> CubicBez {
        let (s0, c0) = self.th0.sin_cos();
        let (s1, c1) = self.th1.sin_cos();
        let s01 = s0 * c1 + s1 * c0;
        let b0 = 2. * s0 - e1 * s01;
        let b1 = 2. * s1 - e0 * s01;
        // Maybe normalize here, or not necessary?
        let a0 = 0.15 * (2. * e0 * s0 + 2. * e1 * s1 - e0 * e1 * s01) - self.area;
        let a1 = 0.15 * (2. * b0 * s0 + 2. * b1 * s1 - (e0 * b1 + e1 * b0) * s01);
        let a2 = -0.15 * (b0 * b1 * s01);
        let roots = solve_quadratic(a0, a1, a2);
        let mut root = roots[0];
        if roots.len() > 1 && roots[1].abs() < root.abs() {
            root = roots[1];
        }
        let d0 = e0 + root * b0;
        let d1 = e1 + root * b1;
        CubicBez::new(
            (0., 0.),
            (d0 * c0, d0 * s0),
            (1. - d1 * c1, d1 * s1),
            (1., 0.),
        )
    }

    fn to_chebs(&self) -> [f64; N_CHEB] {
        let mut sums = [0.0; N_CHEB];
        let q = self.c_normed.deriv();
        let total_arclen = self.arclen(1.0);
        const N: usize = 16;
        let dx = 1.0 / N as f64;
        for k in 0..N {
            let x = (std::f64::consts::PI * dx * (k as f64 + 0.5)).cos();
            let t = self.inv_arclen((0.5 + 0.5 * x) * total_arclen, total_arclen);
            let th = q.eval(t).to_vec2().angle();
            for (sum, cheb) in sums.iter_mut().zip(cheb::chebs(x)) {
                *sum += th * cheb;
            }
        }
        for sum in sums.iter_mut() {
            *sum *= 2.0 * dx;
        }
        sums[0] *= 0.5;
        sums
    }
}

fn balance_ratio_of_cubic(c: CubicBez) -> f64 {
    let d0 = (c.p1 - c.p0).length();
    let d1 = (c.p3 - c.p2).length();
    d0 / (d0 + d1)
}

fn err_plot_main() {
    let bb = BalancedBez {
        th0: 0.2,
        th1: 0.3,
        bulge: 0.3,
        balance: 0.0,
    };
    let c = bb.to_cubic();
    const DIFF_EPS: f64 = 0.001;
    let no_p = NormOffset::new(c, DIFF_EPS);
    let bal_p = balance_ratio_of_cubic(no_p.cubic_fit());
    let no_m = NormOffset::new(c, -DIFF_EPS);
    let bal_m = balance_ratio_of_cubic(no_m.cubic_fit());
    let slope = (bal_p - bal_m) / (2. * DIFF_EPS);
    // Not a very good prediction; iterate
    let predict = 0.375 * (bb.th0.abs() - bb.th1.abs());
    eprintln!("slope = {slope}, predict {predict}");

    const HEIGHT: usize = 256;
    const WIDTH: usize = 256;
    println!("P3");
    println!("{WIDTH} {HEIGHT}");
    println!("255");
    for y in 0..HEIGHT {
        let d = y as f64 / (HEIGHT - 1) as f64 - 0.5;
        let no = NormOffset::new(c, d);
        let predict_bal = bb.balance + d * predict;
        for x in 0..WIDTH {
            let balance = x as f64 / (WIDTH - 1) as f64 - 0.5;
            let bb = BalancedBez {
                th0: no.th0,
                th1: no.th1,
                bulge: 0.0,
                balance,
            };
            let (d0, d1) = bb.get_ds();
            let c = no.approx_from_ds(d0, d1);
            let err = no.cubic_err(c);
            let z = -20. * err.ln() - 50.0;
            let g = z.clamp(0.0, 255.0) as u8;
            let g2 = if y == (HEIGHT - 1) / 2 {
                z.clamp(0.0, 128.0) as u8 + 127
            } else {
                g
            };
            let r2 = if (balance - predict_bal).abs() < 0.005 {
                z.clamp(0.0, 128.0) as u8 + 127
            } else {
                g
            };
            println!("{r2} {g2} {g}");
        }
    }
}

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
    println!("{:?}", no.approx_from_balance_ratio(0.46725));
}

#[derive(Debug)]
struct EsErrorMetric {
    symm: f64,
    asymm: f64,
    aerr: f64,
    dist: f64,
    err: f64,
}

/// Compute ES error metric and its components.
///
/// Input cubic Bézier must be normalized to unit chord.
fn es_err_metric(c: CubicBez) -> EsErrorMetric {
    let th0 = c.p1.y.atan2(c.p1.x);
    let th1 = c.p2.y.atan2(1.0 - c.p2.x);
    let d0 = c.p1.y.hypot(c.p1.x);
    let d1 = c.p2.y.hypot(1.0 - c.p2.x);
    let e0 = (2. / 3.) / (1.0 + th0.cos());
    let e1 = (2. / 3.) / (1.0 + th1.cos());
    let s0 = th0.sin();
    let s1 = th1.sin();
    let s01 = (s0 + s1).sin();
    let amin = 0.15 * (2. * e0 * s0 + 2. * e1 * s1 - e0 * e1 * s01);
    let a = 0.15 * (2. * d0 * s0 + 2. * d1 * s1 - d0 * d1 * s01);
    let aerr = (a - amin).abs();
    let symm = (th0 + th1).abs();
    let asymm = (th0 - th1).abs();
    let dist = (d0 - e0).hypot(d1 - e1);
    let ctr = 3.7e-6 * symm.powi(5) + 6e-3 * asymm * symm.powi(2);
    let halo_symm = 5e-3 * symm * dist;
    let halo_asymm = 7e-2 * asymm * dist;
    let err = 1.25 * ctr + 1.55 * aerr + halo_symm + halo_asymm;
    EsErrorMetric {
        symm,
        asymm,
        aerr,
        dist,
        err,
    }
}

fn err_metric_main() {
    let mut rng = thread_rng();
    let th0 = rng.gen_range(0.0f64..0.5);
    let th1 = rng.gen_range(0.0f64..0.5);
    let d0 = rng.gen_range(0.0..0.6);
    let d1 = rng.gen_range(0.0..0.6);
    let p2 = Point::new(1.0 - d1 * th1.cos(), d1 * th1.sin());
    let p1 = Point::new(d0 * th0.cos(), d0 * th0.sin());
    let c = CubicBez::new(Point::ORIGIN, p1, p2, Point::new(1.0, 0.0));
    println!("{:.5?}", es_err_metric(c));
    let c_subdiv = c.subsegment(0.45..0.55);
    println!("{}", c.to_path(0.1).to_svg());
    println!("{}", c_subdiv.to_path(0.1).to_svg());
    println!(
        "{}",
        norm_bez::normalize_bez(c_subdiv).0.to_path(0.1).to_svg()
    );
    println!("{:.5?}", es_err_metric(norm_bez::normalize_bez(c_subdiv).0));
}

fn arclen_foo() {
    let mut rng = thread_rng();
    let th0 = rng.gen_range(0.0f64..0.5);
    let th1 = rng.gen_range(0.0f64..0.5);
    let d0 = rng.gen_range(0.0..0.6);
    let d1 = rng.gen_range(0.0..0.6);
    let p2 = Point::new(1.0 - d1 * th1.cos(), d1 * th1.sin());
    let p1 = Point::new(d0 * th0.cos(), d0 * th0.sin());
    let c = CubicBez::new(Point::ORIGIN, p1, p2, Point::new(1.0, 0.0));
    let normed = NormOffset::new(c, 0.1);
    let total_arclen = normed.arclen(1.0);
    for i in 0..=10 {
        let arclen = total_arclen * i as f64 / 10.0;
        println!("{arclen} {}", normed.inv_arclen(arclen, total_arclen));
    }
    println!("{:?}", normed.to_chebs());
    let normed_minus = NormOffset::new(c, -0.1);
    println!("{:?}", normed_minus.to_chebs());
}

fn cheb_foo() {
    let p0 = Cheb::new(vec![1., 0.]);
    let p1 = Cheb::new(vec![0., 1.]);
    let p2 = Cheb::new(vec![0., 0., 1.]);
    let p3 = Cheb::new(vec![0., 0., 0., 1.]);
    let p4 = Cheb::new(vec![0., 0., 0., 0., 1.]);
    println!("{:?}", &p0 * &p0);
    println!("{:?}", &p0 * &p1);
    println!("{:?}", &p1 * &p1);
    println!("{:?}", &p1 * &p4);
    println!("{:?}", &p4 * &p1);
    println!("{:?}", &p4 - &p1);
    println!("{:?}", &p1 - &p4);
    println!("{:?}", &p4 + &p1);
}

/// Convert parameters to cubic Bézier.
///
/// There are a bunch of variations to try, but this one seems to be
/// pretty good. The first two parameters are sum and difference of
/// endpoint angles, and the second two are sum and difference of control
/// lengths as a multiple of the two-parabola lengths.
fn params_to_cubic(params: [f64; 4]) -> CubicBez {
    let [a, b, c, d] = params;
    let (th0, th1) = (a + b, a - b);
    let (s0, c0) = th0.sin_cos();
    let (s1, c1) = th1.sin_cos();
    let (x1, x2) = (c + d, c - d);
    let d0 = (1. + x1) * (2. / 3.) / (1. + c0);
    let d1 = (1. + x2) * (2. / 3.) / (1. + c1);
    let y1 = s0 * d0;
    let x1 = c0 * d0 - 1. / 3.;
    let y2 = s1 * d1;
    let x2 = c1 * d1 - 1. / 3.;
    CubicBez::new((0., 0.), (1. / 3. + x1, y1), (2. / 3. - x2, y2), (1., 0.))
}

fn cubic_cheb() {
    for order in 1..=5 {
        println!("    // order {order}");
        let f = |params: [f64; 4]| {
            let c = params_to_cubic(params);
            cheb::cubic_to_chebs(&c)
        };
        cheb::cheb_deriv_order(order, f, cheb::ReportStyle::Derivs);
    }
}

fn offset_partial_cheb(c: CubicBez) -> [f64; 6] {
    let normed_plus = NormOffset::new(c, OFFSET_DERIV_EPS);
    let plus_chebs = normed_plus.to_chebs();
    let normed_minus = NormOffset::new(c, -OFFSET_DERIV_EPS);
    let minus_chebs = normed_minus.to_chebs();
    let mut chebs = [0.0; 6];
    for i in 0..6 {
        chebs[i] = (plus_chebs[i] - minus_chebs[i]) * (0.5 / OFFSET_DERIV_EPS)
    }
    chebs
}

const OFFSET_DERIV_EPS: f64 = 3e-2;
fn offset_cheb() {
    for order in 1..=5 {
        println!("    // order {order}");
        let f = |params: [f64; 4]| offset_partial_cheb(params_to_cubic(params));
        cheb::cheb_deriv_order(order, f, cheb::ReportStyle::Polynom);
    }
}

fn random_cubic() -> CubicBez {
    let mut rng = thread_rng();
    let th0 = rng.gen_range(0.0f64..0.5);
    let th1 = rng.gen_range(-0.5f64..0.5);
    let d0 = rng.gen_range(0.1..0.6);
    let d1 = rng.gen_range(0.1..0.6);
    let p2 = Point::new(1.0 - d1 * th1.cos(), d1 * th1.sin());
    let p1 = Point::new(d0 * th0.cos(), d0 * th0.sin());
    CubicBez::new(Point::ORIGIN, p1, p2, Point::new(1.0, 0.0))
}

fn param_scaling() {
    let c = random_cubic();
    println!("{:?}", norm_bez::normed_bez_to_params(c));
    let c_subdiv = c.subsegment(0.45..0.55);
    let norm_subdiv = norm_bez::normalize_bez(c_subdiv).0;
    println!("{:?}", norm_bez::normed_bez_to_params(norm_subdiv));
}

fn check_offset() {
    let c = random_cubic();
    let params = norm_bez::normed_bez_to_params(c);
    println!("params: {:.4?}", params);
    println!("offset_partial_cheb: {:.4?}", offset_partial_cheb(c));
    println!(
        "          by approx: {:.4?}",
        approx::approx_offset_chebs(params)
    );
    for p in approx::approx_partials(params) {
        println!("{p:.4?}");
    }
    let c_offset = offset::approx_offset(c, 0.1);
    println!("{}", c.to_path(1e-9).to_svg());
    println!("{}", c_offset.to_path(1e-9).to_svg());
}

fn main() {
    // err_metric_main();
    // arclen_foo();
    // cheb_foo();
    //cubic_cheb();
    //offset_cheb();
    //param_scaling();
    check_offset();
}
