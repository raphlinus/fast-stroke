//! Stuff for manipulating Chebyshev polynomials

use std::f64::consts::PI;

use kurbo::{CubicBez, ParamCurve, ParamCurveArclen, ParamCurveDeriv};

pub const N_CHEB: usize = 6;
const DERIV_EPS: f64 = 1e-2;

fn deriv_x(i: usize, order: usize) -> f64 {
    (i as f64 - 0.5 * order as f64) * DERIV_EPS
}

fn binom(i: usize, n: usize) -> i32 {
    let a: &[i32] = match n {
        0 => &[1],
        1 => &[-1, 1],
        2 => &[1, -2, 1],
        3 => &[-1, 3, -3, 1],
        4 => &[1, -4, 6, -4, 1],
        5 => &[-1, 5, -10, 10, -5, 1],
        _ => unimplemented!(),
    };
    a[i]
}

#[derive(Clone, Copy)]
pub enum ReportStyle {
    StupidDebug,
    Polynom,
    Derivs,
}

fn fact(x: usize) -> f64 {
    [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0][x]
}

fn term_expr(orders: &[usize]) -> String {
    use std::fmt::Write;
    let mut term = String::new();
    for (i, order) in orders.iter().enumerate() {
        let add = if term.is_empty() { "" } else { " * " };
        if *order == 1 {
            _ = write!(&mut term, "{add}x{i}");
        } else if *order == 2 {
            _ = write!(&mut term, "{add}x{i} * x{i}");
        } else if *order > 0 {
            _ = write!(&mut term, "{add}x{i}.powi({order})");
        }
    }
    if term.is_empty() {
        _ = write!(&mut term, "1.0");
    }
    term
}

fn cheb_deriv_1(orders: [usize; 4], f: impl Fn([f64; 4]) -> [f64; N_CHEB], report: ReportStyle) {
    let mut sums = [0.0; N_CHEB];
    for i0 in 0..=orders[0] {
        let x0 = deriv_x(i0, orders[0]);
        for i1 in 0..=orders[1] {
            let x1 = deriv_x(i1, orders[1]);
            for i2 in 0..=orders[2] {
                let x2 = deriv_x(i2, orders[2]);
                for i3 in 0..=orders[3] {
                    let x3 = deriv_x(i3, orders[3]);
                    let weight = binom(i0, orders[0])
                        * binom(i1, orders[1])
                        * binom(i2, orders[2])
                        * binom(i3, orders[3]);
                    let ys = f([x0, x1, x2, x3]);
                    for (sum, y) in sums.iter_mut().zip(ys) {
                        *sum += weight as f64 * y;
                    }
                }
            }
        }
    }
    let fact_scale: f64 = orders.iter().map(|order| fact(*order)).product();
    let scale = (1. / DERIV_EPS).powi(orders.iter().sum::<usize>() as i32) / fact_scale;
    for sum in &mut sums {
        *sum *= scale;
    }
    const THRESH: f64 = 1e-3;
    match report {
        ReportStyle::StupidDebug => println!("{orders:?}: {sums:.3?}"),
        ReportStyle::Polynom => {
            if sums[1..].iter().any(|x| x.abs() > THRESH) {
                println!("    let term = {};", term_expr(&orders));
                for (i, sum) in sums.iter().enumerate() {
                    if i > 0 && sum.abs() > THRESH {
                        println!("    c{i} += {} * term;", sum);
                    }
                }
            }
        }
        ReportStyle::Derivs => {
            if sums[1..].iter().any(|x| x.abs() > THRESH) {
                for (j, order) in orders.iter().enumerate() {
                    if *order > 0 {
                        let mut d_orders = orders;
                        d_orders[j] -= 1;
                        println!("    let term_{j} = {};", term_expr(&d_orders));
                        for (i, sum) in sums.iter().enumerate() {
                            if i > 0 && sum.abs() > THRESH {
                                println!("    dc{i}_dx{j} += {} * term_{j};", *order as f64 * sum);
                            }
                        }
                    }
                }
            }
        }
    }
}

pub fn cheb_deriv_order(order: usize, f: impl Fn([f64; 4]) -> [f64; N_CHEB], report: ReportStyle) {
    for o0 in 0..=order {
        let rem = order - o0;
        for o1 in 0..=rem {
            let rem = rem - o1;
            for o2 in 0..=rem {
                let rem = rem - o2;
                let o3 = rem;
                cheb_deriv_1([o0, o1, o2, o3], &f, report);
            }
        }
    }
}

pub fn chebs(x: f64) -> [f64; 6] {
    let x2 = x * x;
    [
        1.,
        x,
        2. * x2 - 1.,
        (4. * x2 - 3.) * x,
        8. * x2 * x2 - 8. * x2 + 1.,
        (16. * x2 * x2 - 20. * x2 + 5.) * x,
    ]
}

pub fn chebs_generic<const N: usize>(x: f64) -> [f64; N] {
    let mut result = [0.0; N];
    let mut a = 1.0;
    let mut b = x;
    for y in &mut result {
        *y = a;
        let bnew = 2.0 * x * b - a;
        (a, b) = (b, bnew)
    }
    result
}

pub fn cubic_to_chebs(c: &CubicBez) -> [f64; N_CHEB] {
    let mut sums = [0.0; N_CHEB];
    let q = c.deriv();
    const ARCLEN_ACCURACY: f64 = 1e-15;
    let arclen = c.arclen(ARCLEN_ACCURACY);
    const N: usize = 16;
    let dx = 1.0 / N as f64;
    for k in 0..N {
        let x = (PI * dx * (k as f64 + 0.5)).cos();
        let s = 0.5 + 0.5 * x;
        let t = c.inv_arclen(s * arclen, ARCLEN_ACCURACY);
        let th = q.eval(t).to_vec2().angle();
        for (sum, cheb) in sums.iter_mut().zip(chebs(x)) {
            *sum += th * cheb;
        }
    }
    for sum in sums.iter_mut() {
        *sum *= 2.0 * dx;
    }
    sums[0] *= 0.5;
    sums
}

// Could be array.
#[derive(Clone, Debug)]
pub struct Cheb(Vec<f64>);

impl Cheb {
    pub fn new(v: Vec<f64>) -> Self {
        Cheb(v)
    }

    fn eval(&self, x: f64) -> f64 {
        let mut a = 1.0;
        let mut b = x;
        let mut sum = a * self.0[0] + b * self.0[1];
        for c in &self.0[2..] {
            let bnew = 2.0 * x * b - a;
            sum += *c * bnew;
            (a, b) = (b, bnew)
        }
        sum
    }

    fn get(&self, i: usize) -> f64 {
        *self.0.get(i).unwrap_or(&0.0)
    }

    fn zget(&self, i: isize) -> f64 {
        if i == 0 {
            self.0[0]
        } else {
            match self.0.get(i.abs() as usize) {
                Some(x) => 0.5 * x,
                None => 0.0,
            }
        }
    }

    // This could be arbitrary multiply easily.
    fn square(&self) -> Self {
        let mut v = Vec::with_capacity(self.0.len() * 2 - 1);
        let mut sum = self.0[0].powi(2);
        for c in &self.0[1..] {
            sum += 0.5 * c * c;
        }
        v.push(sum);
        let n = self.0.len() as isize;
        for i in 1..2 * n - 1 {
            sum = 0.0;
            for j in i - n + 1..n {
                sum += self.zget(i - j) * self.zget(j);
            }
            v.push(2. * sum);
        }
        Cheb(v)
    }

    fn l2_norm(&self) -> f64 {
        let s = self.square();
        let mut sum = s.0[0];
        for i in (2..s.0.len()).step_by(2) {
            sum -= s.0[i] / ((i * i - 1) as f64);
        }
        sum.sqrt()
    }

    fn zero_endpoints(&mut self) {
        let mut even = 0.0;
        for i in (2..self.0.len()).step_by(2) {
            even += self.0[i];
        }
        let mut odd = 0.0;
        for i in (3..self.0.len()).step_by(2) {
            odd += self.0[i];
        }
        self.0[0] = -even;
        self.0[1] = -odd;
    }

    pub fn integrate(&self) -> Cheb {
        let n = self.0.len();
        let mut v = Vec::with_capacity(n + 1);
        v.push(0.25 * self.get(1) - 0.375 * self.get(3));
        v.push(self.0[0] - 0.5 * self.get(2));
        for i in 2..n + 1 {
            v.push((self.0[i - 1] - self.get(i + 1)) * 0.5 / i as f64);
        }
        Cheb(v)
    }
}

impl std::ops::Mul<f64> for &Cheb {
    type Output = Cheb;

    fn mul(self, rhs: f64) -> Cheb {
        Cheb(self.0.iter().map(|x| x * rhs).collect())
    }
}

impl std::ops::Mul for &Cheb {
    type Output = Cheb;

    fn mul(self, rhs: Self) -> Cheb {
        let mut v = Vec::with_capacity(self.0.len() + rhs.0.len() - 1);
        let mut sum = self.0[0] * rhs.0[0];
        for (a, b) in self.0[1..].iter().zip(&rhs.0[1..]) {
            sum += 0.5 * a * b;
        }
        let self_n = self.0.len() as isize;
        let rhs_n = rhs.0.len() as isize;
        v.push(sum);
        for i in 1..self_n + rhs_n - 1 {
            sum = 0.0;
            for j in i - self_n + 1..rhs_n {
                sum += self.zget(i - j) * rhs.zget(j);
            }
            v.push(2. * sum);
        }
        Cheb(v)
    }
}

impl std::ops::Add for &Cheb {
    type Output = Cheb;

    fn add(self, rhs: Self) -> Cheb {
        let mut v = Vec::with_capacity(self.0.len().max(rhs.0.len()));
        for (a, b) in self.0.iter().zip(&rhs.0) {
            v.push(a + b);
        }
        if self.0.len() > rhs.0.len() {
            v.extend(&self.0[rhs.0.len()..]);
        }
        if self.0.len() < rhs.0.len() {
            v.extend(&rhs.0[self.0.len()..]);
        }
        Cheb(v)
    }
}

impl std::ops::Sub for &Cheb {
    type Output = Cheb;

    fn sub(self, rhs: Self) -> Cheb {
        let mut v = Vec::with_capacity(self.0.len().max(rhs.0.len()));
        for (a, b) in self.0.iter().zip(&rhs.0) {
            v.push(a - b);
        }
        if self.0.len() > rhs.0.len() {
            v.extend(&self.0[rhs.0.len()..]);
        }
        if self.0.len() < rhs.0.len() {
            v.extend(rhs.0[self.0.len()..].iter().map(|x| -x));
        }
        Cheb(v)
    }
}
