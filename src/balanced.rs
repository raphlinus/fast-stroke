use kurbo::CubicBez;

#[derive(Debug)]
pub struct BalancedBez {
    pub th0: f64,
    pub th1: f64,
    // these are normalized, but un-normalized could be interesting
    pub bulge: f64,
    pub balance: f64,
}

impl BalancedBez {
    pub fn to_cubic(&self) -> CubicBez {
        let (s0, c0) = self.th0.sin_cos();
        let (s1, c1) = self.th1.sin_cos();
        let e0 = 2. / (3. + 3. * c0);
        let e1 = 2. / (3. + 3. * c1);
        let s01 = s0 * c1 + s1 * c0;
        let b0 = 2. * s0 - e1 * s01;
        let b1 = 2. * s1 - e0 * s01;
        let norm = 1. / b0.hypot(b1);
        let b0 = b0 * norm;
        let b1 = b1 * norm;
        let d0 = e0 + self.bulge * b0 + self.balance * b1;
        let d1 = e1 + self.bulge * b1 - self.balance * b0;
        //println!("d0 = {d0} d1 = {d1}");
        CubicBez::new(
            (0., 0.),
            (d0 * c0, d0 * s0),
            (1. - d1 * c1, d1 * s1),
            (1., 0.),
        )
    }

    pub fn from_cubic(c: CubicBez) -> Self {
        let d01 = c.p1 - c.p0;
        let d23 = c.p3 - c.p2;
        let th0 = d01.atan2();
        let th1 = -d23.atan2();
        let l01 = d01.hypot();
        let l23 = d23.hypot();
        let (s0, c0) = (d01.y / l01, d01.x / l01);
        let (s1, c1) = (-d23.y / l23, d23.x / l23);
        let e0 = 2. / (3. + 3. * c0);
        let e1 = 2. / (3. + 3. * c1);
        let s01 = s0 * c1 + s1 * c0;
        let b0 = 2. * s0 - e1 * s01;
        let b1 = 2. * s1 - e0 * s01;
        let norm = 1. / b0.hypot(b1);
        let b0 = b0 * norm;
        let b1 = b1 * norm;
        let bulge = (l01 - e0) * b0 + (l23 - e1) * b1;
        let balance = (l01 - e0) * b1 - (l23 - e1) * b0;
        BalancedBez {
            th0,
            th1,
            bulge,
            balance,
        }
    }

    pub fn get_ds(&self) -> (f64, f64) {
        let (s0, c0) = self.th0.sin_cos();
        let (s1, c1) = self.th1.sin_cos();
        let e0 = 2. / (3. + 3. * c0);
        let e1 = 2. / (3. + 3. * c1);
        let s01 = s0 * c1 + s1 * c0;
        let b0 = 2. * s0 - e1 * s01;
        let b1 = 2. * s1 - e0 * s01;
        let norm = 1. / b0.hypot(b1);
        let b0 = b0 * norm;
        let b1 = b1 * norm;
        let d0 = e0 + self.bulge * b0 + self.balance * b1;
        let d1 = e1 + self.bulge * b1 - self.balance * b0;
        (d0, d1)
    }
}
