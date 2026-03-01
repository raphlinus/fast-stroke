//! An algorithm for computing the minumum Frechet error approximation.

use kurbo::{CubicBez, ParamCurve, ParamCurveDeriv};

use crate::perturb::{turn, CurveOffset, OffsetSolutionLse};

pub struct MinFrechet {
    /// The t value on the true offset corresponding to 0.5 on the approximation
    t: f64,
    /// An additional perturbation in the unit normal direction
    delta: f64,
}

#[derive(Debug)]
pub enum IterResult {
    // TODO: might want to report error in Done variant
    Done,
    Fail,
    Step,
}

impl MinFrechet {
    pub fn new(t: f64) -> Self {
        Self { t, delta: 0.0 }
    }

    // This gets (a, b) values as in perturb (scaled by endpoint derivs)
    pub fn get_a_b(&self, c: CubicBez, d: f64) -> (f64, f64) {
        let q = c.deriv();
        let b01 = q.p0.to_vec2();
        let b23 = q.p2.to_vec2();
        let n = turn(q.eval(self.t).to_vec2().normalize());
        let target = c.eval(self.t) + (d + self.delta) * n;
        let ca = 0.375 * b01;
        let cb = 0.375 * b23;
        let n0 = turn(b01).normalize();
        let n1 = turn(b23).normalize();
        let cc = 0.5 * (n0 + n1);
        let z = target - (c.eval(0.5) + d * cc);
        let det = d * ca.cross(cb);
        let a = z.cross(cb) / det;
        let b = ca.cross(z) / det;
        (a, b)
    }

    pub fn iterate(&mut self, co: &CurveOffset, d: f64) -> IterResult {
        let (a, b) = self.get_a_b(co.c, d);
        let soln = OffsetSolutionLse::from_a_b(a, b, d);
        let extrema = soln.find_error_extrema(co, self.t);
        if extrema.len() != 3 {
            // TODO: can probably handle 4 extrema...
            return IterResult::Fail;
        }
        //let err = extrema.iter().fold(0.0, |m, e| e.error.abs().max(m));
        let e01 = extrema[0].error + extrema[1].error;
        let e12 = extrema[1].error + extrema[2].error;
        const THRESH: f64 = 1e-12;
        if e01.abs().max(e12.abs()) < THRESH {
            return IterResult::Done;
        }
        let dp_dt = co.q.eval(self.t).to_vec2();
        let dp_dv = turn(dp_dt).normalize();
        let ca = 0.375 * co.q.p0.to_vec2();
        let cb = 0.375 * co.q.p2.to_vec2();
        // The actual derivatives here may have a factor of d, but that cancels
        // out in the second determinant.
        let idet = 1.0 / ca.cross(cb);
        let da_dt = dp_dt.cross(cb) * idet;
        let db_dt = ca.cross(dp_dt) * idet;
        let da_dv = dp_dv.cross(cb) * idet;
        let db_dv = ca.cross(dp_dv) * idet;
        let de_dt_v: [_; 3] = core::array::from_fn(|i| {
            let extremum = &extrema[i];
            let unorm = turn(co.q.eval(extremum.t_offset).to_vec2()).normalize();
            let ta = extremum.t_approx;
            let mta = 1. - ta;
            let de_da = 3. * mta * ta * mta * co.q.p0.to_vec2().dot(unorm);
            let de_db = 3. * mta * ta * ta * co.q.p2.to_vec2().dot(unorm);
            let de_dt = de_da * da_dt + de_db * db_dt;
            let de_dv = de_da * da_dv + de_db * db_dv;
            (de_dt, de_dv)
        });
        let de01_dt = de_dt_v[0].0 + de_dt_v[1].0;
        let de01_dv = de_dt_v[0].1 + de_dt_v[1].1;
        let de12_dt = de_dt_v[1].0 + de_dt_v[2].0;
        let de12_dv = de_dt_v[1].1 + de_dt_v[2].1;
        let e01 = extrema[0].error + extrema[1].error;
        let e12 = extrema[1].error + extrema[2].error;
        let idet2 = 1.0 / (de01_dt * de12_dv - de12_dt * de01_dv);
        let dt = (e01 * de12_dv - e12 * de01_dv) * idet2;
        let dz = (de01_dt * e12 - de12_dt * e01) * idet2;
        self.t += dt;
        self.delta += dz;
        IterResult::Step
    }
}
