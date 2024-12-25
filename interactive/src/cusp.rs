//! Cusp analysis and mitigation.

use std::ops::Range;

use arrayvec::ArrayVec;
use kurbo::{common::solve_cubic, ParamCurve, QuadBez, Vec2};

#[derive(Debug)]
pub struct CuspAnalysis {
    vmin: f64,
    vmax: f64,
    minima: ArrayVec<VelMinimum, 2>,
}

#[derive(Clone, Copy, Debug)]
struct VelMinimum {
    t: f64,
    v: f64,
    utan: Vec2,
}

impl CuspAnalysis {
    /// Do cusp analysis based on hodograph of cubic Bezier.
    pub fn new(q: QuadBez) -> Self {
        let d0 = q.p1 - q.p0;
        let d1 = (q.p2 - q.p1) - d0;
        let d = q.p0.to_vec2();
        let c0 = d.dot(d0);
        let c1 = 2.0 * d0.hypot2() + d.dot(d1);
        let c2 = 3.0 * d1.dot(d0);
        let c3 = d1.hypot2();
        let mut minima = ArrayVec::new();
        let v0 = q.p0.to_vec2().hypot();
        let v2 = q.p2.to_vec2().hypot();
        let mut vmin = v0.min(v2);
        let mut vmax = v0.max(v2);
        for t in solve_cubic(c0, c1, c2, c3) {
            if (0.0..1.0).contains(&t) {
                let tan = q.eval(t).to_vec2();
                let v = tan.hypot();
                // TODO: handle near-zero case
                let utan = tan / v;
                if c1 + 2.0 * c2 * t + 3.0 * c3 * t * t > 0.0 {
                    minima.push(VelMinimum { t, v, utan });
                    vmin = vmin.min(v);
                } else {
                    vmax = vmax.max(v);
                }
            }
        }
        CuspAnalysis { vmin, vmax, minima }
    }

    pub fn get_cusp(&self, range: Range<f64>) -> Option<(f64, Vec2)> {
        let mut result = None;
        for m in &self.minima {
            // should be is_none_or, but msrv 1.82 and not sure we want to bump
            if range.start < m.t
                && range.end > m.t
                && result.map(|r: VelMinimum| r.v > m.v).unwrap_or(true)
            {
                result = Some(*m);
            }
        }
        result.map(|m| (m.t, m.utan))
    }
}
