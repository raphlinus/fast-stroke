//! Stuff for the cubic to Euler spiral error metric

use kurbo::{CubicBez, ParamCurve, ParamCurveArclen};

use crate::{cheb::Cheb, norm_bez, ARCLEN_EPS};

/// Measure distance between two cubics with L2 norm of arc length parametrized distance.
fn measure_cubics(c0: CubicBez, c1: CubicBez) -> f64 {
    const N: usize = 32;
    let dt = 1.0 / N as f64;
    let arclen0 = c0.arclen(ARCLEN_EPS);
    let arclen1 = c1.arclen(ARCLEN_EPS);
    let mut err2 = 0.0;
    for i in 0..N {
        let t = (i as f64 + 0.5) * dt;
        let t0 = c0.inv_arclen(t * arclen0, ARCLEN_EPS);
        let t1 = c1.inv_arclen(t * arclen1, ARCLEN_EPS);
        err2 += (c0.eval(t0) - c1.eval(t1)).length_squared();
    }
    (err2 * dt).sqrt()
}

/// Make double-parabola approximation.
fn approx_euler_of(c: CubicBez) -> CubicBez {
    let d01 = c.p1 - c.p0;
    let d23 = c.p3 - c.p2;
    let th0 = d01.atan2();
    let th1 = -d23.atan2();
    let (s0, c0) = th0.sin_cos();
    let (s1, c1) = th1.sin_cos();
    let d0 = 2. / (3. + 3. * c0);
    let d1 = 2. / (3. + 3. * c1);
    CubicBez::new(
        (0., 0.),
        (d0 * c0, d0 * s0),
        (1. - d1 * c1, d1 * s1),
        (1., 0.),
    )
}

/// Estimate error using error metric.
///
/// This is the estimated error to the double parabola approximation; Euler
/// error is reasonably computed with triangle inequality.
fn est_error(c: CubicBez) -> f64 {
    0.0
}

fn chebs_of(c: CubicBez) -> Cheb {
    let params = norm_bez::normed_bez_to_params(c);
    let [x0, x1, x2, x3] = params;
    let mut c1 = 0.0;
    let mut c2 = 0.0;
    let mut c3 = 0.0;
    let mut c4 = 0.0;
    let mut c5 = 0.0;
    // order 1
    let term = x1;
    c2 += 0.7500006249988824 * term;
    let term = x0;
    c1 += -1.0000000000012204 * term;
    // order 2
    let term = x1 * x3;
    c1 += -0.7500052246339631 * term;
    c3 += 0.7499879393079562 * term;
    let term = x1 * x2;
    c2 += 0.37500029579497046 * term;
    c4 += -0.4687454687057695 * term;
    let term = x0 * x2;
    c1 += -0.4999988542213514 * term;
    c3 += 0.4999942056832605 * term;
    // order 3
    let term = x1 * x3 * x3;
    c2 += -0.5625154189965037 * term;
    c4 += 0.7030798665500154 * term;
    let term = x1 * x2 * x3;
    c1 += 0.2812635799462598 * term;
    c3 += 0.703135262521571 * term;
    c5 += -0.9843254438978398 * term;
    let term = x1 * x2 * x2;
    c2 += -0.1406290655826467 * term;
    c4 += -0.14062653106762263 * term;
    let term = x1.powi(3);
    c2 += 0.02499955303533752 * term;
    c4 += 0.014061679637493382 * term;
    let term = x0 * x3 * x3;
    c1 += -0.1250018620575083 * term;
    c3 += 0.1249888346504472 * term;
    let term = x0 * x2 * x3;
    c2 += -0.3750003908494035 * term;
    c4 += 0.46873041948400707 * term;
    let term = x0 * x2 * x2;
    c1 += 0.21875602788587467 * term;
    c3 += 0.10937615604574621 * term;
    c5 += -0.3281122359713451 * term;
    let term = x0 * x1 * x1;
    c1 += -0.2250039680398372 * term;
    c5 += 0.22500602780690765 * term;
    let term = x0 * x0 * x1;
    c2 += 0.09374762078434573 * term;
    c4 += -0.17969112227995618 * term;
    // order 4
    let term = x1 * x3.powi(3);
    c1 += -0.18752075921065503 * term;
    c3 += -0.46878370539788283 * term;
    c5 += 0.6561483898207509 * term;
    let term = x1 * x2 * x3 * x3;
    c2 += 0.35160571663396534 * term;
    c4 += 0.9844253200694091 * term;
    let term = x1 * x2 * x2 * x3;
    c1 += -0.21096139557474913 * term;
    c3 += -0.3515968293853894 * term;
    c5 += -0.4922211505899686 * term;
    let term = x1 * x2.powi(3);
    c2 += 0.05860062410005297 * term;
    c4 += 0.12305528398559931 * term;
    let term = x1.powi(3) * x3;
    c1 += -0.02148319153503821 * term;
    c3 += -0.013672695567603966 * term;
    c5 += 0.20741180368263196 * term;
    let term = x1.powi(3) * x2;
    c2 += -0.046764685064186 * term;
    c4 += 0.05820282675594289 * term;
    let term = x0 * x3.powi(3);
    c2 += -0.1250059086119771 * term;
    c4 += 0.15622099534123857 * term;
    let term = x0 * x2 * x3 * x3;
    c1 += -0.09375294724583982 * term;
    c3 += -0.2343808595231707 * term;
    c5 += 0.3280987088787046 * term;
    let term = x0 * x2 * x2 * x3;
    c2 += 0.28127178446600565 * term;
    c4 += 0.2812648053590432 * term;
    let term = x0 * x2.powi(3);
    c1 += -0.0781331829820033 * term;
    c3 += -0.15625955835252975 * term;
    let term = x0 * x1 * x1 * x3;
    c2 += -0.18515437403762713 * term;
    c4 += -0.2718990673670346 * term;
    let term = x0 * x1 * x1 * x2;
    c1 += -0.024943950594283337 * term;
    c3 += 0.11133390922504413 * term;
    c5 += 0.3339765583770631 * term;
    let term = x0 * x0 * x1 * x3;
    c1 += 0.15469497774769914 * term;
    c3 += -0.08202602878610438 * term;
    c5 += -0.07267307673610585 * term;
    let term = x0 * x0 * x1 * x2;
    c2 += -0.17422654013943273 * term;
    c4 += -0.2202954070356332 * term;
    let term = x0.powi(3) * x3;
    c2 += 0.09374677691755794 * term;
    c4 += -0.11718671831890148 * term;
    let term = x0.powi(3) * x2;
    c1 += 0.12395345531240887 * term;
    c3 += -0.07552153580599416 * term;
    c5 += -0.048427289056022366 * term;
    // order 5
    let term = x1 * x3.powi(4);
    c2 += -0.14651769110040486 * term;
    c4 += -0.4102155691051238 * term;
    let term = x1 * x2 * x3.powi(3);
    c1 += 0.29304965066280037 * term;
    c3 += 0.41025533592415714 * term;
    c5 += 1.230611682073839 * term;
    let term = x1 * x2 * x2 * x3 * x3;
    c2 += -0.3516502085181006 * term;
    c4 += -0.6153594772411939 * term;
    let term = x1 * x2.powi(3) * x3;
    c1 += 0.16117454412478233 * term;
    c3 += 0.20512528488570345 * term;
    c5 += 0.43952245235995746 * term;
    let term = x1 * x2.powi(4);
    c2 += -0.03479776959202513 * term;
    c4 += -0.05128015498128899 * term;
    let term = x1.powi(3) * x3 * x3;
    c2 += -0.030795684289689856 * term;
    c4 += -0.1646659969475305 * term;
    let term = x1.powi(3) * x2 * x3;
    c1 += 0.08907406893963096 * term;
    c3 += -0.20475693406421525 * term;
    c5 += 0.312260993249001 * term;
    let term = x1.powi(3) * x2 * x2;
    c2 += -0.02345817116944371 * term;
    c4 += 0.14383910747235487 * term;
    let term = x1.powi(5);
    c2 += -0.0017880211200527185 * term;
    c4 += -0.0032820265963874142 * term;
    let term = x0 * x3.powi(4);
    c1 += -0.04688540932810496 * term;
    c3 += -0.11720080924278639 * term;
    c5 += 0.16400750576421944 * term;
    let term = x0 * x2 * x3.powi(3);
    c2 += -0.039059703822036385 * term;
    c4 += -0.10938125730684811 * term;
    let term = x0 * x2 * x2 * x3 * x3;
    c1 += 0.17581569663679963 * term;
    c3 += 0.29301151403208353 * term;
    c5 += 0.41020559236708276 * term;
    let term = x0 * x2.powi(3) * x3;
    c2 += -0.17581783504488452 * term;
    c4 += -0.3691900532394136 * term;
    let term = x0 * x2.powi(4);
    c1 += 0.04444346041361513 * term;
    c3 += 0.06348962109448016 * term;
    c5 += 0.11427334980689796 * term;
    let term = x0 * x1 * x1 * x3 * x3;
    c1 += 0.2775923176219497 * term;
    c3 += -0.5390199674408167 * term;
    c5 += -0.40328559285443544 * term;
    let term = x0 * x1 * x1 * x2 * x3;
    c2 += -0.18724938266601043 * term;
    c4 += 0.300852660926193 * term;
    let term = x0 * x1 * x1 * x2 * x2;
    c1 += 0.10503897448013721 * term;
    c3 += -0.08987415575211627 * term;
    c5 += -0.038751621269511594 * term;
    let term = x0 * x1.powi(4);
    c1 += -0.03744703006881844 * term;
    c3 += 0.030963510293073593 * term;
    c5 += 0.05068647252493934 * term;
    let term = x0 * x0 * x1 * x3 * x3;
    c2 += 0.31291421710710754 * term;
    c4 += -0.43981282699979973 * term;
    let term = x0 * x0 * x1 * x2 * x3;
    c1 += 0.4252144036182015 * term;
    c3 += -0.8906026307056725 * term;
    c5 += -0.1125213127967789 * term;
    let term = x0 * x0 * x1 * x2 * x2;
    c2 += -0.14994177665976494 * term;
    c4 += 0.33438391703992376 * term;
    let term = x0 * x0 * x1.powi(3);
    c2 += -0.04777843196311175 * term;
    c4 += -0.08421982259446589 * term;
    let term = x0.powi(3) * x3 * x3;
    c1 += 0.11302231141726202 * term;
    c3 += 0.022141819266384886 * term;
    c5 += -0.13515523046801298 * term;
    let term = x0.powi(3) * x2 * x3;
    c2 += 0.19919226243852545 * term;
    c4 += -0.3046624288208696 * term;
    let term = x0.powi(3) * x2 * x2;
    c1 += 0.06336380205251846 * term;
    c3 += -0.2076613250208299 * term;
    c5 += 0.05352474286996386 * term;
    let term = x0.powi(3) * x1 * x1;
    c1 += -0.008942537104335754 * term;
    c3 += 0.05030142264089016 * term;
    c5 += 0.038373079309002765 * term;
    let term = x0.powi(4) * x1;
    c2 += -0.011850606939621996 * term;
    c4 += -0.015170628786762935 * term;
    let term = x0.powi(5);
    c1 += -0.0019530407383931228 * term;
    c5 += 0.0029296273884824955 * term;
    let mut chebs = Cheb::new(vec![0.0, c1, c2, c3, c4, c5]).integrate();
    chebs.zero_endpoints();
    chebs
}

fn err_using_chebs(c0: CubicBez, c1: CubicBez, f: &impl Fn(CubicBez) -> Cheb) -> f64 {
    let cheb0 = f(c0);
    let cheb1 = f(c1);
    let err = (&cheb0 - &cheb1).l2_norm();
    0.5 * err
}

pub fn run_one() -> f64 {
    let c = crate::random_cubic();
    let c1 = approx_euler_of(c);
    let dist_err = measure_cubics(c, c1);
    let cheb_err = err_using_chebs(c, c1, &chebs_of);
    cheb_err / dist_err
}

pub fn run_sample() {
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    const N: usize = 1000;
    for _ in 0..N {
        let ratio = run_one();
        s1 += ratio;
        s2 += ratio * ratio;
    }
    let mean = s1 * (1.0 / N as f64);
    let stddev = ((N as f64 * s2 - s1 * s1) / (N as f64 * (N - 1) as f64)).sqrt();
    println!("mean {mean} stddev {stddev}");
}
