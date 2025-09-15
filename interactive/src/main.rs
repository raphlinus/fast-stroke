use kurbo::{BezPath, Circle, CubicBez, Line, PathEl, Point, Shape};
use perturb::CurveOffset;
use xilem_web::{
    elements::{
        html::div,
        svg::{g, svg},
    },
    interfaces::{Element, SvgGeometryElement, SvgPathElement},
    svg::peniko::Color,
    App, DomView, PointerMsg,
};

mod cusp;
mod evolute;
mod offset;
mod perturb;

#[derive(Default)]
struct AppState {
    p0: Point,
    p1: Point,
    p2: Point,
    p3: Point,
    offset: f64,
    tolerance: f64,
    grab: GrabState,
}

#[derive(Default)]
struct GrabState {
    is_down: bool,
    id: i32,
    dx: f64,
    dy: f64,
}

impl GrabState {
    fn handle(&mut self, pt: &mut Point, p: &PointerMsg) {
        match p {
            PointerMsg::Down(e) => {
                if e.button == 0 {
                    self.dx = pt.x - e.position.x;
                    self.dy = pt.y - e.position.y;
                    self.id = e.id;
                    self.is_down = true;
                }
            }
            PointerMsg::Move(e) => {
                if self.is_down && self.id == e.id {
                    pt.x = (self.dx + e.position.x).min(850.).max(8.);
                    pt.y = (self.dy + e.position.y).min(592.).max(8.);
                }
            }
            PointerMsg::Up(e) => {
                if self.id == e.id {
                    self.is_down = false;
                }
            }
        }
    }
}

fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    let c = CubicBez::new(state.p0, state.p1, state.p2, state.p3);
    let path = c.to_path(0.0);
    let stroke = xilem_web::svg::kurbo::Stroke::new(2.0);
    let stroke_thin = xilem_web::svg::kurbo::Stroke::new(2.0);
    let d = 100.0;
    //perturb::scaling_test(c, d);
    let co = CurveOffset::new(c);
    let (a, b) = perturb::draw_arc(c);
    let mut soln_arc_lse = perturb::OffsetSolutionLse::from_a_b(a, b, d);
    let c_arc = soln_arc_lse.apply(&co);
    let path_arc = c_arc.to_path(0.0);
    let err_arc = perturb::plot(&perturb::error_by_rays(c, d, c_arc));

    let (a, b) = perturb::arc_onept_linear(c);
    let soln_aol_lse = perturb::OffsetSolutionLse::from_a_b(a, b, d);
    let c_aol = soln_aol_lse.apply(&co);
    let path_aol = c_aol.to_path(0.0);
    let err_aol = perturb::plot(&perturb::error_by_rays(c, d, c_aol));

    for _ in 0..1 {
        soln_arc_lse.newton_step_rev(&co);
        soln_arc_lse.refine_lse(&co);
    }
    let c_arc_lse = soln_arc_lse.apply(&co);
    let path_arc_lse = c_arc_lse.to_path(0.0);
    let err_arc_lse_refined = perturb::plot(&perturb::error_by_rays(c, d, c_arc_lse));

    /*
    let (a1, b1) = perturb::one_point(c);
    let mut soln_one_point = perturb::OffsetSolution::from_a_b(a1, b1, d);
    let [y3_0, y3_1, y3_2] = soln_one_point.refine_ts(&co);
    let err_one_point = y3_0.max(y3_1).max(y3_2);
    if err_one_point < best_err {
        best_soln = &soln_one_point;
        best_err = err_center;
    }

    let c_offset = best_soln.apply(&co);
    let path_offset = c_offset.to_path(0.0);
    let err_offset = perturb::plot(&perturb::error_by_rays(c, d, c_offset));

    let (a2, b2) = perturb::one_point_at(c, d, 0.5);
    let c_one_point = co.apply(a2, b2, d);
    let path_one_point = c_one_point.to_path(0.0);
    let err_one_point = perturb::plot(&perturb::error_by_rays(c, d, c_one_point));

    let t = perturb::brute_one_point(c, d, 0.5);
    let (a3, b3) = perturb::one_point_at(c, d, t);
    let c_one_point2 = co.apply(a3, b3, d);
    let err_one_point2 = perturb::plot(&perturb::error_by_rays(c, d, c_one_point2));
    */

    let tolerance = 0.25;
    let path_offset = offset::offset_cubic(c, d, tolerance);

    //let evo = evolute::evolute_hacky_approx(c);
    //let evc = evolute::evolute_approx(c, tolerance);

    const NONE: Color = Color::TRANSPARENT;
    const HANDLE_RADIUS: f64 = 6.0;
    let svg_el = svg(g((
        Line::new(state.p0, state.p1).stroke(Color::BLUE, stroke.clone()),
        Line::new(state.p2, state.p3).stroke(Color::BLUE, stroke.clone()),
        //Line::new((100., 200.), (600., 200.)).stroke(Color::GREEN, stroke.clone()),
        /*
        Line::new((100., 200. - y0), (267., 200. - y0)).stroke(Color::LIME, stroke.clone()),
        Line::new((100., 200. + y0), (267., 200. + y0)).stroke(Color::LIME, stroke.clone()),
        Line::new((267., 200. - y1), (433., 200. - y1)).stroke(Color::LIME, stroke.clone()),
        Line::new((267., 200. + y1), (433., 200. + y1)).stroke(Color::LIME, stroke.clone()),
        Line::new((433., 200. - y2), (600., 200. - y2)).stroke(Color::LIME, stroke.clone()),
        Line::new((433., 200. + y2), (600., 200. + y2)).stroke(Color::LIME, stroke.clone()),
        */
        path.stroke(Color::WHITE, stroke_thin.clone()).fill(NONE),
        subdiv_pts(&path_offset),
        path_offset
            .stroke(Color::LIME, stroke_thin.clone())
            .fill(NONE),
        path_arc
            .stroke(Color::ORANGE, stroke_thin.clone())
            .fill(NONE),
        path_aol.stroke(Color::BLUE, stroke_thin.clone()).fill(NONE),
        /*
        path_arc_lse
            .stroke(Color::YELLOW, stroke_thin.clone())
            .fill(NONE),
        //path_lse.stroke(Color::LIME, stroke_thin.clone()).fill(NONE),
        */
        err_arc
            .stroke(Color::ORANGE, stroke_thin.clone())
            .fill(NONE),
        err_aol.stroke(Color::BLUE, stroke_thin.clone()).fill(NONE),
        /*
        err_arc_lse_refined
            .stroke(Color::YELLOW, stroke_thin.clone())
            .fill(NONE),
        */
        //evo.stroke(Color::GREEN, stroke_thin.clone()).fill(NONE),
        //evc.stroke(Color::RED, stroke_thin.clone()).fill(NONE),
        g((
            Circle::new(state.p0, HANDLE_RADIUS)
                .pointer(|s: &mut AppState, msg| s.grab.handle(&mut s.p0, &msg)),
            Circle::new(state.p1, HANDLE_RADIUS)
                .pointer(|s: &mut AppState, msg| s.grab.handle(&mut s.p1, &msg)),
            Circle::new(state.p2, HANDLE_RADIUS)
                .pointer(|s: &mut AppState, msg| s.grab.handle(&mut s.p2, &msg)),
            Circle::new(state.p3, HANDLE_RADIUS)
                .pointer(|s: &mut AppState, msg| s.grab.handle(&mut s.p3, &msg)),
        )),
    )))
    .attr("width", 900)
    .attr("height", 600);
    div((svg_el,)).attr("id", "beztoy-container-inner")
}

fn subdiv_pts(path: &BezPath) -> impl DomView<AppState> {
    let mut circles = vec![];
    for el in path.elements() {
        match el {
            PathEl::MoveTo(p) => circles.push(Circle::new(*p, 5.0)),
            PathEl::CurveTo(_, _, p3) => circles.push(Circle::new(*p3, 5.0)),
            _ => (),
        }
    }
    g(circles).class("subdiv")
}

pub fn main() {
    use tracing_subscriber::{fmt::format::Pretty, prelude::*};
    use tracing_web::{performance_layer, MakeWebConsoleWriter};

    console_error_panic_hook::set_once();

    let fmt_layer = tracing_subscriber::fmt::layer()
        .pretty()
        .with_ansi(false)
        .without_time()
        .with_writer(MakeWebConsoleWriter::new());
    let perf_layer = performance_layer().with_details_from_fields(Pretty::default());

    tracing_subscriber::registry()
        .with(fmt_layer)
        .with(perf_layer)
        .init();

    let mut state = AppState::default();
    state.p0 = Point::new(55.0, 466.0);
    state.p1 = Point::new(350.0, 146.0);
    state.p2 = Point::new(496.0, 537.0);
    state.p3 = Point::new(739.0, 244.0);
    // state.p0 = Point::new(742.483763921753, 245.69451584513587);
    // state.p1 = Point::new(742.3156048952269, 245.37897448914785);
    // state.p2 = Point::new(741.1749624891498, 244.83296435754676);
    // state.p3 = Point::new(739.0, 244.0);
    state.offset = 100.;
    state.tolerance = 1.;

    App::new(xilem_web::document_body(), state, app_logic).run();
}
