use kurbo::{Circle, CubicBez, Line, Point, Shape};
use xilem_web::{
    elements::{
        html::div,
        svg::{g, svg},
    },
    interfaces::{Element, SvgGeometryElement, SvgPathElement},
    svg::peniko::Color,
    App, DomView, PointerMsg,
};

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
    let d = 50.0;
    //perturb::scaling_test(c, d);
    let delta = perturb::two_point_approx(c);
    let c_offset = CubicBez::new(
        c.p0 + d * delta.p0.to_vec2(),
        c.p1 + d * delta.p1.to_vec2(),
        c.p2 + d * delta.p2.to_vec2(),
        c.p3 + d * delta.p3.to_vec2(),
    );
    let path_offset = c_offset.to_path(0.0);
    let err2 = perturb::plot(&perturb::error_by_rays(c, d, c_offset));

    let delta_minmax = perturb::linear_minmax(c);
    let c_minmax = CubicBez::new(
        c.p0 + d * delta_minmax.p0.to_vec2(),
        c.p1 + d * delta_minmax.p1.to_vec2(),
        c.p2 + d * delta_minmax.p2.to_vec2(),
        c.p3 + d * delta_minmax.p3.to_vec2(),
    );
    let path_minmax = c_minmax.to_path(0.0);
    let err_minmax = perturb::plot(&perturb::error_by_rays(c, d, c_minmax));
    let error = perturb::plot_error(c, delta_minmax);
    let spline_error = perturb::spline_error(c, delta_minmax);
    let [y0, y1, y2] = perturb::est_err_bounds(c, delta_minmax);

    const NONE: Color = Color::TRANSPARENT;
    const HANDLE_RADIUS: f64 = 6.0;
    let svg_el = svg(g((
        Line::new(state.p0, state.p1).stroke(Color::BLUE, stroke.clone()),
        Line::new(state.p2, state.p3).stroke(Color::BLUE, stroke.clone()),
        Line::new((100., 200.), (600., 200.)).stroke(Color::GREEN, stroke.clone()),
        Line::new((100., 200. - y0), (267., 200. - y0)).stroke(Color::LIME, stroke.clone()),
        Line::new((100., 200. + y0), (267., 200. + y0)).stroke(Color::LIME, stroke.clone()),
        Line::new((267., 200. - y1), (433., 200. - y1)).stroke(Color::LIME, stroke.clone()),
        Line::new((267., 200. + y1), (433., 200. + y1)).stroke(Color::LIME, stroke.clone()),
        Line::new((433., 200. - y2), (600., 200. - y2)).stroke(Color::LIME, stroke.clone()),
        Line::new((433., 200. + y2), (600., 200. + y2)).stroke(Color::LIME, stroke.clone()),
        path.stroke(Color::WHITE, stroke_thin.clone()).fill(NONE),
        path_minmax
            .stroke(Color::YELLOW, stroke_thin.clone())
            .fill(NONE),
        error.stroke(Color::RED, stroke_thin.clone()).fill(NONE),
        // spline_error
        // .stroke(Color::MEDIUM_ORCHID, stroke_thin.clone())
        // .fill(NONE),
        // err2.stroke(Color::ORANGE, stroke_thin.clone()).fill(NONE),
        // err_minmax
        //     .stroke(Color::LIME, stroke_thin.clone())
        //     .fill(NONE),
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
    state.offset = 100.;
    state.tolerance = 1.;

    App::new(xilem_web::document_body(), state, app_logic).run();
}
