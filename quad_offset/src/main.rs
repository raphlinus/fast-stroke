use xilem_web::{
    elements::{
        html::{div, input, label},
        svg::{g, svg},
    },
    input_event_target_value,
    interfaces::{Element, HtmlInputElement, SvgGeometryElement, SvgPathElement},
    svg::{
        kurbo::{BezPath, Circle, Line, ParamCurve, Point, QuadBez, Shape, Vec2},
        peniko,
    },
    App, DomView, PointerMsg,
};

#[derive(Default)]
struct AppState {
    p0: Point,
    p1: Point,
    p2: Point,
    grab: GrabState,
    offset: f64,
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

/// Unit normal (rotated 90° ccw from tangent) of a QuadBez at parameter t.
fn quad_normal(q: &QuadBez, t: f64) -> Vec2 {
    let d0 = q.p1 - q.p0;
    let d1 = q.p2 - q.p1;
    let tangent = 2.0 * ((1.0 - t) * d0 + t * d1);
    let len = tangent.length();
    if len < 1e-10 {
        Vec2::new(0.0, 1.0)
    } else {
        Vec2::new(-tangent.y / len, tangent.x / len)
    }
}

/// Approximate the offset of a QuadBez by sampling at t=0, 0.5, 1 and
/// fitting a new QuadBez through the three offset points.
fn offset_quad(q: &QuadBez, d: f64) -> QuadBez {
    let p0_off = q.p0 + d * quad_normal(q, 0.0);
    let p_mid_off = q.eval(0.5) + d * quad_normal(q, 0.5);
    let p2_off = q.p2 + d * quad_normal(q, 1.0);
    // B(0.5) = (P0 + 2*P1 + P2)/4  =>  P1 = 2*B(0.5) - (P0+P2)/2
    let p1_off = Point::new(
        2.0 * p_mid_off.x - 0.5 * (p0_off.x + p2_off.x),
        2.0 * p_mid_off.y - 0.5 * (p0_off.y + p2_off.y),
    );
    QuadBez::new(p0_off, p1_off, p2_off)
}

fn app_logic(state: &mut AppState) -> impl DomView<AppState> {
    let q = QuadBez::new(state.p0, state.p1, state.p2);
    let d = state.offset;

    let path: BezPath = q.to_path(0.0);
    let stroke = xilem_web::svg::kurbo::Stroke::new(2.0);
    let stroke_thin = xilem_web::svg::kurbo::Stroke::new(1.5);

    let white: peniko::Color = peniko::Color::WHITE;
    let lime: peniko::Color = peniko::color::palette::css::LIME;
    let gray: peniko::Color = peniko::color::palette::css::GRAY;
    let transparent: peniko::Color = peniko::Color::TRANSPARENT;

    let path_offset: BezPath = offset_quad(&q, d).to_path(0.0);

    let svg_el = svg(g((
        Line::new(state.p0, state.p1).stroke(gray, stroke.clone()),
        Line::new(state.p1, state.p2).stroke(gray, stroke.clone()),
        path.stroke(white, stroke_thin.clone()).fill(transparent),
        path_offset.stroke(lime, stroke_thin.clone()).fill(transparent),
        g((
            Circle::new(state.p0, 6.0)
                .pointer(|s: &mut AppState, msg| s.grab.handle(&mut s.p0, &msg)),
            Circle::new(state.p1, 6.0)
                .pointer(|s: &mut AppState, msg| s.grab.handle(&mut s.p1, &msg)),
            Circle::new(state.p2, 6.0)
                .pointer(|s: &mut AppState, msg| s.grab.handle(&mut s.p2, &msg)),
        )),
    )))
    .attr("width", 900)
    .attr("height", 600);

    let slider = div((
        label(format!("offset: {:.1}", d)).attr("for", "offset-slider"),
        input(())
            .attr("id", "offset-slider")
            .type_("range")
            .attr("min", "-150")
            .attr("max", "150")
            .attr("step", "1")
            .attr("value", d as i32)
            .on_input(|s: &mut AppState, e| {
                if let Some(v) = input_event_target_value(&e) {
                    s.offset = v.parse().unwrap_or(0.0);
                }
            }),
    ))
    .attr("id", "offset-control");

    div((svg_el, slider)).attr("id", "beztoy-container-inner")
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
    state.p0 = Point::new(150.0, 450.0);
    state.p1 = Point::new(450.0, 100.0);
    state.p2 = Point::new(750.0, 450.0);
    state.offset = 10.0;

    App::new(xilem_web::document_body(), state, app_logic).run();
}
