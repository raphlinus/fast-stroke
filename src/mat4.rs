//! 4x4 linear system solver

// Adapted from solution posted to
// <https://stackoverflow.com/questions/23445695/efficiently-solve-ax-b-where-a-is-a-4x4-symmetric-metrix-and-b-is-4x1-vector>
// by Marek Basovník.

fn subdeterminants(mat: &[f64], vec: &[f64]) -> [f64; 10] {
    [
        mat[0] * mat[5] - mat[1] * mat[4],
        mat[0] * mat[6] - mat[2] * mat[4],
        mat[0] * mat[7] - mat[3] * mat[4],
        mat[2] * mat[7] - mat[3] * mat[6],
        mat[1] * mat[6] - mat[2] * mat[5],
        mat[1] * mat[7] - mat[3] * mat[5],
        mat[0] * vec[1] - vec[0] * mat[4],
        mat[1] * vec[1] - vec[0] * mat[5],
        mat[2] * vec[1] - vec[0] * mat[6],
        mat[3] * vec[1] - vec[0] * mat[7],
    ]
}

pub fn solve4x4(mat: &[f64; 16], vec: [f64; 4]) -> [f64; 4] {
    let a = subdeterminants(&mat[0..8], &vec[0..2]);
    let b = subdeterminants(&mat[8..16], &vec[2..4]);
    let r = [
        a[5] * b[8] + a[8] * b[5] - a[7] * b[3] - a[3] * b[7] - a[4] * b[9] - a[9] * b[4],
        a[1] * b[9] + a[9] * b[1] + a[6] * b[3] + a[3] * b[6] - a[2] * b[8] - a[8] * b[2],
        a[2] * b[7] + a[7] * b[2] - a[0] * b[9] - a[9] * b[0] - a[6] * b[5] - a[5] * b[6],
        a[0] * b[8] + a[8] * b[0] + a[6] * b[4] + a[4] * b[6] - a[7] * b[1] - a[1] * b[7],
    ];
    let idet =
        1. / (a[0] * b[3] + a[3] * b[0] + a[2] * b[4] + a[4] * b[2] - a[5] * b[1] - a[1] * b[5]);
    r.map(|x| x * idet)
}

pub fn mat4xv(mat: &[f64; 16], vec: [f64; 4]) -> [f64; 4] {
    fn dot(row: &[f64], vec: [f64; 4]) -> f64 {
        row[0] * vec[0] + row[1] * vec[1] + row[2] * vec[2] + row[3] * vec[3]
    }
    [
        dot(&mat[0..4], vec),
        dot(&mat[4..8], vec),
        dot(&mat[8..12], vec),
        dot(&mat[12..16], vec),
    ]
}

#[test]
fn matrix_solver() {
    let mut mat = [0.0; 16];
    for i in 0..16 {
        mat[i] = (i as f64).powf(0.9).cos();
    }
    let v = [2.0, 3.0, 5.0, 8.0];
    let y = solve4x4(&mat, v);
    println!("{y:?}");
    let z = mat4xv(&mat, y);
    println!("{z:?}");
    for i in 0..4 {
        assert!((z[i] - v[i]) < 2e-12);
    }
}
