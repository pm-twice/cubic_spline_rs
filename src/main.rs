use cubic_spline_rs::{CubicSpline,CubicSplineNa,LinearInterpolation};
use nalgebra::DVector;
use std::io::{BufWriter,Write};
use std::fs;

fn main() {
    let x = vec![0.0,1.0,4.0,5.0,8.0];
    let y = vec![0.0,3.0,4.0,1.0,2.0];

    // Vecを用いたスプライン補間
    let mut sp = CubicSpline::new(x.clone(),y.clone());
    sp.calc_spline();
    let t = (0..1000).map(|t| t as f64 * 0.01 -1.0).collect();
    let intr = sp.interpolate_vec(&t);

    // Nalgebraを用いたスプライン補間
    let mut sp2 = CubicSplineNa::new(x.clone(),y.clone());
    sp2.calc_spline();
    let t2 = DVector::from_fn(1000, |t, _| t as f64 * 0.01 -1.0);
    let intr2 = sp2.interpolate_vec(&t2);

    // 線形補間
    let sp3 = LinearInterpolation::new(x.clone(),y.clone());
    let t3 = (0..1000).map(|t| t as f64 * 0.01 -1.0).collect();
    let intr3 = sp3.interpolate_vec(&t3);

    // 補間結果の保存
    let mut out = BufWriter::new(fs::File::create("out.csv").unwrap());
    writeln!(out, "x,sp(vec),sp(na),linear").unwrap();
    for i in 0..1000 {
        writeln!(out, "{},{},{},{}", t[i], intr[i], intr2[i], intr3[i]).unwrap();
    }

    // 元のデータを保存
    let mut out = BufWriter::new(fs::File::create("org.csv").unwrap());
    writeln!(out, "x,y").unwrap();
    for i in 0..5 {
        writeln!(out, "{},{}", x[i], y[i]).unwrap();
    }

    // スプライン補間により計算された3次式の係数表示
    let (a,b,c,d) = sp.calc_spline_coef();
    println!("a: {:?}\nb: {:?}\nc: {:?}\nd: {:?}", a,b,c,d);
}
