use cubic_spline_rs::{CubicSpline,CubicSplineNa};
use nalgebra::DVector;
use std::io::{BufWriter,Write};
use std::fs;

fn main() {
    let x = vec![0.0,1.0,4.0,5.0,8.0];
    let y = vec![0.0,3.0,4.0,1.0,2.0];
    let mut sp = CubicSpline::new(x.clone(),y.clone());
    sp.calc_spline();
    let t = (0..1000).map(|t| t as f64 * 0.01).collect();
    let intr = sp.interpolate_vec(&t);

    let mut sp2 = CubicSplineNa::new(x.clone(),y.clone());
    sp2.calc_spline();
    let t2 = DVector::from_fn(1000, |t, _| t as f64 * 0.01);

    let intr2 = sp2.interpolate_vec(&t2);

    let mut out = BufWriter::new(fs::File::create("out.csv").unwrap());
    writeln!(out, "x,y,y2").unwrap();
    for i in 0..1000 {
        writeln!(out, "{},{},{}", t[i], intr[i], intr2[i]).unwrap();
    }

    let mut out = BufWriter::new(fs::File::create("org.csv").unwrap());
    writeln!(out, "x,y").unwrap();
    for i in 0..5 {
        writeln!(out, "{},{}", x[i], y[i]).unwrap();
    }

    let (a,b,c,d) = sp.calc_spline_coef();
    println!("a: {:?}\nb: {:?}\nc: {:?}\nd: {:?}", a,b,c,d);
}
