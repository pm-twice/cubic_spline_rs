use cubic_spline_rs::{CubicSpline,CubicSplineNa,LinearInterpolation,BSpline};
use nalgebra::DVector;
use std::io::{BufWriter,Write};
use std::fs;

fn main() {
    let x = vec![0.0,1.0,4.0,5.0,8.0];
    let y = vec![0.5,3.0,4.0,1.0,2.0];
    let n = x.len();

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
    let intr3 = sp3.interpolate_vec(&t);

    // Bスプライン(1次)
    // let sp4 = BSpline::new_uniform(1,x.clone(),y.clone());
    // let intr4 = sp4.interpolate_vec(&t);
    // let bas4: Vec<Vec<f64>> = (0..n).map(|i| sp4.get_basis_vec(&t, i)).collect();

    // Bスプライン(2次)
    let sp5 = BSpline::new_uniform(2,x.clone(),y.clone());
    let intr5 = sp5.interpolate_vec(&t);
    let bas5: Vec<Vec<f64>> = (0..n).map(|i| sp5.get_basis_vec(&t, i)).collect();

    // Bスプライン(3次)
    // let sp6 = BSpline::new_uniform(3,x.clone(),y.clone());
    // let intr6 = sp6.interpolate_vec(&t);
    // let bas6: Vec<Vec<f64>> = (0..n).map(|i| sp6.get_basis_vec(&t, i)).collect();

    // // 補間結果の保存
    // let mut out = BufWriter::new(fs::File::create("out.csv").unwrap());
    // writeln!(out, "x,sp(vec),sp(na),linear,bsp(1d),bsp(2d),bsp(3d)").unwrap();
    // for i in 0..1000 {
    //     writeln!(out, "{},{},{},{},{},{},{}", t[i], intr[i], intr2[i], intr3[i], intr4[i], intr5[i], intr6[i]).unwrap();
    // }

    // // 基底関数の出力
    // let mut out = BufWriter::new(fs::File::create("bas1.csv").unwrap());
    // for i in 0..1000 {
    //     write!(out, "{}", t[i]).unwrap();
    //     for j in 0..n {
    //         write!(out, ",{}", bas4[j][i]).unwrap();
    //     }
    //     write!(out, "\n").unwrap();
    // }
    // // 基底関数の出力
    // let mut out = BufWriter::new(fs::File::create("bas2.csv").unwrap());
    // for i in 0..1000 {
    //     write!(out, "{}", t[i]).unwrap();
    //     for j in 0..n {
    //         write!(out, ",{}", bas5[j][i]).unwrap();
    //     }
    //     write!(out, "\n").unwrap();
    // }
    // // 基底関数の出力
    // let mut out = BufWriter::new(fs::File::create("bas3.csv").unwrap());
    // for i in 0..1000 {
    //     write!(out, "{}", t[i]).unwrap();
    //     for j in 0..n {
    //         write!(out, ",{}", bas6[j][i]).unwrap();
    //     }
    //     write!(out, "\n").unwrap();
    // }

    // 元のデータを保存
    let mut out = BufWriter::new(fs::File::create("org.csv").unwrap());
    writeln!(out, "x,y").unwrap();
    for i in 0..5 {
        writeln!(out, "{},{}", x[i], y[i]).unwrap();
    }

    // // スプライン補間により計算された3次式の係数表示
    // let (a,b,c,d) = sp.calc_spline_coef();
    // println!("a: {:?}\nb: {:?}\nc: {:?}\nd: {:?}", a,b,c,d);
}
