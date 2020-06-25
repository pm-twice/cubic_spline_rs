/// Bスプライン補間を行うクラス
/// 
/// 3次スプライン補間とは異なり、制御点は必ずしも通らない。
/// 
/// ## 基本式
/// S(u) = sigma_i=0^{p-1} Pi bi,n(u)  
/// u ∈ [u0,...,um-1]  
/// Pi: 制御点  
/// bi,n: 基底関数  
/// u: ノットベクトル  
/// 
/// ### 基底関数
/// ノットベクトルuから定義が可能。
/// 
/// j = 0,...,m-2に対して、
/// bj,0(u) =   
///     1   (uj < u < uj+1)
///     0   (otherwise)
/// 
/// j = 0,1,...,m-k-2に対して、
/// bj,k(u) = (u-uj) / (uj+k -uj) bj,k-1(u)
///     + (uj+k+1 -u) / (uj+k+1 -uj+1) bj+1,k-1(u)
/// 
/// ### ノットベクトル
/// 間一様ノットベクトルを用いる。
/// 制御点の次数pと、Bスプラインの次数nに対して、
/// 最初のn+1個を0,最後のn+1個を1とする
/// 残りのm-2(n+1)個は0より大きく1より小さい値で均等間隔とする
/// 
/// ノットの個数mは、p=m-n-1より、p+n+1で表せる。
/// 
/// ## 参考
/// - [pages.mtu.edu](https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/CURVE-APP-global.html)
/// - [wikipedia](https://ja.wikipedia.org/wiki/B-%E3%82%B9%E3%83%97%E3%83%A9%E3%82%A4%E3%83%B3%E6%9B%B2%E7%B7%9A#:~:text=B%2D%E3%82%B9%E3%83%97%E3%83%A9%E3%82%A4%E3%83%B3%E6%9B%B2%E7%B7%9A%EF%BC%88B%E3%82%B9%E3%83%97%E3%83%A9%E3%82%A4%E3%83%B3,%E7%AD%89%E3%81%AE%E6%80%A7%E8%B3%AA%E3%81%8C%E3%81%82%E3%82%8B%E3%80%82&text=%E6%9B%B2%E7%B7%9A%E3%81%AF%E5%BF%85%E3%81%9A%E3%81%97%E3%82%82%E5%88%B6%E5%BE%A1%E7%82%B9%E3%82%92%E9%80%9A%E3%82%89%E3%81%AA%E3%81%84%E3%80%82)
/// - [緑川氏の解説](https://shoichimidorikawa.github.io/Lec/CG-Math/Bspline.pdf)
/// - [旅行好きなソフトエンジニアの備忘録](http://ni4muraano.hatenablog.com/entry/2017/10/30/212326)
/// - [理科系の勉強日記](https://kenbell.hatenablog.com/entry/2019/03/03/221925)
pub struct BSpline {
    x: Vec<f64>,
    y: Vec<f64>,
    degree: usize,
    knot: Vec<f64>,
}

impl BSpline {
    pub fn new_uniform(deg: usize, x: Vec<f64>, y: Vec<f64>) -> BSpline {
        let size = x.len();
        if size != y.len() {
            panic!("x & y must be same size");
        }
        let mut b = BSpline {
            x: x,
            y: y,
            degree: deg,
            knot: vec![],
        };
        b.calc_open_uniform_knot();
        b
    }

    /// 間一様ノットベクトルの計算を行う。(open-uniform knot vector)
    /// deg: スプライン補間の次数
    /// p: 制御点の個数
    /// 曲線の端点が最初と最後の制御点となる特徴がある。
    fn calc_open_uniform_knot(&mut self) {
        let deg = self.degree;
        let p = self.x.len();

        let m = deg + p + 1;    // ノットベクトルの次数
        self.knot = vec![0.0; m];
        // 最初のdeg+1個は0

        // 最後のdeg+1個は1
        let idx = m - (deg+1);
        for i in idx..m {
            self.knot[i] = 1.0;
        }

        // その間は等間隔
        let num = m - 2 *(deg+1) + 1;
        let dx = 1.0/num as f64;
        for i in deg+1..idx {
            self.knot[i] = self.knot[i-1] + dx;
        }
    }

    // 規定関数を計算する。引数uを再帰計算するのが非効率
    // 0<=u<=1
    fn calc_basis(&self, j: usize, k: usize, u: f64) -> f64 { 
        if k == 0 {
            // j = 0,...,m-2に対して、
            // bj,0(u) =   
            //     1   (uj < u < uj+1)
            //     0   (otherwise)

            if self.knot[j] <= u && u < self.knot[j+1] {
                1.0
            } else {
                0.0
            }
        } else {
            // j = 0,1,...,m-k-2に対して、
            // bj,k(u) = (u-uj) / (uj+k -uj) bj,k-1(u)
            //     + (uj+k+1 -u) / (uj+k+1 -uj+1) bj+1,k-1(u)
            let w1 = if self.knot[j+k] - self.knot[j] != 0.0 {  // 0除算の回避
                (u-self.knot[j]) / (self.knot[j+k] - self.knot[j]) * self.calc_basis(j, k-1, u)
            } else {
                0.0
            };
            let w2 = if self.knot[j+k+1] - self.knot[j+1] != 0.0 {  // 0除算の回避
                (self.knot[j+k+1]-u) / (self.knot[j+k+1] - self.knot[j+1]) * self.calc_basis(j+1, k-1, u) 
            } else {
                0.0
            };
            w1+w2
        }
    }

    
    pub fn interpolate_vec(&self, t: &Vec<f64>) -> Vec<f64> {
        t.iter().map(|v| self.interpolate(*v)).collect()
    }

    /// 与えられたtに対して、計算された係数からBスプライン補間を行う関数
    pub fn interpolate(&self, x: f64) -> f64 {
        // S(u) = sigma_i=0^{p-1} Pi bi,n(u)  
        let n = self.x.len();

        // xを0~1のtにスケーリング
        let beg = self.x[0];
        let end = self.x[n-1];
        let t = if x <= beg {
            return self.y[0]
        } else if x >= end {
            return self.y[n-1]
        } else {
            (x - beg)/(end-beg)
        };

        let mut y = 0.0;
        let mut ax = 0.0;
        for i in 0..n {
            y += self.y[i] * self.calc_basis(i, self.degree, t);
            ax += self.x[i] * self.calc_basis(i, self.degree, t);
        }

        assert!((ax-x).abs() < 1e-2);

        y
    }

    pub fn get_basis_vec(&self, t: &Vec<f64>, dim: usize) -> Vec<f64> {
        t.iter().map(|v| self.get_basis(*v, dim)).collect()
    }

    pub fn get_basis(&self, x: f64, dim: usize) -> f64 {
        let n = self.x.len();

        // xを0~1のtにスケーリング
        let beg = self.x[0];
        let end = self.x[n-1];
        let t = if x <= beg {
            return self.y[0]
        } else if x >= end {
            return self.y[n-1]
        } else {
            (x - beg)/(end-beg)
        };

        self.calc_basis(dim, self.degree, t)
    }
}


// #[cfg(test)]
// mod tests {
//     use super::BSpline;

//     #[test]
//     fn it_works() {
//         let k = BSpline::calc_knot(2, 4);
//         assert_eq!(k, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]);

//         let k = BSpline::calc_knot(3, 6);
//         assert_eq!(k, vec![0.0, 0.0, 0.0, 0.0, 0.3333333333333333, 0.6666666666666666, 1.0, 1.0, 1.0, 1.0]);
//     }
// }
