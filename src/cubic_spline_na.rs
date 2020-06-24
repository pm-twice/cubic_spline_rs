use nalgebra::{DMatrix,DVector};
/// # 3次スプライン補間を行うクラス。  
/// Nalgebraを用いる。
/// 
/// 点群(x_i,y_i): (i=1,...,N+1)に対して、
/// 
/// i=1,...,Nの各離散点x_i,x_{i+1}の間を補間する3次多項式は、
/// S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3  
/// ( x_i <= x <= x_{i+1} )
/// として4つのパラメータa_i,b_i,c_i,d_iを用いて定義される。
/// 
/// この時、求める未知数は4Nとなる。  
/// 制約条件としては次の通り。
/// 1. S_i(x)は各点(y_i,y_{i+1})を通る: 2N個の方程式
/// 2. 両端を除く各点で1階微分が連続: N-1個の方程式
/// 3. 両端を除く各点で2階微分が連続: N-1個の方程式
/// 4. 両端の点において、2階の微分係数が0: 2個の方程式
/// 
/// 条件1~3だけでは足りないので、条件4を追加することで、
/// 4N個の制約条件から連立方程式の問題に帰着できる。
/// 
/// ## 1回,2回微分の式
/// 微分下式は次の通り。  
/// 
/// S_i^(1)(x) = b_i + 2c_i(x-x_i) + 3d_i(x-x_i)^2  
/// S_i^(2)(x) = 2c_i + 6d_i(x-x_i)
/// 
/// ## 制約条件の計算
/// h_i = xi+1-x_iとおく。
/// 
/// ### 条件1
/// S_i(x_i) = a_i = y_i    (i=1,...,N)
/// S_i(x_i+1) = a_i + b_i h_i + c_i h_i^2 + d_i h_i ^3 = y_i+1 (i=1,...,N)
/// 
/// 
/// よって、
/// a_i = y_i (i=1,...,N+1)
/// 
/// また、下記の式をb_iについて解くと、(i=1,...,N)
/// a_i + b_i h_i + c_i h_i^2 + d_i h_i ^3 = a_i+1
/// 
/// b_i = (a_i+1 - a_i)/h_i - c_i h_i - d_i h_i^2
/// 
/// ### 条件2
/// S_i^(1)(x_i+1) = S_i+1^(1)(x_i+1) より、
/// b_i + 2c_i h_i + 3d_i h_i^2 = b_i+1
/// 
/// 後ほどb_i,d_iを除去する。
/// 
/// ### 条件3
/// S_i^(2)(x_i+1) = S_i+1^(2)(x_i+1) より、
/// 2c_i + 6d_i h_i = 2c_i+1 から、
/// 
/// i=1,...,N-1に対して、
/// c_i + 3d_i h_i = c_i+1 
/// 
/// d_iについて解くと、
/// d_i = (c_i+1 - c_i)/3h_i
/// 
/// 
/// ### 条件2への代入
/// b_i,d_iを除去する。
/// 
/// b_i + 2c_i h_i + 3d_i h_i^2 = b_i+1
/// へb_iを代入して、
/// 
/// (a_i+1 - a_i)/h_i - c_i h_i - d_i h_i^2 + 2c_i h_i + 3d_i h_i^2 
/// = (a_i+2 - a_i+1)/h_i+1 - c_i+1 h_i+1 - d_i+1 h_i+1^2
/// ⇔
/// (a_i+1 - a_i)/h_i + c_i h_i + 2d_i h_i^2 
/// = (a_i+2 - a_i+1)/h_i+1 - c_i+1 h_i+1 - d_i+1 h_i+1^2
/// さらにd_iを代入して
/// 
/// (a_i+1 - a_i)/h_i + c_i h_i + 2/3 (c_i+1 - c_i) h_i  
/// = (a_i+2 - a_i+1)/h_i+1 - c_i+1 h_i+1 - (c_i+2 - c_i+1) h_i+1 /3
/// ⇔
/// (a_i+1 - a_i)/h_i + (2c_i+1 + c_i) h_i/3  
/// = (a_i+2 - a_i+1)/h_i+1 - (c_i+2 + 2c_i+1) h_i+1 /3
/// ⇔
/// (2c_i+1 + c_i) h_i/3 + (c_i+2 + 2c_i+1) h_i+1 /3
/// = -(a_i+1 - a_i)/h_i + (a_i+2 - a_i+1)/h_i+1
/// ⇔
/// (2c_i+1 + c_i) h_i + (c_i+2 + 2c_i+1) h_i+1
/// = -3(a_i+1 - a_i)/h_i + 3(a_i+2 - a_i+1)/h_i+1
/// ⇔
/// h_i+1c_i+2 + 2(h_i+h_i+1)c_i+1 + h_ic_i
/// = -3(a_i+1 - a_i)/h_i + 3(a_i+2 - a_i+1)/h_i+1
/// 
/// c_iの式に直せた。
/// 
/// 更に修正して、i=2,...,Nに対して、
/// h_ic_i+1 + 2(h_i-1+h_i)c_i + h_i-1c_i-1
/// = -3(a_i - a_i-1)/h_i-1 + 3(a_i+1 - a_i)/h_i
/// 
/// ### 条件4
/// c_1 = 0
/// c_N+1 = 0
/// 
/// ### まとめ
/// a,bについては
/// 
/// 上記をまとめた行列形式は次の通り。
/// c = [c_1, c_2, ..., c_N,c_N+1] (N+1次)
/// k = [0, 
///     -3(a_2 - a_1)/h_1 + 3(a_3 - a_2)/h_2 ,
///     -3(a_3 - a_2)/h_2 + 3(a_4 - a_3)/h_3 ,
///     ...
///     -3(a_N-1 - a_N-2)/h_N-2 + 3(a_N - a_N-1)/h_N-1 ,
///     -3(a_N - a_N-1)/h_N-1 + 3(a_N+1 - a_N)/h_N ,
///     ,0] (N+1次)
/// M = [
/// 1, 0, 0, ...
/// h_1, 2(h_1+h_2), h2, 0... ,0
/// 0, h_2, 2(h_2+h_3), h3, 0, ...,0
/// ...
/// 0, ... h_N-2, 2(h_N-2+h_N-1), h_N-1, 0
/// 0, ... 0, h_N-1, 2(h_N-1+h_N), h_N
/// 0, ...,0,   0 ,         0, 1
/// ] (N+1xN+1次)
/// 
/// Mc = k としてcを求め、そこからd,bを算出する。
/// なお、h_i = xi+1-x_i
/// 
/// ### 参考
/// [TajimaRobotics](https://tajimarobotics.com/cubic-spline-interpolation-2/)  
/// [TajimaRobotics](https://tajimarobotics.com/cubic-spline-interpolation-3/)
/// [fac.ksu.edu.sa](https://fac.ksu.edu.sa/sites/default/files/numerical_analysis_9th.pdf#page=167)
/// [その実装](https://gist.github.com/svdamani/1015c5c4b673c3297309)
/// 
pub struct CubicSplineNa {
    x: Vec<f64>,
    y: Vec<f64>,
    a: DVector<f64>,
    b: DVector<f64>,
    c: DVector<f64>,
    d: DVector<f64>,
}


impl CubicSplineNa {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> CubicSplineNa {
        let size = x.len();
        if size != y.len() {
            panic!("x & y must be same size");
        }
        CubicSplineNa {
            x: x,
            y: y,

            a: DVector::zeros(1),
            b: DVector::zeros(1),
            c: DVector::zeros(1),
            d: DVector::zeros(1),
        }
    }

    pub fn interpolate_vec(&self, t: &DVector<f64>) -> DVector<f64> {
        let n = t.len();
        DVector::from_iterator(n, t.iter().map(|v| self.interpolate(*v)) )
    }

    /// 与えられたtに対して、計算された係数からスプライン補間を行う関数
    pub fn interpolate(&self, t: f64) -> f64 {
        let n = self.x.len() - 1;
        let idx = if t < self.x[1] {
            0
        } else if t >= self.x[n-1] {
            n-1
        } else {
            let mut i = n-1;
            for j in 1..n-1 {
                if self.x[j] <= t && t < self.x[j+1] {
                    i = j;
                    break;
                }
            }
            i
        };

        let h = t-self.x[idx];
        self.a[idx] + self.b[idx]*h + self.c[idx]*h.powi(2) + self.d[idx]*h.powi(3)
    }

    pub fn calc_spline(&mut self) {
        let (a,b,c,d) = self.calc_spline_coef();
        self.a = a;
        self.b = b;
        self.c = c;
        self.d = d;
    }

    // return a,b,c,d
    pub fn calc_spline_coef(&self) -> (DVector<f64>, DVector<f64>, DVector<f64>, DVector<f64>){
        let n = self.x.len();

        let h = DVector::from_fn(n-1, |i, _| self.x[i+1] - self.x[i]);

        // a_i = y_i (i=1,...,N+1)
        let a = DVector::from_fn(n, |i, _| self.y[i]);

        // k = [0, 
        //     -3(a_2 - a_1)/h_1 + 3(a_3 - a_2)/h_2 ,
        //     -3(a_3 - a_2)/h_2 + 3(a_4 - a_3)/h_3 ,
        //     ...
        //     -3(a_N-1 - a_N-2)/h_N-2 + 3(a_N - a_N-1)/h_N-1 ,
        //     -3(a_N - a_N-1)/h_N-1 + 3(a_N+1 - a_N)/h_N ,
        //     ,0] (N+1次)
        let k = DVector::from_fn(n, |i, _| {
            if i == 0 || i == n-1 {
                0.0
            } else {
                -3.0*(a[i]-a[i-1]) / h[i-1] + 3.0*(a[i+1]-a[i]) / h[i]
            }
        });

        // M = [
        // 1, 0, 0, ...
        // h_1, 2(h_1+h_2), h2, 0... ,0
        // 0, h_2, 2(h_2+h_3), h3, 0, ...,0
        // ...
        // 0, ... h_N-2, 2(h_N-2+h_N-1), h_N-1, 0
        // 0, ... 0, h_N-1, 2(h_N-1+h_N), h_N
        // 0, ...,0,   0 ,         0, 1
        // ] (N+1xN+1次)
        let mut m = DMatrix::<f64>::zeros(n,n);
        m[(0, 0)] = 1.0;
        m[(n-1,n-1)] = 1.0;
        for i in 1..n-1 {
            m[(i,i-1)] = h[i-1];
            m[(i,i)] = 2.0 * (h[i-1] + h[i]);
            m[(i,i+1)] = h[i];
        }

        // calc c
        let decomp = m.lu();
        let c = decomp.solve(&k).expect("Linear resolution failed");

        // d_i = (c_i+1 - c_i)/3h_i (N次)
        let d = DVector::from_fn(n-1, |i, _| (c[i+1] - c[i])/(3.0*h[i]));

        // b_i = (a_i+1 - a_i)/h_i - c_i h_i - d_i h_i^2 (N次)
        let b = DVector::from_fn(n-1, |i, _| {
            (a[i+1] - a[i])/h[i] - c[i]*h[i] - d[i]*h[i]*h[i]
        });

        (a,b,c,d)
    }
}

#[cfg(test)]
mod tests {
    use super::CubicSplineNa;
    use nalgebra::DVector;

    #[test]
    fn it_works() {
        let x = vec![0.0,1.0,4.0,5.0,8.0];
        let y = vec![0.0,3.0,4.0,1.0,2.0];

        let sp = CubicSplineNa::new(x,y);

        let (a,b,c,d) = sp.calc_spline_coef();
        assert_eq!(a, DVector::from_vec(vec![0.0, 3.0, 4.0, 1.0, 2.0]));
        assert_eq!(b, DVector::from_vec(vec![3.1805555555555554, 2.6388888888888893, -2.6527777777777777, -2.472222222222222]));
        assert_eq!(c, DVector::from_vec(vec![0.0, -0.5416666666666666, -1.2222222222222223, 1.4027777777777777, 0.0]));
        assert_eq!(d, DVector::from_vec(vec![-0.18055555555555555, -0.0756172839506173, 0.875, -0.1558641975308642]));
    }
}
