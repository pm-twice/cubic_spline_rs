use nalgebra::{DMatrix,DVector};
/// # 3次スプライン補間を行うクラス。  
/// Nalgebraを用いる。
/// 
/// ### 参考
/// [TajimaRobotics](https://tajimarobotics.com/cubic-spline-interpolation-2/)  
/// [TajimaRobotics](https://tajimarobotics.com/cubic-spline-interpolation-3/)
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
