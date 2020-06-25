/// # 3次スプライン補間を行うクラス。  
/// Vecのみで処理。
/// 
/// [fac.ksu.edu.sa](https://fac.ksu.edu.sa/sites/default/files/numerical_analysis_9th.pdf#page=167)
/// [その実装](https://gist.github.com/svdamani/1015c5c4b673c3297309)
/// 
pub struct CubicSpline {
    x: Vec<f64>,
    y: Vec<f64>,

    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}


impl CubicSpline {
    /// xは単調増加を仮定する
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> CubicSpline {
        let size = x.len();
        if size != y.len() {
            panic!("x & y must be same size");
        }
        CubicSpline {
            x: x,
            y: y,
            a: vec![],
            b: vec![],
            c: vec![], 
            d: vec![],
        }
    }

    pub fn interpolate_vec(&self, t: &Vec<f64>) -> Vec<f64> {
        t.iter().map(|v| self.interpolate(*v)).collect()
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

    /// Nalgebraを用いない場合
    /// https://gist.github.com/svdamani/1015c5c4b673c3297309
    /// based on https://fac.ksu.edu.sa/sites/default/files/numerical_analysis_9th.pdf#page=167
    pub fn calc_spline_coef(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>){
        let n = self.x.len() - 1;

        // step 0
        let a: Vec<f64> = (0..n+1).map(|i| self.y[i]).collect();
        let mut l = vec![0.0; n+1];
        let mut u = vec![0.0; n+1];
        let mut z = vec![0.0; n+1];
        let mut c = vec![0.0; n+1];
        let mut b = vec![0.0; n];
        let mut d = vec![0.0; n];

        // Step 1
        let h: Vec<f64> = (0..n).map(|i| self.x[i+1]-self.x[i]).collect();

        // Step 2
        let mut m = vec![0.0; n];
        for i in 1..n {
            m[i] = 3.*(a[i+1] - a[i]) / h[i] - 3.* (a[i] - a[i-1]) / h[i-1];
        }

        // Step 3
        l[0] = 1.0;
        u[0] = 0.0;
        z[0] = 0.0;

        // Step 4
        for i in 1..n {
            l[i] = 2.0 * (self.x[i + 1] - self.x[i - 1]) - h[i - 1] * u[i - 1];
            u[i] = h[i] / l[i];
            z[i] = (m[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        // Step 5
        l[n] = 1.0;
        z[n] = 0.0;
        c[n] = 0.0;

        // Step 6
        for j in (0..n).rev() {
            c[j] = z[j] - u[j] * c[j + 1];
            b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }

        (a,b,c,d)
    }
}

#[cfg(test)]
mod tests {
    use super::CubicSpline;

    #[test]
    fn it_works() {
        let x = vec![0.0,1.0,4.0,5.0,8.0];
        let y = vec![0.0,3.0,4.0,1.0,2.0];

        let sp = CubicSpline::new(x,y);

        let (a,b,c,d) = sp.calc_spline_coef();
        assert_eq!(a, vec![0.0, 3.0, 4.0, 1.0, 2.0]);
        assert_eq!(b, vec![3.1805555555555554, 2.638888888888889, -2.6527777777777777, -2.472222222222222]);
        assert_eq!(c, vec![0.0, -0.5416666666666667, -1.222222222222222, 1.4027777777777777, 0.0]);
        assert_eq!(d, vec![-0.18055555555555558, -0.07561728395061726, 0.875, -0.1558641975308642]);
    }
}
