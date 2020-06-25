/// # 線形補間を行うクラス
/// 
pub struct LinearInterpolation {
    x: Vec<f64>,
    y: Vec<f64>,
}


impl LinearInterpolation {
    /// xは単調増加を仮定する
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> LinearInterpolation {
        let size = x.len();
        if size != y.len() {
            panic!("x & y must be same size");
        }
        LinearInterpolation {
            x: x,
            y: y,
        }
    }

    pub fn interpolate_vec(&self, t: &Vec<f64>) -> Vec<f64> {
        t.iter().map(|v| self.interpolate(*v)).collect()
    }

    /// 与えられたtに対して、計算された係数から線形補間を行う関数
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

        // (f(t)-yi)/(t-xi) = (yi+1-yi)/(xi+1-xi)
        // f(t) = yi + (t-xi)(yi+1-yi)/(xi+1-xi)

        self.y[idx] + (t-self.x[idx])*(self.y[idx+1]-self.y[idx])/(self.x[idx+1]-self.x[idx])
    }
}

#[cfg(test)]
mod tests {
    use super::LinearInterpolation;

    #[test]
    fn it_works() {
        let x = vec![0.0,1.0,4.0,5.0,8.0];
        let y = vec![0.0,3.0,4.0,1.0,2.0];
        let lin = LinearInterpolation::new(x,y);

        assert_eq!(1.5, lin.interpolate(0.5));
    }
}
