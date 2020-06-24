/// 3次スプライン補間を簡易実装したもの
mod cubic_spline_na;
mod cubic_spline;

pub use cubic_spline::CubicSpline;
pub use cubic_spline_na::CubicSplineNa;