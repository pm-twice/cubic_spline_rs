/// 3次スプライン補間を簡易実装したもの
mod cubic_spline_na;
mod cubic_spline;
mod linear_interpolation;
mod bspline;

pub use cubic_spline::CubicSpline;
pub use cubic_spline_na::CubicSplineNa;
pub use linear_interpolation::LinearInterpolation;
pub use bspline::BSpline;