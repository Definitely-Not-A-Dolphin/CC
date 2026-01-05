use crate::ComplexPolar;
use num_traits::{self, Float};
use std::ops::{Div, DivAssign, Mul, MulAssign, Neg};

/// -ComplexPolar<T>
impl<T: Float> Neg for ComplexPolar<T> {
    type Output = ComplexPolar<T>;

    fn neg(self) -> ComplexPolar<T> {
        ComplexPolar::new(-self.radius, self.angle)
    }
}

// Multiplication

/// ComplexPolar<T> * T
impl<T: Float> Mul<T> for ComplexPolar<T> {
    type Output = ComplexPolar<T>;

    fn mul(self, rhs: T) -> ComplexPolar<T> {
        ComplexPolar::new(self.radius * rhs, self.angle)
    }
}

/// Complex<T> *= T
impl<T: Float> MulAssign<T> for ComplexPolar<T> {
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

/// f32 * Complex<f32>
impl Mul<ComplexPolar<f32>> for f32 {
    type Output = ComplexPolar<f32>;

    fn mul(self, rhs: ComplexPolar<f32>) -> ComplexPolar<f32> {
        ComplexPolar::new(rhs.radius * self, rhs.angle)
    }
}

/// f64 * ComplexPolar<f64>
impl Mul<ComplexPolar<f64>> for f64 {
    type Output = ComplexPolar<f64>;

    fn mul(self, rhs: ComplexPolar<f64>) -> ComplexPolar<f64> {
        ComplexPolar::new(rhs.radius * self, rhs.angle)
    }
}

/// ComplexPolar<T> * ComplexPolar<T>
impl<T: Float> Mul<ComplexPolar<T>> for ComplexPolar<T> {
    type Output = ComplexPolar<T>;

    fn mul(self, rhs: ComplexPolar<T>) -> ComplexPolar<T> {
        ComplexPolar::new(self.radius * rhs.radius, self.angle + rhs.angle)
    }
}

/// ComplexPolar<T> *= ComplexPolar<T>
impl<T: Float> MulAssign<ComplexPolar<T>> for ComplexPolar<T> {
    fn mul_assign(&mut self, rhs: ComplexPolar<T>) {
        *self = *self * rhs
    }
}

// Division

/// ComplexPolar<T> / T
impl<T: Float> Div<T> for ComplexPolar<T> {
    type Output = ComplexPolar<T>;

    fn div(self, rhs: T) -> ComplexPolar<T> {
        ComplexPolar::new(self.radius / rhs, self.angle)
    }
}

/// Complex<T> /= T
impl<T: Float> DivAssign<T> for ComplexPolar<T> {
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
    }
}

/// f32 / Complex<f32>
impl Div<ComplexPolar<f32>> for f32 {
    type Output = ComplexPolar<f32>;

    fn div(self, rhs: ComplexPolar<f32>) -> ComplexPolar<f32> {
        ComplexPolar::new(rhs.radius / self, rhs.angle)
    }
}

/// f64 / ComplexPolar<f64>
impl Div<ComplexPolar<f64>> for f64 {
    type Output = ComplexPolar<f64>;

    fn div(self, rhs: ComplexPolar<f64>) -> ComplexPolar<f64> {
        ComplexPolar::new(rhs.radius / self, rhs.angle)
    }
}

/// ComplexPolar<T> / ComplexPolar<T>
impl<T: Float> Div<ComplexPolar<T>> for ComplexPolar<T> {
    type Output = ComplexPolar<T>;

    fn div(self, rhs: ComplexPolar<T>) -> ComplexPolar<T> {
        ComplexPolar::new(self.radius / rhs.radius, self.angle - rhs.angle)
    }
}

/// ComplexPolar<T> /= ComplexPolar<T>
impl<T: Float> DivAssign<ComplexPolar<T>> for ComplexPolar<T> {
    fn div_assign(&mut self, rhs: ComplexPolar<T>) {
        *self = *self / rhs
    }
}
