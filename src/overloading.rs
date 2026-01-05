use crate::Complex;
use num_traits::{self, Float};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

// Addition

/// Complex<T> + T
impl<T: Float> Add<T> for Complex<T> {
    type Output = Complex<T>;

    fn add(self, rhs: T) -> Complex<T> {
        Complex::new(self.real + rhs, self.imag)
    }
}

/// Complex<T> += T
impl<T: Float> AddAssign<T> for Complex<T> {
    fn add_assign(&mut self, rhs: T) {
        *self = *self + rhs;
    }
}

/// f32 + Complex<f32>
impl Add<Complex<f32>> for f32 {
    type Output = Complex<f32>;

    fn add(self, rhs: Complex<f32>) -> Complex<f32> {
        Complex::new(self + rhs.real, rhs.imag)
    }
}

/// f64 + Complex<f64>
impl Add<Complex<f64>> for f64 {
    type Output = Complex<f64>;

    fn add(self, rhs: Complex<f64>) -> Complex<f64> {
        Complex::new(self + rhs.real, rhs.imag)
    }
}

/// Complex<T> + Complex<T>
impl<T: Float> Add<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn add(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(self.real + rhs.real, self.imag + rhs.imag)
    }
}

/// Complex<T> += Complex<T>
impl<T: Float> AddAssign<Complex<T>> for Complex<T> {
    fn add_assign(&mut self, rhs: Complex<T>) {
        *self = *self + rhs;
    }
}

// Subtraction

/// Complex<T> - T
impl<T: Float> Sub<T> for Complex<T> {
    type Output = Complex<T>;

    fn sub(self, rhs: T) -> Complex<T> {
        Complex::new(self.real - rhs, self.imag)
    }
}

/// Complex<T> -= T
impl<T: Float> SubAssign<T> for Complex<T> {
    fn sub_assign(&mut self, rhs: T) {
        *self = *self - rhs
    }
}

/// f32 - Complex<f32>
impl Sub<Complex<f32>> for f32 {
    type Output = Complex<f32>;

    fn sub(self, rhs: Complex<f32>) -> Complex<f32> {
        Complex::new(self - rhs.real, -rhs.imag)
    }
}

/// f64 - Complex<f64>
impl Sub<Complex<f64>> for f64 {
    type Output = Complex<f64>;

    fn sub(self, rhs: Complex<f64>) -> Complex<f64> {
        Complex::new(self - rhs.real, -rhs.imag)
    }
}

/// Complex<T> - Complex<T>
impl<T: Float> Sub<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn sub(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(self.real - rhs.real, self.imag + rhs.imag)
    }
}

/// Complex<T> -= Complex<T>
impl<T: Float> SubAssign<Complex<T>> for Complex<T> {
    fn sub_assign(&mut self, rhs: Complex<T>) {
        *self = *self - rhs;
    }
}

// Negation

/// -Complex<T>
impl<T: Float> Neg for Complex<T> {
    type Output = Complex<T>;

    fn neg(self) -> Complex<T> {
        Complex::new(-self.real, -self.imag)
    }
}

// Multiplication

/// Complex<T> * T
impl<T: Float> Mul<T> for Complex<T> {
    type Output = Complex<T>;

    fn mul(self, rhs: T) -> Complex<T> {
        Complex::new(self.real * rhs, self.imag * rhs)
    }
}

/// Complex<T> *= T
impl<T: Float> MulAssign<T> for Complex<T> {
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

/// f32 * Complex<f32>
impl Mul<Complex<f32>> for f32 {
    type Output = Complex<f32>;

    fn mul(self, rhs: Complex<f32>) -> Complex<f32> {
        Complex::new(rhs.real * self, rhs.imag * self)
    }
}

/// f64 * Complex<f64>
impl Mul<Complex<f64>> for f64 {
    type Output = Complex<f64>;

    fn mul(self, rhs: Complex<f64>) -> Complex<f64> {
        Complex::new(rhs.real * self, rhs.imag * self)
    }
}

/// Complex<T> * Complex<T>
impl<T: Float> Mul<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn mul(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(
            self.real * rhs.real - self.imag * rhs.imag,
            self.real * rhs.imag + self.imag * rhs.real,
        )
    }
}

/// Complex<T> *= Complex<T>
impl<T: Float> MulAssign<Complex<T>> for Complex<T> {
    fn mul_assign(&mut self, rhs: Complex<T>) {
        *self = *self * rhs
    }
}

// Division

/// Complex<T> / T
impl<T: Float> Div<T> for Complex<T> {
    type Output = Complex<T>;

    fn div(self, rhs: T) -> Complex<T> {
        Complex::new(self.real / rhs, self.imag / rhs)
    }
}

/// Complex<T> /= T
impl<T: Float> DivAssign<T> for Complex<T> {
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
    }
}

/// f32 / Complex<f32>
impl Div<Complex<f32>> for f32 {
    type Output = Complex<f32>;

    fn div(self, rhs: Complex<f32>) -> Complex<f32> {
        Complex::conj(rhs) / Complex::square_abs(rhs) * self
    }
}

/// f64 / Complex<f64>
impl Div<Complex<f64>> for f64 {
    type Output = Complex<f64>;

    fn div(self, rhs: Complex<f64>) -> Complex<f64> {
        Complex::conj(rhs) / Complex::square_abs(rhs) * self
    }
}

/// Complex<T> / Complex<T>
impl<T: Float> Div<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn div(self, rhs: Complex<T>) -> Complex<T> {
        let rhsinv = rhs.inv();
        Complex::new(
            self.real * rhsinv.real - self.imag * rhsinv.imag,
            self.real * rhsinv.imag + self.imag * rhsinv.real,
        )
    }
}

/// Complex<T> /= Complex<T>
impl<T: Float> DivAssign<Complex<T>> for Complex<T> {
    fn div_assign(&mut self, rhs: Complex<T>) {
        *self = *self / rhs;
    }
}
