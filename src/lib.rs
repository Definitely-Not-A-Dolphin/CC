/// CCMath: a crate for doing math with complex numbers
//#![warn(missing_docs, clippy::all)]
use num_traits::{self, Float};
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Struct representing a complex number
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Complex<T: Float> {
    real: T,
    imag: T,
}

pub type CC<T> = Complex<T>;

trait Numbers: Float {
    fn two() -> Self;
    fn ten() -> Self;
}

impl<T: Float> Numbers for T {
    fn two() -> T {
        T::one() + T::one()
    }
    fn ten() -> T {
        T::two() * (T::two() * T::two() + T::one())
    }
}

impl<T: Float> Complex<T> {
    /// Creates a new [`Complex`].
    pub fn new(real: T, imag: T) -> Complex<T> {
        Complex { real, imag }
    }

    pub fn i() -> Complex<T> {
        Complex::new(T::zero(), T::one())
    }

    /// Returns the real part of this [`Complex`].
    pub fn real(self) -> T {
        self.real
    }

    /// Returns the imaginary part of this [`Complex`].
    pub fn imag(self) -> T {
        self.imag
    }

    /// Returns the conjugate of this [`Complex`].
    pub fn conj(self) -> Complex<T> {
        Complex::new(self.real, -self.imag)
    }

    /// Returns the square of the absolute value of this [`Complex`].
    pub fn square_abs(self) -> T {
        self.real.powi(2) + self.imag.powi(2)
    }

    /// Returns the absolute value of this [`Complex`].
    pub fn abs(self) -> T {
        Complex::square_abs(self).sqrt()
    }

    /// Returns the arg on the interval (-PI, PI] of this [`Complex`].
    pub fn arg(self) -> T {
        if self == Complex::new(T::zero(), T::zero()) {
            return T::zero();
        };
        self.imag.signum() * T::acos(self.real / Complex::abs(self))
    }

    /// Returns the square root of this [`Complex`].
    pub fn sqrt(self) -> Complex<T> {
        let root_real = ((self.real + self.abs()) / T::two()).sqrt();
        let root_imag = ((-self.real + self.abs()) / T::two()).sqrt();
        Complex::new(root_real, self.imag.signum() * root_imag)
    }

    /// Returns the multiplicative inverse of this [`Complex`].
    pub fn inv(self) -> Complex<T> {
        Complex::conj(self) / Complex::square_abs(self)
    }

    /// Returns this [`Complex`] raised to a power using repeated multiplication.
    pub fn powi(self, exponent: i64) -> Complex<T> {
        match exponent {
            0 => Complex::new(T::one(), T::zero()),
            1 => self,
            -1 => self.inv(),
            _ => {
                let mut result = self;
                for _ in 2..=exponent.abs() {
                    result = result * self;
                }
                if exponent < 0 {
                    return result.inv();
                }
                result
            }
        }
    }

    /// Returns this [`Complex`] raised to a power using De Moivre's formula.
    pub fn powf(self, exponent: T) -> Complex<T> {
        Complex::new(T::cos(self.arg() * exponent), T::sin(self.arg() * exponent))
            * T::powf(Complex::abs(self), exponent)
    }

    /// Returns this [`Complex`] raised to a complex power.
    pub fn powc(self, exponent: Complex<T>) -> Complex<T> {
        Complex::powf(self, exponent.real)
            * Complex::exp(self.ln() * Complex::new(T::zero(), exponent.imag))
    }

    /// Returns e raised to the power of this [`Complex`].
    pub fn exp(self) -> Complex<T> {
        Complex::new(T::cos(self.imag), T::sin(self.imag)) * T::exp(self.real)
    }

    /// Returns base raised to the power of this [`Complex`].
    pub fn expf(self, base: T) -> Complex<T> {
        if base == T::zero() {
            return Complex::new(T::zero(), T::zero());
        };
        Complex::exp(self * base.ln())
    }

    /// Returns the natural logarithm of the absolute value of this [`Complex`].
    pub fn ln_abs(self) -> T {
        Complex::square_abs(self).ln() / (T::one() + T::one())
    }

    /// Returns the natural logarithm of this [`Complex`].
    pub fn ln(self) -> Complex<T> {
        Complex::new(Complex::ln_abs(self), Complex::arg(self))
    }

    /// Returns the logarithm base 10 of this [`Complex`].
    pub fn log(self) -> Complex<T> {
        Complex::ln(self) / T::ln(T::ten())
    }

    /// Returns the logarithm base n of this [`Complex`].
    pub fn logn(self, base: T) -> Complex<T> {
        Complex::ln(self) / T::ln(base)
    }

    pub fn sin(self) -> Complex<T> {
        Complex::new(
            T::sin(self.real) * T::cosh(self.imag),
            T::cos(self.real) * T::sinh(self.imag),
        )
    }

    pub fn cos(self) -> Complex<T> {
        Complex::new(
            T::cos(self.real) * T::cosh(self.imag),
            -T::sin(self.real) * T::sinh(self.imag),
        )
    }

    pub fn tan(self) -> Complex<T> {
        Complex::sin(self) / Complex::cos(self)
    }

    pub fn cot(self) -> Complex<T> {
        Complex::cos(self) / Complex::sin(self)
    }

    pub fn sec(self) -> Complex<T> {
        Complex::cos(self).inv()
    }

    pub fn csc(self) -> Complex<T> {
        Complex::sin(self).inv()
    }
}

/// Operator overloading: Complex<T> + T
impl<T: Float> Add<T> for Complex<T> {
    type Output = Complex<T>;

    fn add(self, rhs: T) -> Complex<T> {
        Complex::new(self.real + rhs, self.imag)
    }
}

/// Operator overloading: f32 + Complex<f32>
impl Add<Complex<f32>> for f32 {
    type Output = Complex<f32>;

    fn add(self, rhs: Complex<f32>) -> Complex<f32> {
        Complex::new(self + rhs.real, rhs.imag)
    }
}

/// Operator overloading: f64 + Complex<f64>
impl Add<Complex<f64>> for f64 {
    type Output = Complex<f64>;

    fn add(self, rhs: Complex<f64>) -> Complex<f64> {
        Complex::new(self + rhs.real, rhs.imag)
    }
}

/// Operator overloading: Complex<T> + Complex<T>
impl<T: Float> Add<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn add(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(self.real + rhs.real, self.imag + rhs.imag)
    }
}

/// Operator overloading: Complex<T> - T
impl<T: Float> Sub<T> for Complex<T> {
    type Output = Complex<T>;

    fn sub(self, rhs: T) -> Complex<T> {
        Complex::new(self.real - rhs, self.imag)
    }
}

/// Operator overloading: f32 - Complex<f32>
impl Sub<Complex<f32>> for f32 {
    type Output = Complex<f32>;

    fn sub(self, rhs: Complex<f32>) -> Complex<f32> {
        Complex::new(self - rhs.real, -rhs.imag)
    }
}

/// Operator overloading: f64 - Complex<f64>
impl Sub<Complex<f64>> for f64 {
    type Output = Complex<f64>;

    fn sub(self, rhs: Complex<f64>) -> Complex<f64> {
        Complex::new(self - rhs.real, -rhs.imag)
    }
}

/// Operator overloading: Complex<T> - Complex<T>
impl<T: Float> Sub<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn sub(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(self.real - rhs.real, self.imag + rhs.imag)
    }
}

/// Operator overloading: -Complex<T>
impl<T: Float> Neg for Complex<T> {
    type Output = Complex<T>;

    fn neg(self) -> Complex<T> {
        Complex::new(-self.real, -self.imag)
    }
}

/// Operator overloading: T * Complex<T>
impl<T: Float> Mul<T> for Complex<T> {
    type Output = Complex<T>;

    fn mul(self, rhs: T) -> Complex<T> {
        Complex::new(self.real * rhs, self.imag * rhs)
    }
}

/// Operator overloading: Complex<f32> * f32
impl Mul<Complex<f32>> for f32 {
    type Output = Complex<f32>;

    fn mul(self, rhs: Complex<f32>) -> Complex<f32> {
        Complex::new(rhs.real * self, rhs.imag * self)
    }
}

/// Operator overloading: Complex<f64> * f64
impl Mul<Complex<f64>> for f64 {
    type Output = Complex<f64>;

    fn mul(self, rhs: Complex<f64>) -> Complex<f64> {
        Complex::new(rhs.real * self, rhs.imag * self)
    }
}

/// Operator overloading: Complex<T> * Complex<T>
impl<T: Float> Mul<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn mul(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(
            self.real * rhs.real - self.imag * rhs.imag,
            self.real * rhs.imag + self.imag * rhs.real,
        )
    }
}

/// Operator overloading: Complex<T> / T
impl<T: Float> Div<T> for Complex<T> {
    type Output = Complex<T>;

    fn div(self, rhs: T) -> Complex<T> {
        Complex::new(self.real / rhs, self.imag / rhs)
    }
}

/// Operator overloading: f32 / Complex<f32>
impl Div<Complex<f32>> for f32 {
    type Output = Complex<f32>;

    fn div(self, rhs: Complex<f32>) -> Complex<f32> {
        Complex::conj(rhs) / Complex::square_abs(rhs) * self
    }
}

/// Operator overloading: f64 / Complex<f64>
impl Div<Complex<f64>> for f64 {
    type Output = Complex<f64>;

    fn div(self, rhs: Complex<f64>) -> Complex<f64> {
        Complex::conj(rhs) / Complex::square_abs(rhs) * self
    }
}

/// Operator overloading: Complex<T> / Complex<T>
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

#[cfg(test)]
mod tests;
