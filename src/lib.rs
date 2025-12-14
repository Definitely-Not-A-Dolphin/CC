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
        let e_to_i_self = Complex::exp(self * Complex::i());
        (e_to_i_self - Complex::inv(e_to_i_self)) / (Complex::i() * T::two())
    }

    pub fn cos(self) -> Complex<T> {
        let e_to_i_self = Complex::exp(self * Complex::i());
        (e_to_i_self + Complex::inv(e_to_i_self)) / T::two()
    }
}

impl<T: Float> Add<T> for Complex<T> {
    type Output = Complex<T>;

    fn add(self, rhs: T) -> Complex<T> {
        Complex::new(self.real + rhs, self.imag)
    }
}

//impl<T: Float> Add<Complex<T>> for T {
//    type Output = Complex<T>;
//
//    fn add(self, rhs: Complex<T>) -> Complex<T> {
//        Complex::new(rhs.real + self, rhs.imag)
//    }
//}

impl<T: Float> Add<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn add(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(self.real + rhs.real, self.imag + rhs.imag)
    }
}

impl<T: Float> Sub<T> for Complex<T> {
    type Output = Complex<T>;

    fn sub(self, rhs: T) -> Complex<T> {
        Complex::new(self.real - rhs, self.imag)
    }
}

//impl<T: Float> Sub<Complex<T>> for T {
//    type Output = Complex<T>;
//
//    fn sub(self, rhs: Complex<T>) -> Complex<T> {
//        Complex::new(self - rhs.real, -rhs.imag)
//    }
//}

impl<T: Float> Sub<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn sub(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(self.real - rhs.real, self.imag + rhs.imag)
    }
}

impl<T: Float> Neg for Complex<T> {
    type Output = Complex<T>;

    fn neg(self) -> Complex<T> {
        Complex::new(-self.real, -self.imag)
    }
}

impl<T: Float> Mul<T> for Complex<T> {
    type Output = Complex<T>;

    fn mul(self, rhs: T) -> Complex<T> {
        Complex::new(self.real * rhs, self.imag * rhs)
    }
}

//impl<T: Float> Mul<Complex<T>> for T {
//    type Output = Complex<T>;
//
//    fn mul(self, rhs: Complex<T>) -> Complex<T> {
//        Complex::new(rhs.real * self, rhs.imag * self)
//    }
//}

impl<T: Float> Mul<Complex<T>> for Complex<T> {
    type Output = Complex<T>;

    fn mul(self, rhs: Complex<T>) -> Complex<T> {
        Complex::new(
            self.real * rhs.real - self.imag * rhs.imag,
            self.real * rhs.imag + self.imag * rhs.real,
        )
    }
}

impl<T: Float> Div<T> for Complex<T> {
    type Output = Complex<T>;

    fn div(self, rhs: T) -> Complex<T> {
        Complex::new(self.real / rhs, self.imag / rhs)
    }
}

//impl<T: Float> Div<Complex<T>> for T {
//    type Output = Complex<T>;
//
//    fn div(self, rhs: Complex<T>) -> Complex<T> {
//        rhs.inv() * self
//    }
//}

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
mod tests {
    use std::f64;

    use super::*;

    #[test]
    fn complex_unary_operators() {
        let z1 = Complex::new(3f64, 4f64);
        let z2 = Complex::new(5.2, -0.9);

        // real and imag
        assert_eq!(z1.real(), 3f64);
        assert_eq!(z1.imag(), 4f64);
        assert_eq!(z2.real(), 5.2);
        assert_eq!(z2.imag(), -0.9);

        // abs
        assert_eq!(z1.abs(), 5f64);
        assert_eq!(z2.abs(), f64::sqrt(27.85));

        // square abs
        assert_eq!(z1.square_abs(), 25f64);
        assert_eq!(z2.square_abs(), 27.85);

        // arg
        assert_eq!(z1.arg(), 0.9272952180016123);
        assert_eq!(z2.arg(), -0.1713791263895069);

        // inv
        assert_eq!(z1.inv(), Complex::new(0.12, -0.16));
        assert_eq!(
            z2.inv(),
            Complex::new(0.1867145421903052, 0.03231597845601436)
        );

        // exp
        assert_eq!(
            Complex::new(0f64, f64::consts::PI).exp(),
            Complex::new(-1f64, 1.2246467991473532e-16)
        );
        assert_eq!(
            z1.exp(),
            Complex::new(-13.128783081462158, -15.200784463067954)
        );

        // ln
        assert_eq!(
            Complex::new(0f64, 0f64).ln(),
            Complex::new(f64::NEG_INFINITY, 0f64)
        );
        assert_eq!(
            Complex::new(-1f64, 0f64).ln(),
            Complex::new(0f64, f64::consts::PI)
        );
    }

    #[test]
    fn simple_binary_operators() {
        let z1 = Complex::new(3f32, 4f32);
        let z2 = Complex::new(-2.5, 6.23);

        // add
        assert_eq!(z1 + 5f32, Complex::new(8f32, 4f32));
        //assert_eq!(5f64 + z1, Complex::new(8f64, 4f64));
        assert_eq!(z1 + z2, Complex::new(0.5, 10.23));

        // subtract
        assert_eq!(z1 - 5f32, Complex::new(-2f32, 4f32));
        //assert_eq!(5f64 - z1, Complex::new(2f64, -4f64));
        assert_eq!(z1 - z2, Complex::new(5.5, 10.23));

        // multiply
        assert_eq!(z1 * -0.5, Complex::new(-1.5, -2f32));
        //assert_eq!(1.8 * z2, Complex::new(-4.5, 11.214));
        assert_eq!(
            z1 * z2,
            Complex::new(-32.42, z1.real * z2.imag + z1.imag * z2.real)
        );

        // dividing
        assert_eq!(z1 / -0.5, Complex::new(-6f32, -8f32));
        // assert_eq!(
        //     1.8 / z2,
        //     Complex::new(-0.0998604173277796, -0.24885215998082677)
        // );
        assert_eq!(
            z1 / z2,
            Complex::new(
                0.38657075,
                z1.real * z2.inv().imag + z1.imag * z2.inv().real
            )
        );
    }

    #[test]
    fn complex_binary_operators() {
        let z1 = Complex::new(3f64, 4f64);
        let z2 = Complex::new(5.2, -0.9);

        // powi
        assert_eq!(z1.powi(6), Complex::new(11753f64, -10296f64));
        assert_eq!(
            z2.powi(-2),
            Complex::new(0.03381799780176568, 0.012067726245692974)
        );

        // powf
        assert_eq!(z1.powf(3f64), Complex::new(-117f64, 43.99999999999998));
        assert_eq!(
            z2.powf(-2.5),
            Complex::new(0.014217542838549325, 0.00649377309897787)
        );

        // powc
        assert_eq!(
            z1.powc(z2),
            Complex::new(-9667.467998399987, -2282.430226186542)
        );
        assert_eq!(
            z2.powc(z1),
            Complex::new(288.7067987011787, -41.7623644411436)
        );
    }
}
