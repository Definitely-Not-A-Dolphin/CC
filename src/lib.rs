use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Complex {
    real: f64,
    imag: f64,
}

fn expimag(exponent: f64) -> Complex {
    Complex::new(f64::cos(exponent), f64::sin(exponent))
}

impl Complex {
    /// Creates a new [`Complex`].
    pub fn new(real: f64, imag: f64) -> Complex {
        Complex {
            real: real,
            imag: imag,
        }
    }

    /// Returns the real part of this [`Complex`].
    pub fn real(self) -> f64 {
        self.real
    }

    /// Returns the imaginary part of this [`Complex`].
    pub fn imag(self) -> f64 {
        self.imag
    }

    /// Returns the conjugate of this [`Complex`].
    pub fn conj(self) -> Complex {
        Complex::new(self.real, -self.imag)
    }

    /// Returns the square of the absolute value of this [`Complex`].
    pub fn square_abs(self) -> f64 {
        self.real.powi(2) + self.imag.powi(2)
    }

    /// Returns the absolute value of this [`Complex`].
    pub fn abs(self) -> f64 {
        self.square_abs().sqrt()
    }

    /// Returns the arg on the interval [0, 2PI) of this [`Complex`].
    pub fn arg(self) -> f64 {
        if self == Complex::new(0f64, 0f64) {
            return 0f64;
        };
        self.imag.signum() * f64::acos(self.real / self.abs())
    }

    /// Returns the square root of this [`Complex`].
    pub fn sqrt(self) -> Complex {
        let root_real = ((self.real + self.abs()) / 2f64).sqrt();
        let root_imag = ((-self.real + self.abs()) / 2f64).sqrt();
        return Complex::new(root_real, self.imag.signum() * root_imag);
    }

    /// Returns the multiplicative inverse of this [`Complex`].
    pub fn inv(self) -> Complex {
        Complex::conj(self) / Complex::square_abs(self)
    }

    /// Returns this [`Complex`] raised to a power using repeated multiplication.
    pub fn powi(self, exponent: i64) -> Complex {
        match exponent {
            0 => Complex::new(1f64, 0f64),
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
                return result;
            }
        }
    }

    /// Returns this [`Complex`] raised to a power using De Moivre's formula.
    pub fn powf(self, exponent: f64) -> Complex {
        Complex::abs(self).powf(exponent)
            * Complex::new(
                f64::cos(exponent * self.arg()),
                f64::sin(exponent * self.arg()),
            )
    }

    /// Returns this [`Complex`] raised to a complex power.
    pub fn powc(self, exponent: Complex) -> Complex {
        Complex::powf(self, exponent.real)
            * Complex::exp(self.ln() * Complex::new(0f64, exponent.imag))
    }

    /// Returns e raised to the power of this [`Complex`].
    pub fn exp(self) -> Complex {
        f64::exp(self.real) * expimag(self.imag)
    }

    /// Returns base raised to the power of this [`Complex`].
    pub fn expf(self, base: f64) -> Complex {
        if base == 0f64 {
            return Complex::new(0f64, 0f64);
        };
        Complex::exp(base.ln() * self)
    }

    /// Returns the natural logarithm of the absolute value of this [`Complex`].
    pub fn ln_abs(self) -> f64 {
        Complex::square_abs(self).ln() / 2f64
    }

    /// Returns the natural logarithm of this [`Complex`].
    pub fn ln(self) -> Complex {
        Complex::new(Complex::ln_abs(self), Complex::arg(self))
    }

    /// Returns the logarithm base 10 of this [`Complex`].
    pub fn log(self) -> Complex {
        Complex::ln(self) / f64::ln(10f64)
    }

    /// Returns the logarithm base n of this [`Complex`].
    pub fn logn(self, base: f64) -> Complex {
        Complex::ln(self) / f64::ln(base)
    }
}

impl Add<f64> for Complex {
    type Output = Complex;

    fn add(self, rhs: f64) -> Complex {
        Complex::new(self.real + rhs, self.imag)
    }
}

impl Add<Complex> for f64 {
    type Output = Complex;

    fn add(self, rhs: Complex) -> Complex {
        Complex::new(rhs.real + self, rhs.imag)
    }
}

impl Add<Complex> for Complex {
    type Output = Complex;

    fn add(self, rhs: Complex) -> Complex {
        Complex::new(self.real + rhs.real, self.imag + rhs.imag)
    }
}

impl Sub<f64> for Complex {
    type Output = Complex;

    fn sub(self, rhs: f64) -> Complex {
        Complex::new(self.real - rhs, self.imag)
    }
}

impl Sub<Complex> for f64 {
    type Output = Complex;

    fn sub(self, rhs: Complex) -> Complex {
        Complex::new(self - rhs.real, -rhs.imag)
    }
}

impl Sub<Complex> for Complex {
    type Output = Complex;

    fn sub(self, rhs: Complex) -> Complex {
        Complex::new(self.real - rhs.real, self.imag + rhs.imag)
    }
}

impl Neg for Complex {
    type Output = Complex;

    fn neg(self) -> Complex {
        Complex::new(-self.real, -self.imag)
    }
}

impl Mul<f64> for Complex {
    type Output = Complex;

    fn mul(self, rhs: f64) -> Complex {
        Complex::new(self.real * rhs, self.imag * rhs)
    }
}

impl Mul<Complex> for f64 {
    type Output = Complex;

    fn mul(self, rhs: Complex) -> Complex {
        Complex::new(rhs.real * self, rhs.imag * self)
    }
}

impl Mul<Complex> for Complex {
    type Output = Complex;

    fn mul(self, rhs: Complex) -> Complex {
        Complex::new(
            self.real * rhs.real - self.imag * rhs.imag,
            self.real * rhs.imag + self.imag * rhs.real,
        )
    }
}

impl Div<f64> for Complex {
    type Output = Complex;

    fn div(self, rhs: f64) -> Complex {
        Complex::new(self.real / rhs, self.imag / rhs)
    }
}

impl Div<Complex> for f64 {
    type Output = Complex;

    fn div(self, rhs: Complex) -> Complex {
        self * rhs.inv()
    }
}

impl Div<Complex> for Complex {
    type Output = Complex;

    fn div(self, rhs: Complex) -> Complex {
        self * rhs.inv()
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
        let z1 = Complex::new(3f64, 4f64);
        let z2 = Complex::new(-2.5, 6.23);

        // add
        assert_eq!(z1 + 5f64, Complex::new(8f64, 4f64));
        assert_eq!(5f64 + z1, Complex::new(8f64, 4f64));
        assert_eq!(z1 + z2, Complex::new(0.5, 10.23));

        // subtract
        assert_eq!(z1 - 5f64, Complex::new(-2f64, 4f64));
        assert_eq!(5f64 - z1, Complex::new(2f64, -4f64));
        assert_eq!(z1 - z2, Complex::new(5.5, 10.23));

        // multiply
        assert_eq!(z1 * -0.5, Complex::new(-1.5, -2f64));
        assert_eq!(1.8 * z2, Complex::new(-4.5, 11.214));
        assert_eq!(z1 * z2, Complex::new(-32.42, 8.690000000000001));

        // dividing
        assert_eq!(z1 / -0.5, Complex::new(-6f64, -8f64));
        assert_eq!(
            1.8 / z2,
            Complex::new(-0.0998604173277796, -0.24885215998082677)
        );
        assert_eq!(
            z1 / z2,
            Complex::new(0.38657077107776017, -0.6366656384742215)
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
