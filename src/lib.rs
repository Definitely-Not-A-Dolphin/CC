/// CCMath: a crate for doing math with complex numbers
use num_traits::{self, Float};

/// Struct representing a complex number
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Complex<T: Float> {
    real: T,
    imag: T,
}

/// Alias for [`Complex`]
pub type CC<T> = Complex<T>;
// Personal tip: use this alias when code is getting a little hard to read, it cleans things up!

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
    pub fn new(real: T, imag: T) -> Self {
        Self { real, imag }
    }

    pub fn i() -> Self {
        Self::new(T::zero(), T::one())
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
    pub fn conj(self) -> Self {
        Self::new(self.real, -self.imag)
    }

    /// Returns the square of the absolute value of this [`Complex`].
    pub fn square_abs(self) -> T {
        self.real.powi(2) + self.imag.powi(2)
    }

    /// Returns the absolute value of this [`Complex`].
    pub fn abs(self) -> T {
        Self::square_abs(self).sqrt()
    }

    /// Returns the arg on the interval (-PI, PI] of this [`Complex`].
    pub fn arg(self) -> T {
        if self == Self::new(T::zero(), T::zero()) {
            return T::zero();
        };
        self.imag.signum() * T::acos(self.real / Self::abs(self))
    }

    /// Returns the square root of this [`Complex`].
    pub fn sqrt(self) -> Self {
        let root_real = ((self.real + self.abs()) / T::two()).sqrt();
        let root_imag = ((-self.real + self.abs()) / T::two()).sqrt();
        Self::new(root_real, self.imag.signum() * root_imag)
    }

    /// Returns the multiplicative inverse of this [`Complex`].
    pub fn inv(self) -> Self {
        Self::conj(self) / Self::square_abs(self)
    }

    /// Returns this [`Complex`] raised to a power using repeated multiplication.
    pub fn powi(self, exponent: i64) -> Self {
        match exponent {
            0 => Self::new(T::one(), T::zero()),
            1 => self,
            -1 => self.inv(),
            _ => {
                let mut result = self;
                for _ in 2..=exponent.abs() {
                    result *=  self;
                }
                if exponent < 0 {
                    return result.inv();
                }
                result
            }
        }
    }

    /// Returns this [`Complex`] raised to a power using De Moivre's formula.
    pub fn powf(self, exponent: T) -> Self {
        Self::new(T::cos(self.arg() * exponent), T::sin(self.arg() * exponent))
            * T::powf(Self::abs(self), exponent)
    }

    /// Returns this [`Complex`] raised to a complex power.
    pub fn powc(self, exponent: Self) -> Self {
        Self::powf(self, exponent.real) * Self::exp(self.ln() * Self::new(T::zero(), exponent.imag))
    }

    /// Returns e raised to the power of this [`Complex`].
    pub fn exp(self) -> Self {
        Self::new(T::cos(self.imag), T::sin(self.imag)) * T::exp(self.real)
    }

    /// Returns base raised to the power of this [`Complex`].
    pub fn expf(self, base: T) -> Self {
        if base == T::zero() {
            return Self::new(T::zero(), T::zero());
        };
        Self::exp(self * base.ln())
    }

    /// Returns the natural logarithm of the absolute value of this [`Complex`].
    pub fn ln_abs(self) -> T {
        T::ln(Self::square_abs(self)) / (T::one() + T::one())
    }

    /// Returns the natural logarithm of this [`Complex`].
    pub fn ln(self) -> Self {
        Self::new(Self::ln_abs(self), Self::arg(self))
    }

    /// Returns the logarithm base 10 of this [`Complex`].
    pub fn log(self) -> Self {
        Self::ln(self) / T::ln(T::ten())
    }

    /// Returns the logarithm base n of this [`Complex`].
    pub fn logn(self, base: T) -> Self {
        Self::ln(self) / T::ln(base)
    }

    /// Returns the sine of this [`Complex`].
    pub fn sin(self) -> Self {
        Self::new(
            T::sin(self.real) * T::cosh(self.imag),
            T::cos(self.real) * T::sinh(self.imag),
        )
    }

    pub fn arcsin(self) -> Self {
        -Self::i() * Self::ln(Self::sqrt(-self.powi(2) + T::one()) + self * Self::i())
    }

    /// Returns the cosine of this [`Complex`].
    pub fn cos(self) -> Self {
        Self::new(
            T::cos(self.real) * T::cosh(self.imag),
            -T::sin(self.real) * T::sinh(self.imag),
        )
    }

    /// Returns the tangent of this [`Complex`].
    pub fn tan(self) -> Self {
        Self::sin(self) / Self::cos(self)
    }

    /// Returns the cotangent of this [`Complex`].
    pub fn cot(self) -> Self {
        Self::cos(self) / Self::sin(self)
    }

    /// Returns the secant of this [`Complex`].
    pub fn sec(self) -> Self {
        Self::cos(self).inv()
    }

    /// Returns the cosecant of this [`Complex`].
    pub fn csc(self) -> Self {
        Self::sin(self).inv()
    }
}

mod overloading;

#[cfg(test)]
mod tests;
