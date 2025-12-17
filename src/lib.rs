/// CCMath: a crate for doing math with complex numbers
use num_traits::Float;
use std::fmt::{self, Debug, Display};

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
    fn half() -> Self;
}

impl<T: Float> Numbers for T {
    fn two() -> T {
        T::one() + T::one()
    }
    fn ten() -> T {
        T::two() * (T::two() * T::two() + T::one())
    }
    fn half() -> T {
        T::one() / T::two()
    }
}

impl<T: Float> Complex<T> {
    /// Creates a new [`Complex`].
    pub fn new(real: T, imag: T) -> Self {
        Self { real, imag }
    }

    /// Returns the imaginary number i
    pub fn i() -> Self {
        Self::new(T::zero(), T::one())
    }

    // Returns two i
    fn two_i() -> Self {
        Self::new(T::zero(), T::two())
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
            T::zero()
        } else {
            self.imag.signum() * T::acos(self.real / Self::abs(self))
        }
    }

    /// Returns the square root of this [`Complex`].
    pub fn sqrt(self) -> Self {
        Self::new(
            ((self.real + self.abs()) / T::two()).sqrt(),
            self.imag.signum() * ((-self.real + self.abs()) / T::two()).sqrt(),
        )
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
                    result *= self;
                }
                if exponent < 0 { result.inv() } else { result }
            }
        }
    }
    // Todo: Implement exponentiating by squaring

    /// Returns this [`Complex`] raised to a power using De Moivre's formula.
    pub fn powf(self, exponent: T) -> Self {
        let arg_exponent = self.arg() * exponent;
        Self::new(T::cos(arg_exponent), T::sin(arg_exponent)) * T::powf(Self::abs(self), exponent)
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
            Self::new(T::zero(), T::zero())
        } else {
            Self::exp(self * T::ln(base))
        }
    }

    /// Returns the natural logarithm of the absolute value of this [`Complex`].
    pub fn ln_abs(self) -> T {
        T::ln(Self::square_abs(self)) / T::two()
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
}

// Trig
impl<T: Float + Debug> Complex<T> {
    /// Returns the sine of this [`Complex`].
    pub fn sin(self) -> Self {
        Self::new(
            T::sin(self.real) * T::cosh(self.imag),
            T::cos(self.real) * T::sinh(self.imag),
        )
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

    // Inverse trig

    /// Returns the arcsine of this [`Complex`].
    pub fn arcsin(self) -> Self {
        -Self::i() * Self::ln(Self::sqrt(-self.powi(2) + T::one()) + Self::i() * self)
    }

    /// Returns the arccosine of this [`Complex`].
    pub fn arccos(self) -> Self {
        Self::i() * Self::ln(Self::sqrt(-self.powi(2) + T::one()) / Self::i() + self)
    }

    /// Returns the arctangent of this [`Complex`].
    pub fn arctan(self) -> Self {
        Self::two_i().inv()
            * Self::ln((Self::i() * self + T::one()) / (-Self::i() * self + T::one()))
    }

    /// Returns the arccotangent of this [`Complex`].
    pub fn arccot(self) -> Self {
        Self::arctan(Self::inv(self))
    }

    /// Returns the arcsecant of this [`Complex`].
    pub fn arcsec(self) -> Self {
        Self::arccos(Self::inv(self))
    }

    // Returns the arccosecant of this [`Complex`].
    pub fn arccsc(self) -> Self {
        Self::arcsin(Self::inv(self))
    }

    // Hyperbolic trig

    /// Returns the hyperbolic sine of this [`Complex`].
    pub fn sinh(self) -> Self {
        Self::new(
            T::sinh(self.real) * T::cos(self.imag),
            T::cosh(self.real) * T::sin(self.imag),
        )
    }

    /// Returns the hyperbolic cosine of this [`Complex`].
    pub fn cosh(self) -> Self {
        Self::new(
            T::cosh(self.real) * T::cos(self.imag),
            T::sinh(self.real) * T::sin(self.imag),
        )
    }

    /// Returns the hyperbolic tangent of this [`Complex`].
    pub fn tanh(self) -> Self {
        Self::sinh(self) / Self::cosh(self)
    }

    /// Returns the hyperbolic cotangent of this [`Complex`].
    pub fn coth(self) -> Self {
        Self::cosh(self) / Self::sinh(self)
    }

    /// Returns the hyperbolic secant of this [`Complex`].
    pub fn sech(self) -> Self {
        Self::cosh(self).inv()
    }

    /// Returns the hyperbolic cosecant of this [`Complex`].
    pub fn csch(self) -> Self {
        Self::sinh(self).inv()
    }

    // Inverse hyperbolic trig

    /// Returns the hyperbolic arcsine of this [`Complex`].
    pub fn arcsinh(self) -> Self {
        Self::ln(Self::sqrt(self.powi(2) + T::one()) + self)
    }

    /// Returns the hyperbolic arccosine of this [`Complex`].
    pub fn arccosh(self) -> Self {
        Self::ln(Self::sqrt(self.powi(2) - T::one()) + self)
    }

    /// Returns the hyperbolic arctangent of this [`Complex`].
    pub fn arctanh(self) -> Self {
        Self::ln((self + T::one()) / (-self + T::one())) * T::half()
    }

    /// Returns the hyperbolic arccotangent of this [`Complex`].
    pub fn arccoth(self) -> Self {
        Self::arctanh(Self::inv(self))
    }

    /// Returns the hyperbolic arcsecant of this [`Complex`].
    pub fn arcsech(self) -> Self {
        Self::arccosh(Self::inv(self))
    }

    /// Returns the hyperbolic arccosecant of this [`Complex`].
    pub fn arccsch(self) -> Self {
        Self::arcsinh(Self::inv(self))
    }
}

// Implements display
impl<T: Float + Display> fmt::Display for Complex<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // The `f` value implements the `Write` trait, which is what the
        // write! macro is expecting. Note that this formatting ignores the
        // various flags provided to format strings.
        if self.imag >= T::zero() {
            write!(f, "{} + {}i", self.real, self.imag)
        } else {
            write!(f, "{} - {}i", self.real, self.imag)
        }
    }
}

mod overloading;

#[cfg(test)]
mod tests;
