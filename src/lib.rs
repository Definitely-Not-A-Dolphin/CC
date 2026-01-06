/// CCMath: a crate for doing math with complex numbers
use num_traits::Float;
use std::fmt::{Display, Formatter, Result};

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
    fn num(n: u32) -> Self;
    fn pi() -> Self;
    fn powc(self, exponent: Complex<Self>) -> Complex<Self>;
    fn powcp(self, exponent: ComplexPolar<Self>) -> ComplexPolar<Self>;
}

impl<T: Float> Numbers for T {
    /// Returns any whole positive number as T
    fn num(n: u32) -> T {
        (0..n).fold(T::zero(), |sum, _| sum + T::one())
    }
    fn pi() -> T {
        [
            1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3,
        ]
        .iter()
        .enumerate()
        .fold(T::num(3), |pi, (k, item)| {
            pi + T::num(*item) / T::num(10).powi(k as i32 + 1)
        })
    }
    /// Returns base raised to the power of this [`Complex`].
    fn powc(self, exponent: Complex<T>) -> Complex<T> {
        Complex::exp(exponent * T::ln(self))
    }
    /// Returns base raised to the power of this [`ComplexPolar`].
    fn powcp(self, exponent: ComplexPolar<T>) -> ComplexPolar<T> {
        Complex::exp(exponent.unpolarize() * T::ln(self)).polarize()
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

    /// Returns the real part of this [`Complex`].
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmath::Complex;
    ///
    /// let z = Complex::new(-1.4, 21.6);
    /// assert_eq!(z.real(), -1.4);
    /// ```
    pub fn real(self) -> T {
        self.real
    }

    /// Returns the imaginary part of this [`Complex`].
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmath::Complex;
    ///
    /// let z = Complex::new(-1.4, 21.6);
    /// assert_eq!(z.real(), -1.4);
    /// ```
    pub fn imag(self) -> T {
        self.imag
    }

    /// Returns the conjugate of this [`Complex`].
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmath::Complex;
    ///
    /// let z = Complex::new(4.0, 5.0);
    /// let z_conjugate = z.conj();
    ///
    /// assert_eq!(z_conjugate, Complex::new(4.0, -5.0));
    /// ```
    pub fn conj(self) -> Self {
        Self::new(self.real, -self.imag)
    }

    /// Returns the square of the absolute value of this [`Complex`].
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmath::Complex;
    ///
    /// let z1 = Complex::new(3.0, 4.0);
    /// let z2 = Complex::new(4.2, 2.1);
    ///
    /// assert_eq!(Complex::square_abs(z1), 25.0);
    /// assert_eq!(Complex::square_abs(z2), 22.05);
    /// ```
    pub fn square_abs(self) -> T {
        self.real.powi(2) + self.imag.powi(2)
    }

    /// Returns the absolute value of this [`Complex`].
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmath::Complex;
    ///
    /// let z1 = Complex::new(3.0, 4.0);
    /// let z2 = Complex::new(4.2, 2.1);
    ///
    /// assert_eq!(Complex::abs(z1), 5.0);
    /// assert_eq!(Complex::abs(z2), f32::sqrt(22.05));
    /// ```
    pub fn abs(self) -> T {
        T::sqrt(Self::square_abs(self))
    }

    /// Returns the argument on the interval (-PI, PI] of this [`Complex`].
    pub fn arg(self) -> T {
        if Self::abs(self) == T::zero() {
            T::zero()
        } else {
            self.imag.signum() * T::acos(self.real / Self::abs(self))
        }
    }

    /// Returns this [`Complex`] in polar form as a [`ComplexPolar`]
    pub fn polarize(self) -> ComplexPolar<T> {
        ComplexPolar::new(self.abs(), self.arg())
    }

    /// Returns the square root of this [`Complex`].
    pub fn sqrt(self) -> Self {
        Self::new(
            T::sqrt((self.real + Self::abs(self)) / T::num(2)),
            self.imag.signum() * T::sqrt((-self.real + Self::abs(self)) / T::num(2)),
        )
    }

    /// Returns the multiplicative inverse of this [`Complex`].
    pub fn inv(self) -> Self {
        Self::conj(self) / Self::square_abs(self)
    }

    /// Returns this [`Complex`] raised to a power using exponentiation by squaring.
    pub fn powi(self, exponent: i64) -> Self {
        match exponent {
            0 => Self::new(T::one(), T::zero()),
            1 => self,
            -1 => Self::inv(self),
            _ => {
                if exponent < 0 {
                    Self::inv(Self::powi(self, -exponent))
                } else if exponent.rem_euclid(2) == 0 {
                    Self::powi(self * self, exponent / 2)
                } else {
                    self * Self::powi(self * self, (exponent - 1) / 2)
                }
            }
        }
    }

    /// Returns this [`Complex`] raised to a power using De Moivre's formula.
    pub fn powf(self, exponent: T) -> Self {
        let arg_exponent = self.arg() * exponent;
        Self::new(T::cos(arg_exponent), T::sin(arg_exponent)) * T::powf(Self::abs(self), exponent)
    }

    /// Returns this [`Complex`] raised to a complex power.
    pub fn powc(self, exponent: Self) -> Self {
        Self::powf(self, exponent.real) * Self::exp(Self::ln(self) * Self::i() * exponent.imag)
    }

    /// Returns this [`Complex`] raised to a [`ComplexPolar`] power.
    pub fn powcp(self, exponent: ComplexPolar<T>) -> Self {
        Self::powc(self, exponent.unpolarize())
    }

    /// Returns e raised to the power of this [`Complex`].
    pub fn exp(self) -> Self {
        Self::new(T::cos(self.imag), T::sin(self.imag)) * T::exp(self.real)
    }

    /// Returns the natural logarithm of the absolute value of this [`Complex`].
    pub fn ln_abs(self) -> T {
        T::ln(Self::square_abs(self)) / T::num(2)
    }

    /// Returns the natural logarithm of this [`Complex`].
    pub fn ln(self) -> Self {
        Self::new(Self::ln_abs(self), Self::arg(self))
    }

    /// Returns the logarithm base 10 of this [`Complex`].
    pub fn log(self) -> Self {
        Self::ln(self) / T::ln(T::num(10))
    }

    /// Returns the logarithm base n of this [`Complex`].
    pub fn logn(self, base: T) -> Self {
        Self::ln(self) / T::ln(base)
    }
}

// Trig
impl<T: Float> Complex<T> {
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
        Self::inv(Self::tan(self))
    }

    /// Returns the secant of this [`Complex`].
    pub fn sec(self) -> Self {
        Self::inv(Self::cos(self))
    }

    /// Returns the cosecant of this [`Complex`].
    pub fn csc(self) -> Self {
        Self::inv(Self::sin(self))
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
        Self::arcsin(self / Self::sqrt(self.powi(2) + T::one()))
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
        Self::inv(Self::tanh(self))
    }

    /// Returns the hyperbolic secant of this [`Complex`].
    pub fn sech(self) -> Self {
        Self::inv(Self::cosh(self))
    }

    /// Returns the hyperbolic cosecant of this [`Complex`].
    pub fn csch(self) -> Self {
        Self::inv(Self::sinh(self))
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
        Self::ln((self + T::one()) / (-self + T::one())) / T::num(2)
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
impl<T: Float + Display> Display for Complex<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.imag >= T::zero() {
            write!(f, "{} + {}i", self.real, self.imag)
        } else {
            write!(f, "{} - {}i", self.real, self.imag)
        }
    }
}

/// Struct representing a complex number in polar form
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct ComplexPolar<T: Float> {
    radius: T,
    angle: T,
}

impl<T: Float> ComplexPolar<T> {
    /// Creates a new [`ComplexPolar`].
    pub fn new(radius: T, angle: T) -> Self {
        Self {
            radius: radius.abs(),
            angle: (angle
                + if radius < T::zero() {
                    T::pi()
                } else {
                    T::zero()
                })
                % (T::num(2) * T::pi()),
        }
    }

    /// Returns the imaginary number i
    pub fn i() -> Self {
        Self::new(T::one(), T::pi())
    }

    /// Returns the real part of this [`ComplexPolar`]
    pub fn real(self) -> T {
        self.radius * T::cos(self.angle)
    }

    /// Returns the imaginary part of this [`ComplexPolar`]
    pub fn imag(self) -> T {
        self.radius * T::sin(self.angle)
    }

    /// Returns the conjugate of this [`ComplexPolar`].
    pub fn conj(self) -> Self {
        Self::new(self.radius, self.angle)
    }

    /// Returns the absolute value of this [`ComplexPolar`].
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmath::ComplexPolar;
    /// use std::f32;
    ///
    /// let z1 = ComplexPolar::new(4.0, f32::consts::PI);
    /// let z2 = ComplexPolar::new(5.0, f32::consts::PI / 2.0);
    ///
    /// assert_eq!(ComplexPolar::abs(z1), 4.0);
    /// assert_eq!(ComplexPolar::abs(z2), 5.0);
    /// ```
    pub fn square_abs(self) -> T {
        self.radius.powi(2)
    }

    /// Returns the absolute value of this [`ComplexPolar`].
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmath::ComplexPolar;
    /// use std::f32;
    ///
    /// let z1 = ComplexPolar::new(4.0, f32::consts::PI);
    /// let z2 = ComplexPolar::new(5.0, f32::consts::PI / 2.0);
    ///
    /// assert_eq!(ComplexPolar::abs(z1), 4.0);
    /// assert_eq!(ComplexPolar::abs(z2), 5.0);
    /// ```
    pub fn abs(self) -> T {
        self.radius
    }

    /// Returns the argument of this [`ComplexPolar`].
    pub fn arg(self) -> T {
        self.angle
    }

    /// Returns this [`ComplexPolar`] in standard form as a [`Complex`]
    pub fn unpolarize(self) -> Complex<T> {
        Complex::exp(Complex::i() * self.angle) * self.radius
    }

    /// Returns the square root of this [`ComplexPolar`].
    pub fn sqrt(self) -> Self {
        Self::new(T::sqrt(self.radius), self.angle / T::num(2))
    }

    /// Returns the multiplicative inverse of this [`ComplexPolar`].
    pub fn inv(self) -> Self {
        Self::new(T::one() / self.radius, -self.angle)
    }

    /// Returns this [`ComplexPolar`] raised to float.
    pub fn powf(self, exponent: T) -> Self {
        Self::new(self.radius.powf(exponent), self.angle * exponent)
    }

    /// Returns this [`ComplexPolar`] raised to a [`Complex`] power.
    pub fn powc(self, exponent: Complex<T>) -> Self {
        Complex::powc(self.unpolarize(), exponent).polarize()
    }

    /// Returns this [`ComplexPolar`] raised to a complex polar power.
    pub fn powcp(self, exponent: ComplexPolar<T>) -> Self {
        Complex::powc(self.unpolarize(), exponent.unpolarize()).polarize()
    }

    /// Returns e raised to the power of this [`ComplexPolar`].
    pub fn exp(self) -> Self {
        self.unpolarize().exp().polarize()
    }

    /// Returns the natural logarithm of the absolute value of this [`ComplexPolar`].
    pub fn ln_abs(self) -> T {
        T::ln(self.abs())
    }

    /// Returns the natural logarithm of this [`ComplexPolar`].
    pub fn ln(self) -> Self {
        self.unpolarize().ln().polarize()
    }

    /// Returns the logarithm base 10 of this [`ComplexPolar`].
    pub fn log(self) -> Self {
        self.unpolarize().log().polarize()
    }

    /// Returns the logarithm base n of this [`Complex`].
    pub fn logn(self, base: T) -> Self {
        self.unpolarize().logn(base).polarize()
    }
}

mod overloading;
mod polaroverloading;

#[cfg(test)]
mod tests;
