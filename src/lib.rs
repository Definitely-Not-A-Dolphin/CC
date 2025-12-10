use std::f64::consts::TAU;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Complex {
    real: f64,
    imag: f64,
}

fn factorial(n: u32) -> u64 {
    if n == 0 {
        return 1;
    };
    let mut prod = 1u64;
    for i in 1..=(n as u64) {
        prod *= i;
    }
    prod
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
        if self.imag > 0f64 {
            return f64::acos(self.real / self.abs());
        }
        TAU - f64::acos(self.real / self.abs())
    }

    /// Returns the half-argument square root of this [`Complex`].
    pub fn sqrt(self) -> Complex {
        let root_real = ((self.real + self.abs()) / 2f64).sqrt();
        let root_imag = ((-self.real + self.abs()) / 2f64).sqrt();
        if self.imag >= 0f64 {
            return Complex::new(root_real, root_imag);
        };
        Complex::new(-root_real, root_imag)
    }

    /// Returns the multiplicative inverse of this [`Complex`].
    pub fn inv(self) -> Complex {
        Complex::new(
            self.real / self.square_abs(),
            -self.imag / self.square_abs(),
        )
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
                    result = multiply(&result, &self);
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
        Complex::new(
            self.abs().powf(exponent) * f64::cos(exponent * self.arg()),
            self.abs().powf(exponent) * f64::sin(exponent * self.arg()),
        )
    }

    pub fn exp(self) -> Complex {
        let mut result = Complex::new(0f64, 0f64);
        for i in 0..=216 {
            result = Complex::new(
                self.powi(i).real / factorial(i as u32) as f64,
                self.powi(i).imag / factorial(i as u32) as f64,
            );
        }
        return result;
    }
}

/// Returns the sum to two [`Complex`] numbers.
pub fn add(a: &Complex, b: &Complex) -> Complex {
    Complex::new(a.real() + b.real(), a.imag() + b.imag())
}

/// Returns the difference of two [`Complex`] numbers.
pub fn subtract(a: &Complex, b: &Complex) -> Complex {
    Complex::new(a.real() - b.real(), a.imag() - b.imag())
}

/// Returns the product of two [`Complex`] numbers.
pub fn multiply(z1: &Complex, z2: &Complex) -> Complex {
    Complex::new(
        z1.real() * z2.real() - z1.imag() * z2.imag(),
        z1.real() * z2.imag() + z1.imag() * z2.real(),
    )
}

/// Returns the quotient of two [`Complex`] numbers.
pub fn divide(z1: &Complex, z2: &Complex) -> Complex {
    let common_denominator = z2.real().powi(2) + z2.imag().powi(2);
    Complex::new(
        (z1.real() * z2.real() + z1.imag() * z2.imag()) / common_denominator,
        (z1.imag() * z2.real() - z1.real() * z2.imag()) / common_denominator,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_binary_operators() {
        let z1 = Complex::new(3f64, 4f64);
        let z2 = Complex::new(-2.5, 6.23);
        assert_eq!(add(&z1, &z2), Complex::new(0.5, 10.23));
        // God I love floating point arithmetic
        assert_eq!(subtract(&z1, &z2), Complex::new(5.5, -2.2300000000000004));
        assert_eq!(multiply(&z1, &z2), Complex::new(-32.42, 8.690000000000001));
        assert_eq!(
            divide(&z1, &z2),
            Complex::new(0.38657077107776017, -0.6366656384742215)
        );
    }

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
        assert_eq!(z2.arg(), 6.111806180790079);

        // inv
        assert_eq!(z1.inv(), Complex::new(0.12, -0.16));
        assert_eq!(
            z2.inv(),
            Complex::new(0.1867145421903052, 0.03231597845601436)
        );
    }

    #[test]
    fn complex_binary_operators() {
        let z1 = Complex::new(3f64, 4f64);
        let z2 = Complex::new(5.2, -0.9);

        // powf
        assert_eq!(z1.powf(3f64), Complex::new(-117f64, 43.99999999999998));
        assert_eq!(
            z2.powf(-2.5),
            Complex::new(-0.014217542838549323, -0.006493773098977875)
        );
    }
}
