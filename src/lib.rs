#[derive(Debug, PartialEq)]
struct Complex {
    real: f64,
    imag: f64,
}

impl Complex {
    fn new(real: f64, imag: f64) -> Complex {
        Complex {
            real: real,
            imag: imag,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = Complex::new(2f64, 2f64);
        assert_eq!(
            result,
            Complex {
                real: 2f64,
                imag: 2f64
            }
        );
    }
}
