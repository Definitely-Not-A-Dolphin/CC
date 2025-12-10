#[derive(Debug, PartialEq)]
pub struct Rational {
    numerator: i64,
    denominator: i64,
}

impl Rational {
    /// Creates a new [`Rational`].
    pub fn new(numerator: i64, denominator: i64) -> Rational {
        Rational {
            numerator: numerator,
            denominator: denominator,
        }
    }

    pub fn numerator(&self) -> i64 {
        self.numerator
    }

    pub fn denominator(&self) -> i64 {
        self.denominator
    }
}

pub fn reduce(q: Rational) {
    let denominator = q.denominator();
    let mut denominator_factors: Vec<i64> = vec![];
    for i in 1..=denominator {
        if denominator as f64 / i as f64 == (denominator as f64 / i as f64).floor() {
            denominator_factors.push(i as i64);
        }
    }
    println!("{:?}", denominator_factors);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thing() {
        assert_eq!(5.9 as i32, 5 as i32);
        assert_eq!(Rational::new(2, 3), Rational::new(2, 3))
    }
}
