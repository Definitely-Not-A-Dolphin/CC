#![allow(clippy::excessive_precision)]
use std::f64;

use super::*;

#[test]
fn unary_operators() {
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
fn binary_operators() {
    let z1 = Complex::new(3f32, 4f32);
    let z2 = Complex::new(-2.5, 6.23);

    // add
    assert_eq!(z1 + 5f32, Complex::new(8f32, 4f32));
    assert_eq!(5f32 + z1, Complex::new(8f32, 4f32));
    assert_eq!(z1 + z2, Complex::new(0.5, 10.23));

    // subtract
    assert_eq!(z1 - 5f32, Complex::new(-2f32, 4f32));
    assert_eq!(5f32 - z1, Complex::new(2f32, -4f32));
    assert_eq!(z1 - z2, Complex::new(5.5, 10.23));

    // multiply
    assert_eq!(z1 * -0.5, Complex::new(-1.5, -2f32));
    assert_eq!(1.8 * z2, Complex::new(-4.5, 11.214));
    assert_eq!(
        z1 * z2,
        Complex::new(-32.42, z1.real * z2.imag + z1.imag * z2.real)
    );

    // dividing
    assert_eq!(z1 / -0.5, Complex::new(-6f32, -8f32));
    assert_eq!(1.8 / z2, Complex::new(-0.099860415, -0.24885215));
    assert_eq!(
        z1 / z2,
        Complex::new(
            0.38657075,
            z1.real * z2.inv().imag + z1.imag * z2.inv().real
        )
    );

    let z3 = Complex::new(3f64, 4f64);
    let z4 = Complex::new(5.2, -0.9);

    // powi
    assert_eq!(z3.powi(6), Complex::new(11753f64, -10296f64));
    assert_eq!(
        z4.powi(-2),
        Complex::new(0.03381799780176568, 0.012067726245692974)
    );

    // powf
    assert_eq!(z3.powf(3f64), Complex::new(-117f64, 43.99999999999998));
    assert_eq!(
        z4.powf(-2.5),
        Complex::new(0.014217542838549325, 0.00649377309897787)
    );

    // powc
    assert_eq!(
        z3.powc(z4),
        Complex::new(-9667.467998399987, -2282.430226186542)
    );
    assert_eq!(
        z4.powc(z3),
        Complex::new(288.7067987011787, -41.7623644411436)
    );
}

#[test]
fn complex_trig() {
    let z1 = Complex::new(3f64, 4f64);

    // sin
    assert_eq!(
        Complex::sin(z1),
        Complex::new(3.853738037919377, -27.016813258003932)
    );

    // cos
    assert_eq!(
        Complex::cos(z1),
        Complex::new(-27.034945603074224, -3.851153334811777)
    );

    // tan
    assert_eq!(
        Complex::tan(z1),
        Complex::new(-0.00018734620462948492, 0.999355987381473)
    );

    // cot
    assert_eq!(
        Complex::cot(z1),
        Complex::new(-0.00018758773798366324, -1.0006443924715591)
    );

    // sec
    assert_eq!(
        Complex::sec(z1),
        Complex::new(-0.03625349691586887, 0.005164344607753179)
    );

    // csc
    assert_eq!(
        Complex::csc(z1),
        Complex::new(0.005174473184019398, 0.03627588962862602)
    );
}

#[test]
fn complex_inverse_trig() {
    let z1 = Complex::new(3f64, 4f64);

    // arcsin
    assert_eq!(
        Complex::arcsin(z1),
        Complex::new(0.6339838656391737, 2.305509031243471)
    );

    // arccos
    assert_eq!(
        Complex::arccos(z1),
        Complex::new(0.9368124611557229, -2.305509031243471)
    );

    // arctan
    assert_eq!(
        Complex::arctan(z1),
        Complex::new(1.4483069952314647, 0.15899719167999923)
    );

    // arccot
    assert_eq!(
        Complex::arccot(z1),
        Complex::new(0.12248933156343211, -0.15899719167999918)
    );

    // arcsec
    assert_eq!(
        Complex::arcsec(z1),
        Complex::new(1.4520455954874867, 0.1604455337745046)
    );

    // arccsc
    assert_eq!(
        Complex::arccsc(z1),
        Complex::new(0.11875073130740989, -0.1604455337745046)
    );
}

#[test]
fn complex_hyperbolic_trig() {
    let z1 = Complex::new(3f64, 4f64);

    // sinh
    assert_eq!(
        Complex::sinh(z1),
        Complex::new(-6.5481200409110025, -7.61923172032141)
    );

    // cosh
    assert_eq!(
        Complex::cosh(z1),
        Complex::new(-6.580663040551157, -7.581552742746545)
    );

    // tanh
    assert_eq!(
        Complex::tanh(z1),
        Complex::new(1.000709536067233, 0.004908258067495952)
    );

    // coth
    assert_eq!(
        Complex::coth(z1),
        Complex::new(0.9992669278059017, -0.004901182394304371)
    );

    // sech
    assert_eq!(
        Complex::sech(z1),
        Complex::new(-0.06529402785794704, 0.07522496030277322)
    );

    // csch
    assert_eq!(
        Complex::csch(z1),
        Complex::new(-0.0648774713706355, 0.0754898329158637)
    );
}

#[test]
fn complex_inverse_hyperbolic_trig() {
    let z1 = Complex::new(3f64, 4f64);

    // arcsinh
    assert_eq!(
        Complex::arcsinh(z1),
        Complex::new(2.2999140408792695, 0.9176168533514786)
    );

    // arccosh
    assert_eq!(
        Complex::arccosh(z1),
        Complex::new(2.305509031243477, 0.9368124611557199)
    );

    // arctanh
    assert_eq!(
        Complex::arctanh(z1),
        Complex::new(0.11750090731143398, 1.4099210495965755)
    );

    // arccoth
    assert_eq!(
        Complex::arccoth(z1),
        Complex::new(0.11750090731143398, -0.16087527719832115)
    );

    // arcsech
    assert_eq!(
        Complex::arcsech(z1),
        Complex::new(0.1604455337745046, -1.4520455954874867)
    );

    // arccsch
    assert_eq!(
        Complex::arccsch(z1),
        Complex::new(0.12124561370968728, -0.15950663187736328)
    );
}
