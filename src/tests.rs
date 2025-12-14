#![allow(clippy::excessive_precision)]
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
    assert_eq!(
        1.8 / z2,
        Complex::new(-0.099860415, -0.24885215)
    );
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
