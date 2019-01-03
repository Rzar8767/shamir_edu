extern crate rand;

use rand::{
    distributions::{Distribution, Uniform},
    thread_rng,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn general_test() {
        general_usage(125, 3);
    }

    #[test]
    #[should_panic]
    fn general_failure() {
        general_usage(6326214, 3);
    }

    fn general_usage(secret: i64, threshold: u8) {
        let shamir_test = SecretData::from_secret(secret, threshold, None);
        let mut points: Vec<i64> = vec![];
        let mut values: Vec<i64> = vec![];
        for x in 1..=3 {
            points.push(x);
            values.push(shamir_test.get_share(x as u8).1);
        }
        assert_eq!(secret, SecretData::recover_secret(threshold, &points, &values, shamir_test.prime).unwrap());
    }
}

pub struct SecretData {
    pub secret_data: Option<i64>,
    pub coefficients: Vec<i64>,
    pub prime: i64,
}

impl SecretData {
    pub fn from_secret(secret: i64, threshold: u8, prime: Option<i64>) -> SecretData {
        let prime = match prime {
            Some(num) => num,
            None => 6326213,
        };
        if secret >= prime {
            panic!("Secret will be lost during modulo operation")
        }

        let mut coefficients: Vec<i64> = vec![secret];
        let between = Uniform::from(1..prime);
        let mut rng = thread_rng();
        coefficients.extend(between.sample_iter(&mut rng).take((threshold - 1) as usize));

        SecretData {
            secret_data: None,
            coefficients: coefficients,
            // all coefficients have to be lower than it and it has to be known when trying to retrieve the secret
            prime: prime,
        }
    }

    pub fn get_share(&self, id: u8) -> (i64, i64) {
        if id < 1 {
            // Don't need to check for id > 255 because u8 can't be
            panic!("share id must be between 1 and 255");
        }
        (
            id as i64,
            poly_curve_eval(&self.coefficients, id, self.prime),
        )
    }

    pub fn recover_secret(
        threshold: u8,
        points: &[i64],
        values: &[i64],
        prime: i64,
    ) -> Option<i64> {
        if threshold as usize > values.len() {
            println!("Number of shares is below the threshold");
            return None;
        }
        if points.len() != values.len() {
            return None;
        }

        Some(lagrange_interpolation_at_zero(points, values, prime))
    }
}

pub fn poly_curve_eval(coefficients: &Vec<i64>, x: u8, modulo: i64) -> i64 {
    let mut accu: i64 = 0;
    for coefficient in coefficients.iter().rev() {
        accu *= x as i64;
        accu += coefficient;
        accu %= modulo;
    }
    accu
}

// mathematical functions below borrowed from https://github.com/snipsco/rust-threshold-secret-sharing/blob/master/src/numtheory.rs
// used only for educational purposes
pub fn lagrange_interpolation_at_zero(points: &[i64], values: &[i64], prime: i64) -> i64 {
    assert_eq!(points.len(), values.len());
    // Lagrange interpolation for point 0
    let mut acc = 0i64;
    for i in 0..values.len() {
        let xi = points[i];
        let yi = values[i];
        let mut num = 1i64;
        let mut denum = 1i64;
        for j in 0..values.len() {
            if j != i {
                let xj = points[j];
                num = (num * xj) % prime;
                denum = (denum * (xj - xi)) % prime;
            }
        }
        acc = (acc + yi * num * mod_inverse(denum, prime)) % prime;
    }
    acc
}

/// Inverse of `k` in the *Zp* field defined by `prime`.
pub fn mod_inverse(k: i64, prime: i64) -> i64 {
    let k2 = k % prime;
    let r = if k2 < 0 {
        -gcd(prime, -k2).2
    } else {
        gcd(prime, k2).2
    };
    (prime + r) % prime
}

/// Euclidean GCD implementation (recursive). The first member of the returned
/// triplet is the GCD of `a` and `b`.
pub fn gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let n = a / b;
        let c = a % b;
        let r = gcd(b, c);
        (r.0, r.2, r.1 - r.2 * n)
    }
}
