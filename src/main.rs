extern crate shamir_naive;

fn main() {
    let secret = 125;
    let threshold = 3;
    let shamir_test = shamir_naive::SecretData::from_secret(secret, threshold, None);
    println!("prime: {}, secret: {:?}, coefficients: {:?}", shamir_test.prime, shamir_test.secret_data, shamir_test.coefficients);
    let mut points: Vec<i64> = vec![];
    let mut values: Vec<i64> = vec![];
    
    for x in 1..=3
    {
        points.push(x);
        values.push(shamir_test.get_share(x as u8).1);
    }

    println!("recovered secret is: {:?}", shamir_naive::SecretData::recover_secret(threshold, &points, &values, shamir_test.prime));
}
