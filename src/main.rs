use gcd::Gcd;
use itertools::Itertools;
use log::{debug, error, info};

fn isqrt(num: u64) -> u64 {
    let num2: f64 = num as f64;
    (num2.sqrt().round().floor() as i64) as u64
}

fn equal_rational(num1: u64, denom1: u64, num2: u64, denom2: u64) -> bool {
    num1 * denom2 == denom1 * num2
}

fn normalize(num: u64, denom: u64) -> (u64, u64) {
    (num / num.gcd(denom), denom / num.gcd(denom))
}

fn is_prime(num: u64) -> bool {
    if num == 2 { return true; }
    if num % 2 == 0 {
        return false;
    }
    for i in 1..isqrt(num) / 2 {
        if num % (2 * i + 1) == 0 {
            return false;
        }
    }
    true
}

fn is_prime_with_list(num: u64, primes: &[u64]) -> bool {
    for i in primes {
        if i * i > num {
            return true;
        }
        if num % i == 0 {
            return false;
        }
    }
    true
}
fn factorize_with_list(num2: u64, primes: &[u64]) -> Vec<u64> {
    let mut vec = Vec::new();
    let mut num = num2;
    while num != 1 {
        for p in primes {
            if num % p == 0 {
               num /= p;
               vec.push(*p);
               break;
            }
        }
    }
    return vec;
}

fn euler_totient(num: u64, primes: &[u64]) -> u64 {
  
  let mut factors = factorize_with_list(num, primes);
  factors.dedup();
  let mut count = num;
  for f in factors {
      count *= f - 1;
      count /= f;
  }
  return count;
}

fn primes_up_to(num: u64) -> Vec<u64> {
    (2..num).filter(|num| is_prime(*num)).collect()
}

fn check_quadratic(a: i64, b: i64) -> usize {
    (0..)
        .take_while(|n| (n * n + a * n + b).try_into().is_ok_and(is_prime))
        .count()
}

fn problem45() -> u64 {
    let triangle = |n| (n as u64) * (n + 1) / 2;
    let penta = |n| n * (3 * n - 1) / 2;
    let hexa = |n| n * (2 * n - 1);
    for i in 286u64.. {
        let t: u64 = triangle(i);
        let approx: u64 = ((i as f32) / 3_f32.sqrt()).floor() as u64;
        for p in approx - 2..approx + 2 {
            let p2 = penta(p);
            for h in i / 2 - 2..i / 2 + 2 {
                let h2 = hexa(h);
                if t == p2 && t == h2 {
                    log::debug!("i: {}", i);
                    return t;
                }
            }
        }
    }
    return 0;
}

fn problem50() -> u64 {
    let primes = primes_up_to(10001);
    log::debug!("primes: {}", primes.len());
    // treat sequences starting at 2 differently due to parity
    let i = primes.iter().scan(0 as u64, |sum, x| Some(*sum + x));
    let max = i
        .filter(|s| is_prime_with_list(*s, &primes))
        .max()
        .expect("wah");
    let (a, b) = (0..primes.len() - 1)
        .cartesian_product(1..primes.len())
        .max_by_key(|(i, j)| {
            if i >= j {
                return 0;
            };
            let s = &primes[*i as usize..=*j as usize].iter().sum();
            if is_prime_with_list(*s, &primes) {
                *s
            } else {
                0
            }
        })
        .expect("nonempty");
    max
}

fn problem72() -> u64 {
  let primes = primes_up_to(50003);
  (2..50000).map(|n| euler_totient(n, &primes)).sum()
}

fn problem33() -> u64 {
    let mut res_i: u64 = 1;
    let mut res_j: u64 = 1;
    for j in 1..100 {
        for i in 10..j {
            if i % 10 == 0 && j % 10 == 0 {
                continue;
            }
            if i % 10 == j / 10 && equal_rational(i, j, i / 10, j % 10) {
                let (norm_i, norm_j) = normalize(i, j);

                log::debug!("{} / {} = {} / {}", i, j, norm_i, norm_j);
                res_i *= norm_i;
                res_j *= norm_j;
            }
            if i / 10 == j % 10 && equal_rational(i, j, i % 10, j / 10) {
                log::debug!("{}, {}", i, j);
            }
        }
    }
    let (_, res) = normalize(res_i, res_j);
    res
}

fn problem27() -> i64 {
    log::debug!("n^2 + n + 41: {}", check_quadratic(1, 41));
    log::debug!("n^2 - 79n + 1601: {}", check_quadratic(-79, 1601));
    log::debug!("n^2 - 61n + 971: {}", check_quadratic(-61, 971));
    let (a, b) = ((-1000 as i64)..1000)
        .cartesian_product((-1000 as i64)..1000)
        .max_by_key(|(a, b)| check_quadratic(*a, *b))
        .expect("nonempty");
    a * b
}

fn main() {
    env_logger::init();
    use std::time::Instant;
    println!("27: {}", problem27());
    println!("33: {}", problem33());
    println!("45: {}", problem45());
    let now = Instant::now();
    println!("72: {}", problem72());
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    //  println!("50: {}",problem50());
}
