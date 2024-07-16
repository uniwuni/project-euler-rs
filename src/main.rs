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

fn is_semiprime_with_list_hashset(num: u64, primes: &[u64], primes2: &std::collections::HashSet<u64>) -> bool {
    if num % 1000000 == 0 { log::debug!("is_semiprime {num}"); }
    if num % 4 == 0 { return num == 4; }
    if num % 9 == 0 { return num == 9; }
    if num % 25 == 0 { return num == 25; }
    for i in primes {
        if &num % i == 0 {
            return is_prime_with_hashset(num / i, primes2);
        }
    }
    true
}

fn problem187() -> usize {
  let primes = primes_up_to(11000);
  let primes2: std::collections::HashSet<u64> = primes.clone().into_iter().collect();
  (2..100000000).filter(|n| is_semiprime_with_list_hashset(*n, &primes, &primes2)).count()
}

fn is_prime_with_list(num: u64, primes: &[u64]) -> bool {
    for i in primes {
        if num % i == 0 {
            return &num == i;
        }
    }
    true
}

fn is_prime_with_hashset(num: u64, primes: &std::collections::HashSet<u64>) -> bool {
    for i in [2,3,5,7,11,13,17,19,23,29,31,37,41,43] {
        if num == i {
            return true;
        }
        if num % i == 0 {
            return false;
        }
    }
    primes.contains(&num)
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
  if num % 10000 == 0 
    { log::debug!("totient: {}", num); }
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
  let primes = primes_up_to(1000003);
  (2..1000001).map(|n| euler_totient(n, &primes)).sum()
}

fn problem69() -> u64 {
  let primes = primes_up_to(1000003);
  let Some((n,k)) = (2..1000001).map(|n| (n, euler_totient(n, &primes)))
                                .max_by(|(n,k), (m,p)| (n * p).cmp(&(m * k)))
    else { panic!("how") };
  n
}

fn problem216() -> u64 {
  let primes = primes_up_to((17 * 10000)/10);
  0
}

fn pow_mod(base: u64, exponent: u64, m: u64) -> u64 {
    let e = exponent % m;
    let b = base % m;
    if b == 1 || e == 0 { return 1; }
    if b == 0 {return 0;}
    if e % 2 == 0 {
        pow_mod((base * base) % m, exponent/2, m)
    } else {
        (pow_mod((base * base) % m, exponent/2, m) * base) % m
    }
}

fn hyper_mod(base: u64, exponent: u64, m: u64) -> u64 {
  if(exponent == 1) { return base % m; }
  pow_mod(base, hyper_mod(base, exponent - 1, m),m)
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
    //println!("69: {}", problem69());
    //println!("13 ** 67 mod 123: {}", pow_mod(13,67,123))
  //  let now = Instant::now();
  //  println!("72: {}", problem72());
  //  let elapsed = now.elapsed();
    println!("187: {}", problem187());
    println!("188: {}", hyper_mod(1777, 1855, 100000000));
  //  println!("Elapsed: {:.2?}", elapsed);
    //  println!("50: {}",problem50());
}
