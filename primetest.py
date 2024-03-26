import concurrent.futures
from gmpy2 import mpz, f_divmod_2exp, f_mod_2exp, bit_length, gcd
import gmpy2
import itertools
import random
import time
import math

def sq_adj(x, adj=0):
    return x * x + adj

def sq_mod_n(x, n, adj=0):
    return sq_adj(x, adj) % n

def mul_mod_n(x, a, n):
    return (x * a) % n


def powmod(a, x, n):
    """Compute a^x (mod n)"""
    running_total = a
    binary_representation = []
    i = x
    while i > 0:
        binary_representation.insert(0, i % 2)
        i = i // 2
    prev_total = 0
    for b in binary_representation[1:]:
        running_total = sq_mod_n(running_total, n)
        prev_total = running_total
        if b == 1:
            running_total = mul_mod_n(running_total, a, n)

    return prev_total


def isprime(n):
    # uses fermat's little theorem
    if n <= 3:
        return True
    elif n % 2 == 0:
        return False
    aa = range(2, min(n, 10))
    prev = [powmod(a, n, n) for a in aa]
    return all([p in (-1, 1) for p in prev])

def reduce_modulus(s, p):
    s_bottom, s_top = f_divmod_2exp(s, p)
    return s_bottom + s_top

def lucaslehmer(p):
    # uses lucas lehmer primality test for candidate primes 2^p - 1
    start = time.time()

    s = 4
    for i in range(p - 2):
        s = gmpy2.square(s)
        s = gmpy2.sub(s, 2)
        while bit_length(s) > p:
            s = reduce_modulus(s, p)
        if gmpy2.bit_scan0(s) == p:
            s = 0
    return s == 0, time.time() -  start, False

def prime_candidates(start, stop):
    n = 1 << (start - 1)
    for i in range(start, stop):
        n = n * 2
        if isprime(i):
            yield i


def mersennes(start, stop, use_threads=True):
    cpu_elapsed = 0
    start_time = time.time()
    results = []
    skips = 0
    count = 0
    last_output = time.time()
    if use_threads:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for i in prime_candidates(start, stop):
                future = executor.submit(lucaslehmer, i)
                results.append((i, future))
            print(f"Submitted (start={start}, stop={stop})")
            for i, future in results:
                q, elapsed, skipped = future.result()
                cpu_elapsed += elapsed
                if q:
                    print(f"2^{i}-1", q)
                    last_output = time.time()
                elif time.time() > last_output + 10:
                    print(f" ({i})")
                    last_output = time.time()
                count += 1
                if skipped:
                    skips += 1
            print("Wall clock", time.time() - start_time, "cpu time", cpu_elapsed)
    else:
        for i in prime_candidates(start, stop):
            q, elapsed, skipped = lucaslehmer(i)
            if q:
                print(f"2^{i}-1", q)
            elif time.time() > last_output + 10:
                print(f" ({i})")
                last_output = time.time()
            count += 1
            if skipped:
                skips += 1
    print(f"Skipped: {skips}/{count}")

if __name__ == '__main__':
    import cProfile
    start = time.time()
    cProfile.run("mersennes(1000, 5000, use_threads=False)", sort="cumtime")
#    mersennes(50000, 51000)
    print("elapsed", time.time() - start)
