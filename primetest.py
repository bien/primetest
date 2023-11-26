import concurrent.futures
from gmpy2 import mpz
import random
import time

def sq_adj(x, adj=0):
    return x * x + adj

def just_mod(x, n):
    return x % n

def just_div(x, n):
    return x // n

def just_mult(a, b):
    return a * b

def just_mod_2(x, n):
    """n must be power of 2"""
    return x & (n - 1)

def div_bits(x, nbits):
    return x >> nbits

def sq_mod_n(x, n, adj=0):
    return just_mod(sq_adj(x, adj), n)

def mul_mod_n(x, a, n):
    return (x * a) % n

def powmod(a, n):
    """Compute a^n (mod n)"""
    running_total = a
    binary_representation = []
    i = n
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
    aa = random.choices(range(1, min(10000, n - 1)), k=5)
    prev = [powmod(a, n) for a in aa]
    return all([p in (-1, 1) for p in prev])

def lucaslehmer_montgomery(p):
    # uses lucas lehmer primality test for candidate primes 2^p - 1
    start = time.time()
    n = (1 << p) - 1
    r, rbits = next_power_2(n)
    s = Montgomery(4, r, n, rbits)
    for i in range(p - 2):
        s.mult(s)
        s.add(-2)
#        print(n, i, s)
    return s.to_int() == 0, time.time() - start

def extended_gcd(a, b):
    s = 0
    old_s = 1
    r = b
    old_r = a
    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
    if old_s < 0:
        old_s += b
    return old_s

def to_montgomery_form(x, r, n):
    return (x * r) % n

def from_montgomery_form(m, rprime, n):
    return (m * rprime) % n

class Montgomery:
    def __init__(self, x, r, n, rbits=0):
        self.n = n
        self.m = mpz(to_montgomery_form(x, r, n))
        self.r = r
        self.rbits = rbits
        self.rprime = extended_gcd(r, n)
        self.nprime = r + -extended_gcd(n, r)

    def to_int(self):
        return from_montgomery_form(self.m, self.rprime, self.n)
    
    def redc(self, t):
        m = just_mod_2(just_mult(just_mod_2(t, self.r), self.nprime), self.r)
        t = div_bits((t + just_mult(m, self.n)), self.rbits)
        if t >= self.n:
            return t - self.n
        else:
            return t

    def mult(self, a):
        self.m = self.redc(just_mult(self.m, a.m))

    def add(self, x):
        self.m += x
        if self.m >= self.n:
            self.m -= self.n
        if self.m < 0:
            self.m += self.n

def prime_candidates(start, stop):
    n = 1 << (start - 1)
    for i in range(start, stop):
        n = n * 2
        if isprime(i):
            yield i

class FakeFuture:
    def __init__(self, i):
        self.i = i

    def result(self):
        start = time.time()
        return lucaslehmer_montgomery(self.i)

def mersennes(start, stop, use_threads=True):
    cpu_elapsed = 0
    start_time = time.time()
    results = []
    if use_threads:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for i in prime_candidates(start, stop):
                future = executor.submit(lucaslehmer_montgomery, i)
                results.append((i, future))
            print("Submitted")
    else:
        for i in prime_candidates(start, stop):
            future = FakeFuture(i)
            results.append((i, future))

    for i, future in results:
        q, elapsed = future.result()
        cpu_elapsed += elapsed
        if q:
            print(f"2^{i}-1", q)
        elif q % 100 == 0:
            print(f" (i)")
    print("Wall clock", time.time() - start_time, "cpu time", cpu_elapsed)

def next_power_2(n):
    x = 1
    xbits = 0
    while x < n:
        x = x * 2
        xbits += 1
    return x, xbits

if __name__ == '__main__':
    import cProfile
    start = time.time()
    cProfile.run("mersennes(30000, 80000, use_threads=True)", sort="cumtime")
    print("elapsed", time.time() - start)

