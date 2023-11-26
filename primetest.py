from itertools import zip_longest
import numpy as np
import random

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
#    return (x * x + adj) % n
#    return np.mod(np.power(x, 2) + adj, n)

def mul_mod_n(x, a, n):
    return np.mod(x * a, n)

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
    prev = powmod(np.array(aa), n)
    return np.all(np.isin(prev, (-1, 1)))

def lucaslehmer(p):
    # uses lucas lehmer primality test for candidate primes 2^p - 1
    s = 4
    n = (1 << p) - 1
    for i in range(p - 2):
        s = sq_mod_n(s, n, adj=-2)
#        print(n, i, s)
    return s == 0

def lucaslehmer_montgomery(p):
    # uses lucas lehmer primality test for candidate primes 2^p - 1
    n = (1 << p) - 1
    r, rbits = next_power_2(n)
    s = Montgomery(4, r, n, rbits)
    for i in range(p - 2):
        s.mult(s)
        s.add(-2)
#        print(n, i, s)
    return s.to_int() == 0

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
        self.m = to_montgomery_form(x, r, n)
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
        while self.m < 0:
            self.m += self.n

def prime_candidates(k, start=2):
    n = 2 ** (start - 1)
    for i in range(start, k):
        n = n * 2
        if isprime(i):
            yield i

def grouper(iterable, n, *, incomplete='fill', fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def mersennes(k, start=2):
    for i in prime_candidates(k, start):
        q = lucaslehmer_montgomery(i)
        if q:
            print(f"2^{i}-1", q)

def next_power_2(n):
    x = 1
    xbits = 0
    while x < n:
        x = x * 2
        xbits += 1
    return x, xbits

#print(isprime(17))
import time
import cProfile
start = time.time()
cProfile.run("mersennes(2000)", sort='cumtime')
print("elapsed", time.time() - start)

#print(to_montgomery_form(3, 100, 17))
#n = 17
#r = 32
#mm = Montgomery(7, r, n)
#mm.mult(Montgomery(15, r, n))
#print(mm.to_int())