import concurrent.futures
from gmpy2 import mpz, f_divmod_2exp, f_mod_2exp, bit_length, gcd
import gmpy2
import itertools
import random
import time
import math

def sq_adj(x, adj=0):
    return x * x + adj

def just_mod(x, n):
    return x % n

def just_div(x, n):
    return x // n

def just_mult(a, b):
    return a * b

def just_mult1(a, b):
    return a * b

def just_mult2(a, b):
    return a * b

def just_mult3(a, b):
    return a * b

def just_mod_bits(x, nbits):
    return f_mod_2exp(x, nbits)

def div_bits(x, nbits):
    return x >> nbits

def sq_mod_n(x, n, adj=0):
    return just_mod(sq_adj(x, adj), n)

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

def powmod2(a, x, n):
    return powmod(a, x, n)

KNOWN_M = {}
def pollard_factor(n, a=3, b=100):
#    else:
#        m = 1
#        for q in range(3, b + 1):
#            if isprime(q):
#                m *= pow(q, math.floor(math.log(b, q)))
#                m *= q
#        k = a ** m - 1
    if b in KNOWN_M:
        m = KNOWN_M[b]
    else:
        m = 1
        factors = []
        for q in itertools.chain([2], range(3, b + 1, 2)):
            if isprime(q):
                factors.append(pow(q, int(math.log(b, q))))
        for factor in factors:
            m *= factor
        KNOWN_M[b] = m
    k = powmod2(a, m, n) - 1
    g = gcd(k, n)
    if g == n and b > 5:
        return pollard_factor(n, b=b//2)
    else:
        return g

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

def just_eq(a, b):
    return a == b

def just_sub(a, b):
    return a - b

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

def lucaslehmer_montgomery(p):
    # uses lucas lehmer primality test for candidate primes 2^p - 1
    start = time.time()

    n = (1 << p) - 1
#    factor = pollard_factor(mpz(n), b=2 + p)
#    if n > factor > 1:
#        print(f"factor 2^{p}-1: {factor}")
 #       return False, time.time() - start, True
    r, rbits = next_power_2(n)
    s = Montgomery(4, r, n, rbits)
    for i in range(p - 2):
        s.mult(s)
        s.add(-2)
    return s.to_int() == 0, time.time() -  start, False

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
        m = just_mod_bits(just_mult2(just_mod_bits(t, self.rbits), self.nprime), self.rbits)
        t = div_bits((t + just_mult3(m, self.n)), self.rbits)
        if t >= self.n:
            return t - self.n
        else:
            return t
        
#    def multi_redc(self, t, word_length=8):
#        """T is little-endian"""
#        rwords = math.ceil(self.rbits / word_length)
 #       nwords = math.ceil(self.nbits / word_length)
  #      word_bits = (1 << word_length) - 1
   #     t_array = [0] * (rwords + nwords)
    #    for i in range(nwords):
     #       t_array[i] = nwords >> (i * word_length) & ((1 << word_length) - 1)
      #  for i in range(rwords):
#            c = 0
#            m = just_mod_bits(t_array[i] * self.nprime, self.word_bits)
#            for j in range(nwords):
#                x = t_array[i + j] + m * 


    def mult(self, a):
        self.m = self.redc(just_mult1(self.m, a.m))

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
        return lucaslehmer_montgomery(self.i)

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
    cProfile.run("mersennes(1000, 5000, use_threads=False)", sort="cumtime")
#    mersennes(50000, 51000)
    print("elapsed", time.time() - start)

