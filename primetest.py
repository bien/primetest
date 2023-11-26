import numpy as np
import random

def sq_adj(x, adj=0):
    return x * x + adj

def just_mod(x, n):
    return x % n

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

def mersennes(k, start=2):
    n = 2 ** (start - 1)
    for i in range(start, k):
        n = n * 2
        if isprime(i):
            q = lucaslehmer(i)
#            q = isprime(n - 1)
            if q:
                print(f"2^{i}-1", q)

#print(isprime(17))
import time
import cProfile
start = time.time()
cProfile.run("mersennes(10000)")
print("elapsed", time.time() - start)
