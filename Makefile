LDFLAGS = /Users/michael/primetest/libgmp/gmp-6.3.0/.libs/libgmp.a
#LDFLAGS = -lgmp
CFLAGS = -O2 -Wall -g -I/Users/michael/primetest/libgmp/gmp-6.3.0

primetest: primetest.o

clean:
	rm -f primetest *.o

# $ time ./primetest 10000 15000
# brew libgmp: real=1m11.434s user=1m9.895
# stock libgmp-6.3.0: real=1m10.636s user=1m9.582s
