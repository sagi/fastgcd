CFLAGS=-std=c99 -Wall -O4 -pthread -I./gmp-patched/include/ -Dmpz_raw_64
LDFLAGS=-static -L./gmp-patched/lib -lgmp

fastgcd: fastgcd.c
	$(CC) $(CFLAGS) $? $(LDFLAGS) -o $@

