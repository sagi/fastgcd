//
// fastgcd -- Efficient implementation of all-pairs GCD
// NH,AH 2012/03
//

// Usage: 
//   fastgcd INPUT
//   (where INPUT is a file containing hex-encoded RSA moduli)
//

// Notes:
// When operating on a large set of moduli, requires extreme amounts
// of RAM and disk space.  For our dataset, requires ~150 GB diskspace
// and 30 GB RAM.
// We require a patched version of GMP because the current version
// doesn't support storing integers larger than 2^31 bytes.

// TODO:
//   - better resume ability
//   - skip some intermediate tree levels for better parallelism?
//   - test kernel tuning for better memory/cache tradeoffs

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <pthread.h>
#include <gmp.h>

#define NTHREADS 4 // Get from compile-time argument?

#ifdef mpz_raw_64 // if patched gmp, use large format int i/o
#define __mpz_inp_raw mpz_inp_raw_64
#define __mpz_out_raw mpz_out_raw_64
#else // otherwise use normal i/o...beware 2^31 byte size limit
#define __mpz_inp_raw mpz_inp_raw
#define __mpz_out_raw mpz_out_raw
#endif

#define INPUT_FN    "input.mpz"
#define OUTPUT_FN   "output.mpz"

//#define min(x,y) (((x)<(y)) ? (x) : (y))

// return current time as a double, with usec precision
double now()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec / 1000000.;
}

// test if file exists
int file_exists(char *filename)
{
  FILE *in;
  if ((in = fopen(filename, "r"))) {
    fclose(in);
    return 1;
  }
  return 0;
}

// read hex strings from infile and write final count followed by gmp
// binary format values to output
void prep_hex_input(char *infile, char *outfile)
{
  double start = now();
  fprintf(stderr, "preprocessing input from %s\n", infile);

  FILE *in = fopen(infile, "r");
  assert(in);
  FILE *out = fopen(outfile, "wb");
  assert(out);

  int count=0;
  fwrite(&count, sizeof(count), 1, out);
  mpz_t x;
  mpz_init(x);
  for (;;) {
    int res = gmp_fscanf(in, "%Zx", x);
    if (res == EOF)
      break;
    if (res != 1) {
      fprintf(stderr, "invalid input\n");
      exit(1);
    }
    __mpz_out_raw(out, x);
    count++;
  }
  fclose(in);
  rewind(out);
  fwrite(&count, sizeof(count), 1, out);
  fclose(out);

  fprintf(stderr, "preprocessing %d elements took %0.3fs\n", count, now()-start);
}

typedef struct vec_ {
  mpz_t *el;
  int count; 
} vec_t;

// init vector v to contain count mpzs
void init_vec(vec_t *v, int count)
{
  assert(v);
  v->count = count;
  v->el = malloc(count * sizeof(mpz_t));
  assert(v->el);
  for (int i=0; i < v->count; i++)
    mpz_init(v->el[i]);
}

// free the vector v
void free_vec(vec_t *v)
{
  assert(v);
  for (int i=0; i < v->count; i++)
    mpz_clear(v->el[i]);
  free(v->el);
}

// initializes vec_t *v and fills it with contents of named binary format file
void input_bin_array(vec_t *v, char *filename)
{
  double start = now();
  fprintf(stderr, "reading %s...", filename);
  FILE *in = fopen(filename, "rb");
  assert(in);
  int count;
  int ret = fread(&count, sizeof(count), 1, in);
  assert(ret == 1);
  assert(count >= 0);
  init_vec(v, count);
  size_t bytes = 0;
  for (int i=0; i < count; i++)
    bytes += __mpz_inp_raw(v->el[i], in);
  fclose(in);
  fprintf(stderr, "%d elements, %zu bytes (%0.3fs)\n", v->count, bytes, now()-start);
}

// writes vec_t *v to the named file in binary format
void output_bin_array(vec_t *v, char *filename)
{
  double start = now();
  fprintf(stderr, "writing %s...", filename);
  FILE *out = fopen(filename, "wb");
  assert(out);
  fwrite(&v->count, sizeof(v->count), 1, out);
  size_t bytes = 0;
  for (int i=0; i < v->count; i++)
    bytes += __mpz_out_raw(out, v->el[i]);
  fclose(out);
  fprintf(stderr, "%d elements, %zu bytes (%0.3fs)\n",v->count, bytes, now()-start);
}

// writes vec_t *v to the named file as lines of hex strings
void output_hex_array(vec_t *v, char *filename)
{
  fprintf(stderr, "writing %s...", filename);
  FILE *out = fopen(filename,"w");
  assert(out);
  for (int i=0; i<v->count; i++)
    gmp_fprintf(out,"%Zx\n",v->el[i]);
  fclose(out);
  fprintf(stderr, "ok\n");
}

// sort and uniqify v
void uniq(vec_t *v)
{
  qsort(v->el, v->count, sizeof(mpz_t), (__compar_fn_t)mpz_cmp);

  int size = 0;
  for (int i=0; i < v->count; i++) {
    if (mpz_cmp(v->el[size],v->el[i]))
      mpz_set(v->el[++size],v->el[i]);
  }
  for (int i=size+1; i < v->count; i++)
    mpz_clear(v->el[i]);
  v->count = size+1;
}

// Executes func(n) over the range [start,end) using NTHREADS
// worker threads.  This is essentially a parallel version of:
//    for (int n=start; n < end; n++) { func(n); } 
// You are responsible for ensuring that func() is thread-safe!
void iter_threads(int start, int end, void (*func)(int n))
{
  int n = start;
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

  void *thread_body(void *ptr)
  {
    double start = now();
    for (;;) {
      pthread_mutex_lock( &mutex );
      int i = (n)++;
      pthread_mutex_unlock( &mutex );
      if (i >= end)
        break;
      func(i);
    }
    fprintf(stderr, "(ok %0.3fs)", now()-start);
    return NULL;
  }

  pthread_t thread_id[NTHREADS];
  for (int i=0; i<NTHREADS; i++)
    pthread_create(&thread_id[i], NULL, thread_body, NULL);
  for (int i=0; i<NTHREADS; i++)
    pthread_join(thread_id[i], NULL);
  fprintf(stderr, "\n");
}


int product_tree() {
  vec_t v;

  double tstart = now();
  fprintf(stderr, "multiplying numbers...\n");
    
  input_bin_array(&v, INPUT_FN);
  
  int level=0;
  while (v.count > 1) {
    double start = now();
    fprintf(stderr, "level %d\n", level);
    vec_t w;
    init_vec(&w,(v.count+1)/2);

    void mul_job(int i) {
      mpz_mul(w.el[i], v.el[2*i], v.el[2*i+1]);
    }
    iter_threads(0, v.count/2, mul_job);
    if (v.count & 1)
      mpz_set(w.el[v.count/2],v.el[v.count-1]); 
    
    char name[255];
    snprintf(name, sizeof(name)-1, "p%d.mpz", level);
    output_bin_array(&w, name);
      
    free_vec(&v);
    v = w;
    level++;
    fprintf(stderr, "%0.3fs\n", now()-start);
  }

  free_vec(&v);
  fprintf(stderr, "product tree took %0.3fs\n", now()-tstart);
  return level;
}

void remainder_tree(int level)
{
  double tstart = now();
  fprintf(stderr, "computing remainder tree\n");

  char name[255];
  snprintf(name, sizeof(name)-1, "p%d.mpz", level);
  vec_t P, v;
  input_bin_array(&P, name);
  
  /* Potential speedup:
  init_vec(&v,2);
  mpz_init(v.el[0],P.el[0]);
  mpz_init(v.el[1],P.el[0]);
  level--;
  P = v;
  */

  while (level > 0) {
    double start = now();
    fprintf(stderr, "level %d\n", level);
    level--;
    snprintf(name, sizeof(name)-1, "p%d.mpz", level);
    input_bin_array(&v, name);

    void mul_job(int i) {
      mpz_t s;
      mpz_init(s);
      mpz_mul(s, v.el[i], v.el[i]);
      mpz_mod(v.el[i], P.el[i/2], s);
      mpz_clear(s);
    }
    iter_threads(0, v.count, mul_job);

    //for (int i=0; i < v.count; i++) {
    //  mpz_mul(s, v.el[i], v.el[i]);
    //  mpz_mod(v.el[i], P.el[i/2], s);
    //}

    free_vec(&P);
#ifdef OUTPUT_REMAINDER_LEVELS
    snprintf(name, sizeof(name)-1, "r%d.mpz", level);
    output_bin_array(&v, name);
#endif
    P = v;
    fprintf(stderr, "%0.3fs\n", now()-start);
  }  

  // final round

  double start = now();
  fprintf(stderr, "output\n");
  input_bin_array(&v, INPUT_FN);

  vec_t w;
  init_vec(&w, v.count);
  void muldiv_job(int i) {
    mpz_t s;
    mpz_init(s);
    mpz_mul(s, v.el[i], v.el[i]);
    mpz_mod(w.el[i], P.el[i/2], s);
    mpz_divexact(w.el[i], w.el[i], v.el[i]);
    mpz_gcd(w.el[i], w.el[i], v.el[i]);
    mpz_clear(s);
  }
  iter_threads(0, v.count, muldiv_job);

  output_bin_array(&w, OUTPUT_FN);
  fprintf(stderr, "%0.3fs\n", now()-start);

  free_vec(&w);
  free_vec(&v);
  free_vec(&P);
  fprintf(stderr, "remainder tree took %0.3fs\n", now()-tstart);
}

void emit_results() {
  vec_t moduli;
  input_bin_array(&moduli, INPUT_FN);

  vec_t gcds;
  input_bin_array(&gcds, OUTPUT_FN);

  double pstart = now();
  fprintf(stderr, "emitting results\n");

  // find elements of w that aren't 1
  int size;
  size = 0;
  for (int i=0; i < gcds.count; i++) {
    if (mpz_cmp_ui(gcds.el[i],1)) {
      mpz_set(moduli.el[size],moduli.el[i]);
      mpz_set(gcds.el[size],gcds.el[i]);
      size++;
    }
  }
  for (int i=size; i < gcds.count; i++) {
    mpz_clear(moduli.el[i]);
    mpz_clear(gcds.el[i]);
  }
  moduli.count = size;
  gcds.count = size;

  output_hex_array(&moduli,"vulnerable_moduli");
  output_hex_array(&gcds,"gcds");

  free_vec(&moduli);
  free_vec(&gcds);

  fprintf(stderr, "emitting %d results took %0.3fs\n", size, now()-pstart);
  fprintf(stderr, "mrow!\n");
}

int main(int argc, char *argv[])
{
  double start = now();

  if (argc != 2) {
    fprintf(stderr, "usage: %s INPUT | --resume\n", argv[0]);
    exit(1);
  }

  if (strcmp(argv[1], "--resume")) {
    char *input = argv[1];
    remove(INPUT_FN);
    prep_hex_input(input, "input.tmp");
    if (rename("input.tmp", INPUT_FN)) {
      fprintf(stderr, "can't rename input temp file\n");
      exit(1);
    }
  } else if (!file_exists(INPUT_FN)) {
    fprintf(stderr, "unable to resume\n");
    exit(1);
  }

  int level = product_tree();
  remainder_tree(level-1);
  emit_results();
  fprintf(stderr, "run took %0.3fs\n", now()-start);
  return 0;
}



