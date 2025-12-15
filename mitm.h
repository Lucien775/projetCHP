#ifndef MITM_H
#define MITM_H 

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

#define ROTL32(x,r) (((x)<<(r)) | (x>>(32-(r))))
#define ROTR32(x,r) (((x)>>(r)) | ((x)<<(32-(r))))
#define ER32(x,y,k) (x=ROTR32(x,8), x+=y, x^=k, y=ROTL32(y,3), y^=x)
#define DR32(x,y,k) (y^=x, y=ROTR32(y,3), x^=k, x-=y, x=ROTL32(x,8))

typedef uint64_t u64;       /* portable 64-bit integer */
typedef uint32_t u32;       /* portable 32-bit integer */
typedef uint8_t u8;
struct __attribute__ ((packed)) entry { u32 k; u64 v; };  /* hash table entry */
struct __attribute__ ((packed)) pair_zx { u64 z; u64 x;}; /*Pour l'envoi*/
struct __attribute__ ((packed)) z_dest { u64 z; u8 dest;};/*Utile*/

/***************************** global variables ******************************/

extern u64 n;         /* block size (in bits) */
extern u64 mask;          /* this is 2**n - 1 */

extern u64 dict_size;     /* number of slots in the hash table */
extern struct entry *A_local;   /* the hash table */

/* (P, C) : two plaintext-ciphertext pairs */
extern u32 P[2][2];
extern u32 C[2][2];


/******************************** dictionary ********************************/
/*
 * "classic" hash table for 64-bit key-value pairs, with linear probing.  
 * It operates under the assumption that the keys are somewhat random 64-bit integers.
 * The keys are only stored modulo 2**32 - 5 (a prime number), and this can lead 
 * to some false positives.
 */
static const u32 EMPTY = 0xffffffff;
static const u64 PRIME = 0xfffffffb;

/******************************** function ********************************/
double wtime();
u64 murmur64(u64 x);
void human_format(u64 n, char *target);
void Speck64128KeySchedule(const u32 K[],u32 rk[]);
void Speck64128Encrypt(const u32 Pt[], u32 Ct[], const u32 rk[]);
void Speck64128Decrypt(u32 Pt[], const u32 Ct[], u32 const rk[]);
void dict_shard_setup(u64 size);
void dict_insert_sharded(u64 key, u64 value);
int dict_probe(u64 key, int maxval, u64 values[]);
u64 f(u64 k);
u64 g(u64 k);
bool is_good_pair(u64 k1, u64 k2);
int golden_claw_search(int maxres, u64 k1[], u64 k2[]);
void usage(char **argv);
void process_command_line_options(int argc, char ** argv);
int main(int argc, char **argv);

#endif
