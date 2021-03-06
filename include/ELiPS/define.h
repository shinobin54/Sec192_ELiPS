#ifndef DEFINE_H
#define DEFINE_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

/*************OPTION*************/
//bit
#define X64

//debug
#define DEBUG_COST_A

/************************************/

#ifdef X64
#define FPLIMB_BITS 512
#define FPLIMB 8
#define FPLIMB2 16
#endif

#ifdef X32
#define FPLIMB_BITS 480
#define FPLIMB 15
#define FPLIMB2 30
#endif

#define scalar_t mpz_t


#define bls12_X_length 48
#define bls12_X2_length 47

//#define bls12_X_length 76
//#define bls12_X2_length 75

#define BLS12_G1_SCM bls12_g1_scm_2split_5naf_interleaving_mixture_lazy_montgomery
#define BLS12_G2_SCM bls12_g2_scm_4split_5naf_interleaving_mixture_lazy_montgomery
#define BLS12_G3_EXP bls12_g3_exp_4split_5naf_interleaving_GS_lazy_montgomery
#define BLS12_OPTATE_PAIRING bls12_optate_pairing_compress_lazy_montgomery

#define bls12_g1_scm bls12_g1_scm_2split_5naf_interleaving_mixture_lazy_montgomery
#define bls12_g2_scm bls12_g2_scm_4split_5naf_interleaving_mixture_lazy_montgomery
#define bls12_g3_exp bls12_g3_exp_4split_5naf_interleaving_GS_lazy_montgomery
#define bls12_optate_pairing bls12_optate_pairing_compress_lazy_montgomery

//temporary
#define bls12_g1_scm_plain bls12_g1_scm_basic
#define bls12_g2_scm_plain bls12_g2_scm_basic
#define bls12_g3_exp_plain bls12_g3_exp_basic

extern int bls12_X_binary[bls12_X_length + 1];
extern int bls12_X2_binary[bls12_X2_length + 1];

extern int cost_add, cost_add_ui, cost_sub, cost_sub_ui, cost_mul, cost_mul_ui, cost_sqr, cost_inv, cost_mod;

typedef struct
{
	int add;
	int add_ui;
	int sub;
	int sub_ui;
	int mul;
	int mul_ui;
	int sqr;
	int inv;
	int mod;
} cost;

/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/

typedef struct
{
	mp_limb_t x0[FPLIMB];
} fp_t;
typedef struct
{
	fp_t x0;
	fp_t x1;
} fp2_t;
typedef struct
{
	fp2_t x0;
	fp2_t x1;
} fp4_t;
typedef struct
{
	fp4_t x0, x1, x2;
} fp12_t;
typedef struct
{
	fp12_t x0, x1;
} fp24_t;

//tmp finite field
extern mp_limb_t buf[FPLIMB];

/*============================================================================*/
/* Elliptic Curve                                                             */
/*============================================================================*/
typedef struct
{
	fp_t x, y;
	int infinity;
} efp_t;

typedef struct
{
	fp2_t x, y;
	int infinity;
} efp2_t;

typedef struct
{
	fp4_t x, y;
	int infinity;
} efp4_t;

typedef struct
{
	fp12_t x, y;
	int infinity;
} efp12_t;

typedef struct
{
	fp24_t x, y;
	int infinity;
} efp24_t;

/*============================================================================*/
/* Jacobian Elliptic Curve                                                   */
/*============================================================================*/
typedef struct
{
	fp_t x, y, z;
	int infinity;
} efp_projective_t;

typedef struct
{
	fp2_t x, y, z;
	int infinity;
} efp2_projective_t;

typedef struct
{
	fp4_t x, y, z;
	int infinity;
} efp4_projective_t;

typedef struct
{
	fp12_t x, y, z;
	int infinity;
} efp12_projective_t;

typedef struct
{
	fp24_t x, y, z;
	int infinity;
} efp24_projective_t;

/*============================================================================*/
/* Jacobian Elliptic Curve                                                   */
/*============================================================================*/
typedef struct
{
	fp_t x, y, z;
	int infinity;
} efp_jacobian_t;

typedef struct
{
	fp2_t x, y, z;
	int infinity;
} efp2_jacobian_t;

typedef struct
{
	fp4_t x, y, z;
	int infinity;
} efp4_jacobian_t;

typedef struct
{
	fp12_t x, y, z;
	int infinity;
} efp12_jacobian_t;

typedef struct
{
	fp24_t x, y, z;
	int infinity;
} efp24_jacobian_t;

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/
enum f_state
{
	f_p1,
	f_p2,
	f_p3,
	f_p4,
	f_p5,
	f_p6,
	f_p7,
	f_p8,
	f_p9,
	f_p10,
	f_p11,
	f_p12
};
extern gmp_randstate_t state;
extern mpz_t X_z, prime_z, order_z, trace_z;
extern mp_limb_t X, prime[FPLIMB], order[FPLIMB], trace[FPLIMB];
extern mp_limb_t prime2[FPLIMB2];
extern fp2_t Alpha_1, Alpha_1_inv;
extern mp_limb_t epsilon1[FPLIMB], epsilon2[FPLIMB];
extern mp_limb_t Two_inv[FPLIMB];
extern mpz_t Two_inv_z;
extern mpz_t root_2, root_X;
extern mpz_t efp_total, efp12_total;
extern fp2_t frobenius_constant[12][6];
extern fp2_t skew_frobenius_constant[12][2];
extern mp_limb_t curve_b[FPLIMB];

//montgomery
extern mp_limb_t R[FPLIMB], Ri[FPLIMB], R1[FPLIMB], RR[FPLIMB], Ni[FPLIMB];
extern int m;

extern mp_limb_t u[FPLIMB + 1];
extern mp_limb_t N[FPLIMB2], R2[FPLIMB], R3[FPLIMB], RmodP[FPLIMB];
/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
extern struct timeval tv_start, tv_end;
extern float MILLER_OPT;
extern float FINALEXP_OPT, FINALEXP_OPT_EASY, FINALEXP_OPT_HARD;
extern float MILLER_OPT_PROJECTIVE, FINALEXP_OPT_PROJECTIVE;
extern float MILLER_OPT_MONTGOMERY, FINALEXP_OPT_MONTGOMERY;
extern float MILLER_OPT_PROJECTIVE_MONTGOMERY, FINALEXP_OPT_PROJECTIVE_MONTGOMERY;
extern float MILLER_OPT_PROJECTIVE_7SPARSE_MONTGOMERY, FINALEXP_OPT_PROJECTIVE_7SPARSE_MONTGOMERY;

extern cost MILLER_OPT_COST, FINALEXP_OPT_COST;
extern cost MILLER_OPT_PROJECTIVE_COST, FINALEXP_OPT_PROJECTIVE_COST;
extern cost MILLER_OPT_MONTGOMERY_COST, FINALEXP_OPT_MONTGOMERY_COST;
extern cost MILLER_OPT_PROJECTIVE_MONTGOMERY_COST, FINALEXP_OPT_PROJECTIVE_MONTGOMERY_COST;
extern cost MILLER_OPT_PROJECTIVE_7SPARSE_MONTGOMERY_COST, FINALEXP_OPT_PROJECTIVE_7SPARSE_MONTGOMERY_COST;

#endif
