#include <ELiPS/define.h>

int bls12_X_binary[bls12_X_length + 1];
int bls12_X2_binary[bls12_X2_length + 1];

int cost_add, cost_add_ui, cost_sub, cost_sub_ui, cost_mul, cost_mul_ui, cost_sqr, cost_inv, cost_mod;

mp_limb_t buf[FPLIMB], tmp_mul[FPLIMB2], tmp1[FPLIMB], tmp2[FPLIMB];

/*============================================================================*/
/* Pairing functions                                                          */
/*============================================================================*/

gmp_randstate_t state;
mpz_t X_z, prime_z, order_z, trace_z;
mp_limb_t X, prime[FPLIMB], order[FPLIMB], trace[FPLIMB];
mp_limb_t prime2[FPLIMB2];
fp2_t Alpha_1, Alpha_1_inv;
mp_limb_t epsilon1[FPLIMB], epsilon2[FPLIMB];
mp_limb_t Two_inv[FPLIMB];
mpz_t Two_inv_z;
mpz_t root_2, root_X;
mpz_t efp_total, efp12_total;
fp2_t frobenius_constant[12][6];
fp2_t skew_frobenius_constant[12][2];
mp_limb_t curve_b[FPLIMB];

//montgomery
mp_limb_t R[FPLIMB], Ri[FPLIMB], R1[FPLIMB], RR[FPLIMB], Ni[FPLIMB];
int m;
mp_limb_t u[FPLIMB + 1];
mp_limb_t N[FPLIMB2], R2[FPLIMB], R3[FPLIMB], RmodP[FPLIMB];
/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
struct timeval tv_start, tv_end;
float MILLER_TATE, MILLER_PLAINATE, MILLER_OPT, MILLER_XATE;
float FINALEXP_PLAIN, FINALEXP_OPT, FINALEXP_OPT_EASY, FINALEXP_OPT_HARD;
float MILLER_OPT_PROJECTIVE, FINALEXP_OPT_PROJECTIVE;
float MILLER_OPT_MONTGOMERY, FINALEXP_OPT_MONTGOMERY;
float MILLER_OPT_PROJECTIVE_MONTGOMERY, FINALEXP_OPT_PROJECTIVE_MONTGOMERY;
float MILLER_OPT_PROJECTIVE_7SPARSE_MONTGOMERY, FINALEXP_OPT_PROJECTIVE_7SPARSE_MONTGOMERY;

cost MILLER_OPT_COST, FINALEXP_OPT_COST;
cost MILLER_OPT_PROJECTIVE_COST, FINALEXP_OPT_PROJECTIVE_COST;
cost MILLER_OPT_MONTGOMERY_COST, FINALEXP_OPT_MONTGOMERY_COST;
cost MILLER_OPT_PROJECTIVE_MONTGOMERY_COST, FINALEXP_OPT_PROJECTIVE_MONTGOMERY_COST;
cost MILLER_OPT_PROJECTIVE_7SPARSE_MONTGOMERY_COST, FINALEXP_OPT_PROJECTIVE_7SPARSE_MONTGOMERY_COST;

float G1SCM_PLAIN, G1SCM_2SPLIT, G1SCM_2SPLIT_JSF;
float G2SCM_PLAIN, G2SCM_2SPLIT, G2SCM_2SPLIT_JSF, G2SCM_4SPLIT;
float G3EXP_PLAIN, G3EXP_2SPLIT, G3EXP_2SPLIT_JSF, G3EXP_4SPLIT;
