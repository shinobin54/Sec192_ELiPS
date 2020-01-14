#ifndef efp4_H
#define efp4_H

#include <ELiPS/efp2.h>

/**
 * @brief Initializes a efp4_t type struct
 *
 * @param[in]P --a pointer to be initialized.
 */
extern void efp4_init(efp4_t *P);

/**
 * @brief Print a efp4_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp4_printf(char *str,efp4_t *P);

/**
 * @brief Print a efp4_t type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]P --a pointer to be printed.
 */
extern void efp4_println(char *str,efp4_t *P);

/**
 * @brief Set a efp4_t type struct to a efp4_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp4_set(efp4_t *ANS,efp4_t *A);

/**
 * @brief Set an unsigned int to a efp4_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]UI1 --an unsigned long int to set.
 * @param[in]UI2 --an unsigned long int to set.
 */
extern void efp4_set_ui(efp4_t *ANS,unsigned long int UI1,unsigned long int UI2);

/**
 * @brief Set a mpn type struct to a efp4_t type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void efp4_set_mpn(efp4_t *ANS,mp_limb_t *A);

/**
 * @brief Negate efp4_t type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void efp4_set_neg(efp4_t *ANS,efp4_t *A);

/**
 * @brief Compare efp4_t type construct and efp4_t type construct
 *
 * @param[in]A --a pointer in efp4_t.
 * @param[in]B --a pointer in efp4_t.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  efp4_cmp(efp4_t *A,efp4_t *B);

/**
 * @brief Generate random rational point.
 *
 * @param[out]P --a pointer in efp4_t.
 */
extern void efp4_rational_point(efp4_t *P);

/**
 * @brief Generate rational point on G1 (BN12).
 *
 * @param[out]P --a pointer in efp4_t.
 */
extern void bn12_generate_g1(efp4_t *P);

/**
 * @brief Doubling a efp4_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp4_t.
 */
extern void efp4_ecd(efp4_t *ANS,efp4_t *P);

/**
 * @brief Doubling a efp4_t type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp4_t.
 */
extern void efp4_ecd_lazy(efp4_t *ANS,efp4_t *P);

/**
 * @brief Addition a efp4_t type struct and a efp4_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp4_t.
 * @param[in]P2 --a pointer in efp4_t.
 */
extern void efp4_eca(efp4_t *ANS,efp4_t *P1,efp4_t *P2);

/**
 * @brief Addition a efp4_t type struct and a efp4_t type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P1 --a pointer in efp4_t.
 * @param[in]P2 --a pointer in efp4_t.
 */
extern void efp4_eca_lazy(efp4_t *ANS,efp4_t *P1,efp4_t *P2);

/**
 * @brief Scalar multiplication a efp4_t type struct
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp4_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void efp4_scm(efp4_t *ANS,efp4_t *P,mpz_t scalar);

/**
 * @brief Scalar multiplication a efp4_t type struct (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]P --a pointer in efp4_t.
 * @param[in]scalar --a pointer in mpz.
 */
extern void efp4_scm_lazy(efp4_t *ANS,efp4_t *P,mpz_t scalar);

#endif

