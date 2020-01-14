#ifndef fp4_H
#define fp4_H

#include <ELiPS/fp2.h>


/**
 * @brief Initializes a fp4 type struct
 *
 * @param[in]A --a pointer to be initialized.
 */
extern void fp4_init(fp4_t *A);

/**
 * @brief Print a fp4 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp4_printf(char *str,fp4_t *A);

/**
 * @brief Print a fp4 type struct
 *
 * @param[in]str --a pointer to be printed.
 * @param[in]A --a pointer to be printed.
 */
extern void fp4_println(char *str,fp4_t *A);
extern void fp4_printf_montgomery(char *str,fp4_t *A);
/**
 * @brief Set a fp4 type struct to a fp4 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp4_set(fp4_t *ANS,fp4_t *A);

/**
 * @brief Set an unsigned int to a fp4 type struct (x0=UI,x1=0,x2=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp4_set_ui(fp4_t *ANS,unsigned long int UI);

/**
 * @brief Set an unsigned int to a fp4 type struct (x0=UI,x1=0,x2=0)
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --an unsigned long int to set.
 */
extern void fp4_set_ui_ui(fp4_t *ANS,unsigned long int UI);

/**
 * @brief Set a mpn type struct to a fp4 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a pointer to set.
 */
extern void fp4_set_mpn(fp4_t *ANS,mp_limb_t *A);

/**
 * @brief Negate fp4 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be negated.
 */
extern void fp4_set_neg(fp4_t *ANS,fp4_t *A);

extern void fp4_to_montgomery(fp4_t *ANS,fp4_t *A);
extern void fp4_mod_montgomery(fp4_t *ANS,fp4_t *A);
/**
 * @brief Set a random number to a fp4 type struct
 *
 * @param[out]ANS --a pointer to be setted.
 * @param[in]A --a random seed.
 */
extern void fp4_set_random(fp4_t *ANS,gmp_randstate_t state);

/**
 * @brief Multiplication a fp4 type struct and a fp4 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 */
extern void fp4_mul(fp4_t *ANS,fp4_t *A,fp4_t *B);

/**
 * @brief Multiplication a fp4 type struct and a fp4 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 */
extern void fp4_mul_lazy(fp4_t *ANS,fp4_t *A,fp4_t *B);
extern void fp4_mul_lazy_montgomery(fp4_t *ANS,fp4_t *A,fp4_t *B);

/**
 * @brief Multiplication a fp4 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --an unsigned long int.
 */
extern void fp4_mul_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);

/**
 * @brief Multiplication a fp4 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in mpn.
 */
extern void fp4_mul_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B);

/**
 * @brief Squaring a fp4 type struct and a fp4 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 */
extern void fp4_sqr(fp4_t *ANS,fp4_t *A);

/**
 * @brief Squaring a fp4 type struct and a fp4 type struct on prime field (Lazy Reduciton)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 */
extern void fp4_sqr_lazy(fp4_t *ANS,fp4_t *A);
extern void fp4_sqr_lazy_montgomery(fp4_t *ANS,fp4_t *A);

/**
 * @brief Addition a fp4 type struct and a fp4 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 */
extern void fp4_add(fp4_t *ANS,fp4_t *A,fp4_t *B);

/**
 * @brief Addition a fp4 type struct and a fp4 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 */
extern void fp4_add_lazy(fp4_t *ANS,fp4_t *A,fp4_t *B);

/**
 * @brief Addition a fp4 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --an unsigned long int.
 */
extern void fp4_add_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);

/**
 * @brief Addition a fp4 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --an unsigned long int.
 */
extern void fp4_add_ui_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);

/**
 * @brief Addition a fp4 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in mpn.
 */
extern void fp4_add_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B);

/**
 * @brief Subtraction a fp4 type struct and a fp4 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 */
extern void fp4_sub(fp4_t *ANS,fp4_t *A,fp4_t *B);

/**
 * @brief Subtraction a fp4 type struct and a fp4 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 */
extern void fp4_sub_lazy(fp4_t *ANS,fp4_t *A,fp4_t *B);

/**
 * @brief Subtraction a fp4 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --an unsigned long int.
 */
extern void fp4_sub_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp4 type struct and an unsigned long int on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --an unsigned long int.
 */
extern void fp4_sub_ui_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI);

/**
 * @brief Subtraction a fp4 type struct and a mpn type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in mpn.
 */
extern void fp4_sub_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B);

/**
 * @brief Invert a fp4 type struct on prime field
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp4_inv(fp4_t *ANS,fp4_t *A);

/**
 * @brief Invert a fp4 type struct on prime field (Lazy Reduction)
 *
 * @param[out]ANS --a pointer of answer.
 * @param[in]A --a pointer to be inverted.
 */
extern void fp4_inv_lazy(fp4_t *ANS,fp4_t *A);
extern void fp4_inv_lazy_montgomery(fp4_t *ANS,fp4_t *A);

/**
 * @brief LegendreSymbol on prime field
 *
 * @param[in]A --a pointer in fp4.
 * 
 * @return int --a LegendreSymbol (0 or 1 or -1)
 */
extern int  fp4_legendre(fp4_t *A);

/**
 * @brief Whether A is a Cubic non residure on prime field
 *
 * @param[in]A --a pointer in fp4.
 * 
 * @return int --a CNR (0 or 1 or -1)
 */
extern int  fp4_isCNR(fp4_t *A);

/**
 * @brief Sqrt on prime field
 *
 * @param[in]A --a pointer in fp4.
 * @param[out]ANS --a pointer of answer.
 */
extern void fp4_sqrt(fp4_t *ANS,fp4_t *A);

/**
 * @brief Power A by mpz type struct
 *
 * @param[in]scalar --a pointer in fp4.
 * @param[in]A --a pointer in fp4.
 * @param[out]ANS --a pointer in fp4.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern void fp4_pow(fp4_t *ANS,fp4_t *A,mpz_t scalar);

/**
 * @brief Compare fp4 type construct and fp4 type construct
 *
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in fp4.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  fp4_cmp(fp4_t *A,fp4_t *B);

/**
 * @brief Compare fp4 type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp4.
 * @param[in]UI --an unsigned long int.
 * 
 * @return int --(A=UI 0 or other 1)
 */
extern int  fp4_cmp_ui(fp4_t *A,unsigned long int UI);

/**
 * @brief Compare fp4 type construct and mpn type construct
 *
 * @param[in]A --a pointer in fp4.
 * @param[in]B --a pointer in mpn.
 * 
 * @return int --(A=B 0 or other 1)
 */
extern int  fp4_cmp_mpn(fp4_t *A,mp_limb_t *B);

/**
 * @brief Compare fp4 type struct and zero
 *
 * @param[in]A --a pointer in fp4.
 * 
 * @return int --(one 0 or other 1)
 */
extern int  fp4_cmp_zero(fp4_t *A);

/**
 * @brief Compare fp4 type struct and one
 *
 * @param[in]A --a pointer in fp4.
 * 
 * @return int --(zero 0 or other 1)
 */
extern int  fp4_cmp_one(fp4_t *A);

extern int fp4_montgomery_trick(fp4_t *A_inv,fp4_t *A,int n);

extern void fp4_mul_basis(fp4_t *ANS,fp4_t *A);
extern void fp4_mul_basis_lazy(fp4_t *ANS,fp4_t *A);
#endif
