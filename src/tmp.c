#include <ELiPS/line_calculate.h>
void fp12_6_sparse_mul_lazy_montgomery(fp12_t *C, fp12_t *A, fp12_t *B) {
    fp6_t tmp0_fp6, tmp1_fp6, tmp2_fp6;
    fp2_t v0, v1, v2, v3, v4;

    fp2_mul_lazy_montgomery(&tmp0_fp6.x0, &A->x0.x0, &B->x0.x0);
    fp2_mul_lazy_montgomery(&tmp0_fp6.x1, &A->x0.x1, &B->x0.x0);//
    fp2_mul_lazy_montgomery(&tmp0_fp6.x2, &A->x0.x2, &B->x0.x0);//
    fp2_add(&tmp2_fp6.x0, &B->x0.x0, &B->x1.x0);//
    fp2_set(&tmp2_fp6.x1, &B->x1.x1);//

    // fp6_6_sparse_mul(&tmp1_fp6, &A->x1, &B->x1);
    fp2_mul_lazy_montgomery(&v0, &A->x1.x0, &B->x1.x0);
    fp2_mul_lazy_montgomery(&v1, &A->x1.x1, &B->x1.x1);
    fp2_add(&v2, &A->x1.x1, &A->x1.x2);
    fp2_mul_lazy_montgomery(&v2, &v2, &B->x1.x1);
    fp2_sub(&v2, &v2, &v1);
    fp2_mul_basis(&v4, &v2);
    fp2_add(&v4, &v4, &v0);
    fp2_add(&v2, &A->x1.x0, &A->x1.x1);
    fp2_add(&v3, &B->x1.x0, &B->x1.x1);
    fp2_mul_lazy_montgomery(&tmp1_fp6.x1, &v2, &v3);
    fp2_sub(&tmp1_fp6.x1, &tmp1_fp6.x1, &v0);
    fp2_sub(&tmp1_fp6.x1, &tmp1_fp6.x1, &v1);
    fp2_add(&v2, &A->x1.x0, &A->x1.x2);
    fp2_mul_lazy_montgomery(&tmp1_fp6.x2, &v2, &B->x1.x0);
    fp2_sub(&tmp1_fp6.x2, &tmp1_fp6.x2, &v0);
    fp2_add(&tmp1_fp6.x2, &tmp1_fp6.x2, &v1);
    fp2_set(&tmp1_fp6.x0, &v4);
    // fp6_sparse_mul end

    fp6_add(&C->x1, &A->x0, &A->x1);

    // fp6_6_sparse_mul(&C->x1, &C->x1, &tmp2_fp6);
    fp2_mul_lazy_montgomery(&v0, &C->x1.x0, &tmp2_fp6.x0);
    fp2_mul_lazy_montgomery(&v1, &C->x1.x1, &tmp2_fp6.x1);
    fp2_add(&v2, &C->x1.x1, &C->x1.x2);
    fp2_mul_lazy_montgomery(&v2, &v2, &tmp2_fp6.x1);
    fp2_sub(&v2, &v2, &v1);
    fp2_mul_basis(&v4, &v2);
    fp2_add(&v4, &v4, &v0);
    fp2_add(&v2, &C->x1.x0, &C->x1.x1);
    fp2_add(&v3, &tmp2_fp6.x0, &tmp2_fp6.x1);
    fp2_mul_lazy_montgomery(&C->x1.x1, &v2, &v3);
    fp2_sub(&C->x1.x1, &C->x1.x1, &v0);
    fp2_sub(&C->x1.x1, &C->x1.x1, &v1);
    fp2_add(&v2, &C->x1.x0, &C->x1.x2);
    fp2_mul_lazy_montgomery(&C->x1.x2, &v2, &tmp2_fp6.x0);
    fp2_sub(&C->x1.x2, &C->x1.x2, &v0);
    fp2_add(&C->x1.x2, &C->x1.x2, &v1);
    fp2_set(&C->x1.x0, &v4);
    // fp6_sparse_mul end

    fp6_sub(&C->x1, &C->x1, &tmp0_fp6);
    fp6_sub(&C->x1, &C->x1, &tmp1_fp6);
    fp6_mul_basis(&tmp1_fp6, &tmp1_fp6);
    fp6_add(&C->x0, &tmp0_fp6, &tmp1_fp6);
}