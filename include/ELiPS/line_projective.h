#ifndef LINE_PROJECTIVE_H
#define LINE_PROJECTIVE_H

#include <ELiPS/EFp12.h>

//Pseudo 8-sparse
extern void FpJ12_init(FpJ12 *A);
extern void FpJ12_set_ui(FpJ12 *A,unsigned int UI);
extern void FpJ12_to_Fp12(Fp12 *ANS,FpJ12 *A);
//extern void FpJ12_clear(FpJ12 *A);
extern void FpJ12_mul(FpJ12 *ANS,FpJ12 *P,FpJ12 *Q);
extern void FpJ12_sqr(FpJ12 *ANS,FpJ12 *A);
extern void ff_ltq_projective_projective(FpJ12 *f,EFpJ12 *T,EFp12 *Q,EFpJ12 *P);
extern void ff_ltt_projective_projective(FpJ12 *f,EFpJ12 *T,EFp12 *Q);

extern void ff_ltq_projective_plain(FpJ12 *f,EFp12 *T,EFp12 *Q,EFp12 *P);
extern void ff_ltt_projective_plain(FpJ12 *f,EFp12 *T,EFp12 *Q);
//extern void ff_ltq_projective(FpJ2 *f,EFpJ2 *T,EFp2 *Q,EFpJ2 *P);
//extern void ff_ltt_projective(FpJ2 *f,EFpJ2 *T,EFp2 *Q);

extern void EFp12_lineTP(Fp12 *rop,EFp12 *q,EFp12 *p,EFp12 *t);
extern void EFp12_lineTT(Fp12 *rop,EFp12 *q,EFp12 *t);


#endif
