#include "ELiPS/bn12.h"
#include "ELiPS/bls12.h"
/*============================================================================*/
/* main                                                                       */
/*============================================================================*/

int main(void){
  int error;
  BLS12_init();
  BLS12_print_parameters();

  EFp12 P,Q,T;
  EFp12_init(&P);
  EFp12_init(&Q);
  EFp12_init(&T);

  EFpJ12 JP,JT;
  EFpJ12_init(&JP);
  EFpJ12_init(&JT);

  Fp12 Z,Z2;
  Fp12_init(&Z);
    Fp12_init(&Z2);

  FpJ12 JZ,JZ2;
  FpJ12_init(&JZ);
  FpJ12_init(&JZ2);



  gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
/*
  BLS12_EFp12_generate_G1(&P);
  EFp12_generate_G2(&Q);
  EFp12_to_EFpJ12(&JP,&P);
  EFp12_to_EFpJ12(&JT,&Q);
  EFp12_set(&T,&Q);

  EFp12_lineTT(&Z2,&Q,&T);
  ff_ltt_projective_Fp12(&JZ2,&JT,&Q);

  EFp12_lineTP(&Z,&Q,&P,&T);
  ff_ltq_projective_Fp12(&JZ,&JT,&Q,&JP);


  Fp12_inv(&JZ.z,&JZ.z);
  Fp12_mul(&JZ.x,&JZ.x,&JZ.z);
  Fp12_inv(&JZ2.z,&JZ2.z);
  Fp12_mul(&JZ2.x,&JZ2.x,&JZ2.z);
  if(Fp12_cmp(&Z,&JZ.x)==0) printf("ADD:ok\n");
  if(Fp12_cmp(&Z2,&JZ2.x)==0) printf("DBL:ok\n");
*/
  printf("homo start\n");
  if(BLS12_test_opt_ate_pairing_projective_plain(2)==1){
    printf("ng\n");
  }
  }
