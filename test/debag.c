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



  printf("debag start\n");
  if(BLS12_test_opt_ate_pairing_debag(2)==1){
    printf("ng\n");
  }


  }
