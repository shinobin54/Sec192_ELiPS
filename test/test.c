#include "ELiPS/bn12.h"
#include "ELiPS/bls12.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void)
{
  //BN12_init();
  //BN12_print_parameters();
  bls12_init();
  bls12_print_parameters();

  //test_Field(100000,10000,10000,10000,10);
  //test_Frobenius_map();
  //test_skew_frobenius_map();
  //test_twist();
  //test_mod(10000000,100000);

  //test_EFp(10,10,10);
  //test_EFp2(10,10,10);
  //test_EFp12(10,10,10);

  //BN12_test_rational_point();
  //BN12_test_plain_ate_pairing();
  //BN12_test_opt_ate_pairing(10);
  //BN12_test_x_ate_pairing();
  //BN12_test_G1_SCM(100);
  //BN12_test_G2_SCM(10);
  //BN12_test_G3_EXP(10);

  //bls12_test_rational_point();
  //bls12_test_plain_ate_pairing();
  bls12_test_opt_ate_pairing(100);
  //bls12_test_x_ate_pairing();
  //bls12_test_G1_SCM(100);
  //bls12_test_G2_SCM(10);
  //bls12_test_G3_EXP(10);

  //Fast_test_bls12(100,100,100,100);
  //test_All();
  return 0;
}
