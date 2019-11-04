#ifndef BLS12_TEST_H
#define BLS12_TEST_H

#include <ELiPS/bls12_pairing.h>
#include <ELiPS/G1_SCM.h>
#include <ELiPS/G2_SCM.h>
#include <ELiPS/G3_exp.h>


extern int BLS12_test_rational_point();
extern void BLS12_test_plain_ate_pairing();
extern int BLS12_test_opt_ate_pairing();
extern int BLS12_test_G1_SCM(int scm);
extern int BLS12_test_G2_SCM(int scm);
extern int BLS12_test_G3_EXP(int exp);

/*mine*/
extern int BLS12_test_opt_ate_pairing_projective_projective(int pairing);
extern int BLS12_test_opt_ate_pairing_projective_plain(int pairing);
extern int BLS12_test_opt_ate_pairing_debag(int pairing);
#endif
