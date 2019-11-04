#include "ELiPS/bn12.h"
#include "ELiPS/bls12.h"
int main(){
  BLS12_init();
  BLS12_print_parameters();
  Fp_two_inv_set();

  Fp check,rop,rop2;
  Fp_init(&check);
  Fp_init(&rop);
  Fp_init(&rop2);

  Fp_set_ui(&check,11111111111);
  Fp_println("c=",&check);

  Fp_lshift2(&rop2,&check);
  Fp_println("c=",&rop2);

  Fp_rshift2(&rop,&rop2);
  Fp_println("(1/2)*c=",&rop);



  return 0;
}
