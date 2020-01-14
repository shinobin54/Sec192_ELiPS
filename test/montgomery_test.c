#include "ELiPS/efp4.h"
#include "ELiPS/bls12_init.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
void fp_montgomery_test(){
  fp_t a,b,c,am,bm,cm;
  fp_init(&a);
  fp_init(&b);
  fp_init(&c);
  fp_init(&am);
  fp_init(&bm);
  fp_init(&cm);

  fp_set_random(&a,state);
  fp_set_random(&b,state);
  fp_to_montgomery(&am,&a);
  fp_to_montgomery(&bm,&b);

  fp_mul(&c,&a,&b);
  fp_mulmod_montgomery(&cm,&am,&bm);

  fp_println("c=",&c);
  fp_printf_montgomery("cm=",&cm);
}
void fp4_montgomery_test(){
  fp4_t a,b,c,am,bm,cm;
  fp4_init(&a);
  fp4_init(&b);
  fp4_init(&c);
  fp4_init(&am);
  fp4_init(&bm);
  fp4_init(&cm);

  fp4_set_random(&a,state);
  fp4_set_random(&b,state);
  fp4_to_montgomery(&am,&a);
  fp4_to_montgomery(&bm,&b);

  fp4_mul(&c,&a,&b);
  fp4_mul_lazy_montgomery(&cm,&am,&bm);

  fp4_println("c=",&c);
  fp4_printf_montgomery("cm=",&cm);
}
int main(void)
{
  bls12_init();
  bls12_print_parameters();
  //fp_montgomery_test();
  fp4_montgomery_test();
  return 0;
}
