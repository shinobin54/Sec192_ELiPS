
#include "ELiPS/efp4.h"
#include "ELiPS/bls12_init.h"

/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(void)
{
  bls12_init();
  bls12_print_parameters();

  efp4_t A,B,C;
  
  efp4_init(&A);
  efp4_init(&B);
  efp4_init(&C);

  efp4_rational_point(&A);
  efp4_println("random A=",&A);


  
  return 0;
}
