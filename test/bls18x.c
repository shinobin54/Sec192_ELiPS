#include "ELiPS/bn12.h"
#include "ELiPS/bls12.h"
#include "ELiPS/hello.h"
#include "ELiPS/bls18_init.h"
/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
void mpz_print(mpz_t a){
    mpz_out_str(stdout,10,a);
}
int lenx;
void search_4(){
    int b[4],bitcnt,v1;

    for(b[0]=0;b[0]<lenx;b[0]++){
        for(b[1]=0;b[1]<lenx;b[1]++){
            for(b[2]=0;b[2]<lenx;b[2]++){
                for(b[3]=0;b[3]<lenx;b[3]++){
                mpz_set_ui(X_z,0);
                mpz_setbit(X_z,b[0]);
                mpz_setbit(X_z,b[1]);
                mpz_setbit(X_z,b[2]);
                mpz_setbit(X_z,lenx);
                mpz_setbit(X_z,b[3]);
                if(BLS18_generate_prime()==0) continue;
                if(BLS18_generate_order()==0) continue;

                bitcnt=(int)mpz_sizeinbase(prime_z,2);
                printf("bit=%d\n",bitcnt);
                printf("x=");
                mpz_print(X_z);
                printf(",");
                printf("prime_z=");
                mpz_print(prime_z);
                printf(",");
                printf("order_z=");
                mpz_print(order_z);
                v1=mpz_popcount(X_z);
                printf("popcount=%d\n",v1);
                printf("\n");
                }
            }
        }
    }
}
int main(void){
    mpz_init(X_z);
    mpz_init(prime_z);
    mpz_init(order_z);
    lenx=80;
    search_4();


  }
