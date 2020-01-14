#include <ELiPS/fp4.h>
//fp4_t
void fp4_init(fp4_t *A){
    fp2_init(&A->x0);
    fp2_init(&A->x1);
}

void fp4_printf(char *str,fp4_t *A){
    gmp_printf("%s(",str);
    fp2_printf("",&A->x0);
    gmp_printf(",");
    fp2_printf("",&A->x1);
    gmp_printf(")");
}

void fp4_println(char *str,fp4_t *A){
    gmp_printf("%s(",str);
    fp2_printf("",&A->x0);
    gmp_printf(",");
    fp2_printf("",&A->x1);
    gmp_printf(")\n");
}
void fp4_printf_montgomery(char *str,fp4_t *A){
    gmp_printf("%s(",str);
    fp2_printf_montgomery("",&A->x0);
    gmp_printf(",");
    fp2_printf_montgomery("",&A->x1);
    gmp_printf(")");
}
void fp4_set(fp4_t *ANS,fp4_t *A){
    fp2_set(&ANS->x0,&A->x0);
    fp2_set(&ANS->x1,&A->x1);
}
void fp4_set_ui(fp4_t *ANS,unsigned long int UI){
    fp2_set_ui(&ANS->x0,UI);
    fp2_set_ui(&ANS->x1,0);
}
void fp4_set_ui_ui(fp4_t *ANS,unsigned long int UI){
    fp2_set_ui(&ANS->x0,UI);
    fp2_set_ui(&ANS->x1,UI);
}

void fp4_set_mpn(fp4_t *ANS,mp_limb_t *A){
    fp2_set_mpn(&ANS->x0,A);
    fp2_set_ui(&ANS->x1,0);    
}

void fp4_set_neg(fp4_t *ANS,fp4_t *A){
    fp2_set_neg(&ANS->x0,&A->x0);
    fp2_set_neg(&ANS->x1,&A->x1);
}
void fp4_to_montgomery(fp4_t *ANS,fp4_t *A){
    fp2_to_montgomery(&ANS->x0,&A->x0);
    fp2_to_montgomery(&ANS->x1,&A->x1);
}
void fp4_mod_montgomery(fp4_t *ANS,fp4_t *A){
    fp2_mod_montgomery(&ANS->x0,&A->x0);
    fp2_mod_montgomery(&ANS->x1,&A->x1);
}
void fp4_set_random(fp4_t *ANS,gmp_randstate_t state){
    fp2_set_random(&ANS->x0,state);
    fp2_set_random(&ANS->x1,state);
}
void fp4_mul(fp4_t *ANS,fp4_t *A,fp4_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2;
    
    //set
    fp2_mul(&tmp2_fp2,&A->x1,&B->x1);//b*d
    fp2_add(&tmp1_fp2,&A->x0,&A->x1);//a+b
    fp2_add(&ANS->x1,&B->x0,&B->x1);//c+d
    fp2_mul(&ANS->x1,&tmp1_fp2,&ANS->x1);//(a+b)(c+d)
    fp2_mul(&tmp1_fp2,&A->x0,&B->x0);//a*c
    
    //x0
    fp2_mul_basis(&ANS->x0,&tmp2_fp2);//b*d*v
    fp2_add(&ANS->x0,&ANS->x0,&tmp1_fp2);//a*c+b*d*v
    
    //x1
    fp2_sub(&ANS->x1,&ANS->x1,&tmp1_fp2);
    fp2_sub(&ANS->x1,&ANS->x1,&tmp2_fp2);
}
void fp4_mul_lazy(fp4_t *ANS,fp4_t *A,fp4_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2;
    //set
    fp2_mul_lazy(&tmp2_fp2,&A->x1,&B->x1);//b*d
    fp2_add_lazy(&tmp1_fp2,&A->x0,&A->x1);//a+b
    fp2_add_lazy(&ANS->x1,&B->x0,&B->x1);//c+d
    fp2_mul_lazy(&ANS->x1,&tmp1_fp2,&ANS->x1);//(a+b)(c+d)
    fp2_mul_lazy(&tmp1_fp2,&A->x0,&B->x0);//a*c
    
    //x0
    fp2_mul_basis(&ANS->x0,&tmp2_fp2);//b*d*v
    fp2_add(&ANS->x0,&ANS->x0,&tmp1_fp2);//a*c+b*d*v
    
    //x1
    fp2_sub(&ANS->x1,&ANS->x1,&tmp1_fp2);
    fp2_sub(&ANS->x1,&ANS->x1,&tmp2_fp2);
}
void fp4_mul_lazy_montgomery(fp4_t *ANS,fp4_t *A,fp4_t *B){
    static fp2_t tmp1_fp2,tmp2_fp2;
    //set
    fp2_mul_lazy_montgomery(&tmp2_fp2,&A->x1,&B->x1);//b*d
    fp2_add_lazy(&tmp1_fp2,&A->x0,&A->x1);//a+b
    fp2_add_lazy(&ANS->x1,&B->x0,&B->x1);//c+d
    fp2_mul_lazy_montgomery(&ANS->x1,&tmp1_fp2,&ANS->x1);//(a+b)(c+d)
    fp2_mul_lazy_montgomery(&tmp1_fp2,&A->x0,&B->x0);//a*c
    
    //x0
    fp2_mul_basis(&ANS->x0,&tmp2_fp2);//b*d*v
    fp2_add(&ANS->x0,&ANS->x0,&tmp1_fp2);//a*c+b*d*v
    
    //x1
    fp2_sub(&ANS->x1,&ANS->x1,&tmp1_fp2);
    fp2_sub(&ANS->x1,&ANS->x1,&tmp2_fp2);
}
void fp4_mul_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
    fp2_mul_ui(&ANS->x0,&A->x0,UI);
    fp2_mul_ui(&ANS->x1,&A->x1,UI);
}

void fp4_mul_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B){
    fp2_mul_mpn(&ANS->x0,&A->x0,B);
    fp2_mul_mpn(&ANS->x1,&A->x1,B);
}
void fp4_mul_basis(fp4_t *ANS,fp4_t *A){
    static fp4_t tmp1_fp4;
    fp_set(&tmp1_fp4.x1,&A->x0);
    
    fp2_mul_basis(&tmp1_fp4.x0,&A->x1);
    fp4_set(ANS,&tmp1_fp4);
}
void fp4_mul_basis_lazy(fp4_t *ANS,fp4_t *A){
    static fp4_t tmp1_fp4;
    fp_set(&tmp1_fp4.x1,&A->x0);
    
    fp2_mul_basis_lazy(&tmp1_fp4.x0,&A->x1);
    fp4_set(ANS,&tmp1_fp4);
}
void fp4_sqr(fp4_t *ANS,fp4_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    fp2_add(&tmp1_fp2,&A->x0,&A->x1);
    fp2_mul_basis(&tmp2_fp2,&A->x1);
    fp2_add(&tmp2_fp2,&tmp2_fp2,&A->x0);
    fp2_mul(&tmp3_fp2,&A->x0,&A->x1);
	
    //x0
    fp2_mul(&ANS->x0,&tmp1_fp2,&tmp2_fp2);
    fp2_sub(&ANS->x0,&ANS->x0,&tmp3_fp2);
    fp2_mul_basis(&tmp1_fp2,&tmp3_fp2);
    fp2_sub(&ANS->x0,&ANS->x0,&tmp1_fp2);

    //x1
    fp2_add(&ANS->x1,&tmp3_fp2,&tmp3_fp2);
}
void fp4_sqr_lazy(fp4_t *ANS,fp4_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    fp2_add_lazy(&tmp1_fp2,&A->x0,&A->x1);
    fp2_mul_basis_lazy(&tmp2_fp2,&A->x1);
    fp2_add_lazy(&tmp2_fp2,&tmp2_fp2,&A->x0);
    fp2_mul_lazy(&tmp3_fp2,&A->x0,&A->x1);
	
    //x0
    fp2_mul_lazy(&ANS->x0,&tmp1_fp2,&tmp2_fp2);
    fp2_sub(&ANS->x0,&ANS->x0,&tmp3_fp2);
    fp2_mul_basis(&tmp1_fp2,&tmp3_fp2);
    fp2_sub(&ANS->x0,&ANS->x0,&tmp1_fp2);
    
    //x1
    fp2_add(&ANS->x1,&tmp3_fp2,&tmp3_fp2);
}
void fp4_sqr_lazy_montgomery(fp4_t *ANS,fp4_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2;
    fp2_add_lazy(&tmp1_fp2,&A->x0,&A->x1);
    fp2_mul_basis_lazy(&tmp2_fp2,&A->x1);
    fp2_add_lazy(&tmp2_fp2,&tmp2_fp2,&A->x0);
    fp2_mul_lazy_montgomery(&tmp3_fp2,&A->x0,&A->x1);
	
    //x0
    fp2_mul_lazy_montgomery(&ANS->x0,&tmp1_fp2,&tmp2_fp2);
    fp2_sub(&ANS->x0,&ANS->x0,&tmp3_fp2);
    fp2_mul_basis(&tmp1_fp2,&tmp3_fp2);
    fp2_sub(&ANS->x0,&ANS->x0,&tmp1_fp2);
    
    //x1
    fp2_add(&ANS->x1,&tmp3_fp2,&tmp3_fp2);
}
void fp4_add(fp4_t *ANS,fp4_t *A,fp4_t *B){
    fp2_add(&ANS->x0,&A->x0,&B->x0);
    fp2_add(&ANS->x1,&A->x1,&B->x1);
}
void fp4_add_lazy(fp4_t *ANS,fp4_t *A,fp4_t *B){
    fp2_add_lazy(&ANS->x0,&A->x0,&B->x0);
    fp2_add_lazy(&ANS->x1,&A->x1,&B->x1);
}

void fp4_add_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
    fp2_add_ui(&ANS->x0,&A->x0,UI);
    fp2_add_ui(&ANS->x1,&A->x1,0);
}

void fp4_add_ui_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
    fp2_add_ui_ui(&ANS->x0,&A->x0,UI);
    fp2_add_ui_ui(&ANS->x1,&A->x1,UI);
}
void fp4_add_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B){
    fp2_add_mpn(&ANS->x0,&ANS->x0,B);
    fp2_add_mpn(&ANS->x1,&ANS->x1,B);
}

void fp4_sub(fp4_t *ANS,fp4_t *A,fp4_t *B){
    fp2_sub(&ANS->x0,&A->x0,&B->x0);
    fp2_sub(&ANS->x1,&A->x1,&B->x1);
}
void fp4_sub_lazy(fp4_t *ANS,fp4_t *A,fp4_t *B){
    fp2_sub_lazy(&ANS->x0,&A->x0,&B->x0);
    fp2_sub_lazy(&ANS->x1,&A->x1,&B->x1);
}

void fp4_sub_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
    fp2_sub_ui(&ANS->x0,&ANS->x0,UI);
    fp2_sub_ui(&ANS->x1,&ANS->x1,0);
}

void fp4_sub_ui_ui(fp4_t *ANS,fp4_t *A,unsigned long int UI){
    fp2_sub_ui_ui(&ANS->x0,&ANS->x0,UI);
    fp2_sub_ui_ui(&ANS->x1,&ANS->x1,UI);
}
void fp4_sub_mpn(fp4_t *ANS,fp4_t *A,mp_limb_t *B){
    fp2_sub_mpn(&ANS->x0,&ANS->x0,B);
    fp2_sub_mpn(&ANS->x1,&ANS->x1,B);
}

void fp4_inv(fp4_t *ANS,fp4_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2;
    fp2_set(&tmp1_fp2,&A->x0);
    fp2_set_neg(&tmp2_fp2,&A->x1);
    
    fp2_mul(&tmp3_fp2,&tmp1_fp2,&A->x0);
    fp2_mul(&tmp4_fp2,&tmp2_fp2,&A->x1);
    fp2_mul_basis(&tmp4_fp2,&tmp4_fp2);
    fp2_add(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);
    fp2_inv(&tmp3_fp2,&tmp3_fp2);
    fp2_mul(&ANS->x0,&tmp1_fp2,&tmp3_fp2);
    fp2_mul(&ANS->x1,&tmp2_fp2,&tmp3_fp2);
}
void fp4_inv_lazy(fp4_t *ANS,fp4_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2;
    fp2_set(&tmp1_fp2,&A->x0);
    fp2_set_neg(&tmp2_fp2,&A->x1);
    
    fp2_mul_lazy(&tmp3_fp2,&tmp1_fp2,&A->x0);
    fp2_mul_lazy(&tmp4_fp2,&tmp2_fp2,&A->x1);
    fp2_mul_basis(&tmp4_fp2,&tmp4_fp2);
    fp2_add(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);
    fp2_inv_lazy(&tmp3_fp2,&tmp3_fp2);
    fp2_mul_lazy(&ANS->x0,&tmp1_fp2,&tmp3_fp2);
    fp2_mul_lazy(&ANS->x1,&tmp2_fp2,&tmp3_fp2);
}
void fp4_inv_lazy_montgomery(fp4_t *ANS,fp4_t *A){
    static fp2_t tmp1_fp2,tmp2_fp2,tmp3_fp2,tmp4_fp2;
    fp2_set(&tmp1_fp2,&A->x0);
    fp2_set_neg(&tmp2_fp2,&A->x1);
    
    fp2_mul_lazy_montgomery(&tmp3_fp2,&tmp1_fp2,&A->x0);
    fp2_mul_lazy_montgomery(&tmp4_fp2,&tmp2_fp2,&A->x1);
    fp2_mul_basis(&tmp4_fp2,&tmp4_fp2);
    fp2_add(&tmp3_fp2,&tmp3_fp2,&tmp4_fp2);
    fp2_inv_lazy_montgomery(&tmp3_fp2,&tmp3_fp2);
    fp2_mul_lazy_montgomery(&ANS->x0,&tmp1_fp2,&tmp3_fp2);
    fp2_mul_lazy_montgomery(&ANS->x1,&tmp2_fp2,&tmp3_fp2);
}
int  fp4_legendre(fp4_t *A){
    mpz_t exp;
    mpz_init(exp);
    fp4_t tmp;
    fp4_init(&tmp);
    
    mpz_pow_ui(exp,prime_z,4);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp4_pow(&tmp,A,exp);
    
    if(fp4_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

void fp4_sqrt(fp4_t *ANS,fp4_t *A){
    fp4_t x,y,t,k,n,tmp;
    fp4_init(&x);
    fp4_init(&y);
    fp4_init(&t);
    fp4_init(&k);
    fp4_init(&n);
    fp4_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    //gmp_randstate_t state;
    //gmp_randinit_default (state);
    //gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    fp4_set_random(&n,state);
    while(fp4_legendre(&n)!=-1){
        fp4_set_random(&n,state);
    }
    mpz_pow_ui(q,prime_z,12);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    fp4_pow(&y,&n,q);
    mpz_set_ui(z,e);    
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp4_pow(&x,A,exp);
    fp4_mul(&tmp,&x,&x);
    fp4_mul(&k,&tmp,A);
    fp4_mul(&x,&x,A);
    while(fp4_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        fp4_pow(&tmp,&k,exp);
        while(fp4_cmp_one(&tmp)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            fp4_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        fp4_pow(&t,&y,result);
        fp4_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        fp4_mul(&x,&x,&t); 
        fp4_mul(&k,&k,&y);
    }
    fp4_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
}

void fp4_pow(fp4_t *ANS,fp4_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp4_t tmp;
    fp4_init(&tmp);
    fp4_set(&tmp,A);
    
    for(i=1;i<length; i++){
        fp4_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            fp4_mul(&tmp,A,&tmp);
        }
    }
    
    fp4_set(ANS,&tmp);
}

int  fp4_cmp(fp4_t *A,fp4_t *B){
    if(fp2_cmp(&A->x0,&B->x0)==0 && fp2_cmp(&A->x1,&B->x1)==0){
        return 0;   
    }
    return 1;
}

int  fp4_cmp_ui(fp4_t *A,unsigned long int UI){
    if(fp2_cmp_ui(&A->x0,UI)==0 && fp2_cmp_ui(&A->x1,UI)==0){
        return 0;
    }
    return 1;
}

int  fp4_cmp_mpn(fp4_t *A,mp_limb_t *B){
    if(fp2_cmp_mpn(&A->x0,B)==0 && fp2_cmp_mpn(&A->x1,B)==0){
        return 0;
    }
    return 1;
}

int  fp4_cmp_zero(fp4_t *A){
    if(fp2_cmp_zero(&A->x0)==0 && fp2_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}

int  fp4_cmp_one(fp4_t *A){
    if(fp2_cmp_one(&A->x0)==0 && fp2_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
int fp4_montgomery_trick(fp4_t *A_inv,fp4_t *A,int n){
    int i;
    fp4_t ANS[n],ALL_inv;
	fp4_set(&ANS[0],&A[0]);
	fp4_t check;
	
	for(i=1;i<n;i++){
	fp4_mul_lazy(&ANS[i],&ANS[i-1],&A[i]);
	}
	fp4_inv_lazy(&ALL_inv,&ANS[n-1]);	
	for(i=n-1;i>0;i--){
    fp4_mul_lazy(&A_inv[i],&ALL_inv,&ANS[i-1]);	
    fp4_mul_lazy(&ALL_inv,&ALL_inv,&A[i]);
    }
    
    fp4_set(&A_inv[0],&ALL_inv);
    /*
    for(i=0;i<n;i++){
    fp4_mul(&check,&A[i],&A_inv[i]);
    printf("check:%d",i);	
	fp4_println("=",&check);
    }
    */
    return 0;
}
int fp4_montgomery_trick_montgomery(fp4_t *A_inv,fp4_t *A,int n){
    int i;
    fp4_t ANS[n],ALL_inv;
	fp4_set(&ANS[0],&A[0]);
	fp4_t check;
	
	for(i=1;i<n;i++){
	fp4_mul_lazy_montgomery(&ANS[i],&ANS[i-1],&A[i]);
	}
	fp4_inv_lazy_montgomery(&ALL_inv,&ANS[n-1]);	
	for(i=n-1;i>0;i--){
    fp4_mul_lazy_montgomery(&A_inv[i],&ALL_inv,&ANS[i-1]);	
    fp4_mul_lazy_montgomery(&ALL_inv,&ALL_inv,&A[i]);
    }
    
    fp4_set(&A_inv[0],&ALL_inv);
    return 0;
}
