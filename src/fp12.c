#include <ELiPS/fp12.h>
//fp12_t
void fp12_init(fp12_t *A){
    fp4_init(&A->x0);
    fp4_init(&A->x1);
    fp4_init(&A->x2);
}

void fp12_printf(char *str,fp12_t *A){
    gmp_printf("%s(",str);
    fp4_printf("",&A->x0);
    gmp_printf(",");
    fp4_printf("",&A->x1);
    gmp_printf(",");
    fp4_printf("",&A->x2);
    gmp_printf(")");
}

void fp12_println(char *str,fp12_t *A){
    gmp_printf("%s(",str);
    fp4_printf("",&A->x0);
    gmp_printf(",");
    fp4_printf("",&A->x1);
    gmp_printf(",");
    fp4_printf("",&A->x2);
    gmp_printf(")\n");
}
void fp12_printf_montgomery(char *str,fp12_t *A){
    gmp_printf("%s(",str);
    fp4_printf_montgomery("",&A->x0);
    gmp_printf(",");
    fp4_printf_montgomery("",&A->x1);
    gmp_printf(",");
    fp4_printf_montgomery("",&A->x2);
    gmp_printf(")");
}
void fp12_set(fp12_t *ANS,fp12_t *A){
    fp4_set(&ANS->x0,&A->x0);
    fp4_set(&ANS->x1,&A->x1);
    fp4_set(&ANS->x2,&A->x2);
}

void fp12_set_ui(fp12_t *ANS,unsigned long int UI){
    fp4_set_ui(&ANS->x0,UI);
    fp4_set_ui(&ANS->x1,0);
    fp4_set_ui(&ANS->x2,0);
}

void fp12_set_ui_ui(fp12_t *ANS,unsigned long int UI){
    fp4_set_ui(&ANS->x0,UI);
    fp4_set_ui(&ANS->x1,UI);
    fp4_set_ui(&ANS->x2,UI);
}
void fp12_set_mpn(fp12_t *ANS,mp_limb_t *A){
    fp4_set_mpn(&ANS->x0,A);
    fp4_set_ui(&ANS->x1,0);
    fp4_set_ui(&ANS->x2,0);
}

void fp12_set_neg(fp12_t *ANS,fp12_t *A){
    fp4_set_neg(&ANS->x0,&A->x0);
    fp4_set_neg(&ANS->x1,&A->x1);
    fp4_set_neg(&ANS->x2,&A->x2);
}
void fp12_to_montgomery(fp12_t *ANS,fp12_t *A){
    fp4_to_montgomery(&ANS->x0,&A->x0);
    fp4_to_montgomery(&ANS->x1,&A->x1);
    fp4_to_montgomery(&ANS->x2,&A->x2);
}
void fp12_mod_montgomery(fp12_t *ANS,fp12_t *A){
    fp4_mod_montgomery(&ANS->x0,&A->x0);
    fp4_mod_montgomery(&ANS->x1,&A->x1);
    fp4_mod_montgomery(&ANS->x2,&A->x2);
}
void fp12_set_random(fp12_t *ANS,gmp_randstate_t state){
    fp4_set_random(&ANS->x0,state);
    fp4_set_random(&ANS->x1,state);
    fp4_set_random(&ANS->x2,state);
}

void fp12_mul(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4,tmp6_fp4,tmp7_fp4;
    //set
    fp4_mul(&tmp1_fp4,&A->x0,&B->x0);//x0*y0
    fp4_mul(&tmp2_fp4,&A->x1,&B->x1);//x1*y1
    fp4_mul(&tmp3_fp4,&A->x2,&B->x2);//x2*y2
    
    fp4_add(&tmp5_fp4,&A->x0,&A->x1);//x0+x1
    fp4_add(&tmp4_fp4,&B->x0,&B->x1);//y0+y1
    fp4_mul(&tmp5_fp4,&tmp5_fp4,&tmp4_fp4);//(x0+x1)(y0+y1)
    
    fp4_add(&tmp6_fp4,&A->x1,&A->x2);//x1+x2
    fp4_add(&tmp4_fp4,&B->x1,&B->x2);//y1+y2
    fp4_mul(&tmp6_fp4,&tmp6_fp4,&tmp4_fp4);//(x1+x2)(y1+y2)
    
    fp4_add(&tmp7_fp4,&B->x0,&B->x2);//y2+y0
    fp4_add(&tmp4_fp4,&A->x0,&A->x2);//x2+x0
    fp4_mul(&tmp7_fp4,&tmp7_fp4,&tmp4_fp4);//(x2+x0)(y2+y0)
    //x0
    fp4_sub(&tmp6_fp4,&tmp6_fp4,&tmp2_fp4);
    fp4_sub(&tmp6_fp4,&tmp6_fp4,&tmp3_fp4);//(x1+x2)(y1+y2)-x1y1-x2y2
    fp4_mul_basis(&tmp4_fp4,&tmp6_fp4);
    fp4_add(&ANS->x0,&tmp1_fp4,&tmp4_fp4);
    //x1
    fp4_sub(&tmp5_fp4,&tmp5_fp4,&tmp1_fp4);
    fp4_sub(&tmp5_fp4,&tmp5_fp4,&tmp2_fp4);
    fp4_mul_basis(&tmp4_fp4,&tmp3_fp4);
    fp4_add(&ANS->x1,&tmp4_fp4,&tmp5_fp4);
    //x2
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp1_fp4);
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp3_fp4);
    fp4_add(&ANS->x2,&tmp2_fp4,&tmp7_fp4);
}
/*
//karat
void fp12_mul(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp4_t v0,v1,v2,a0a1,a0a2,a1a2,b0b1,b0b2,b1b2,buf;
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4,tmp6_fp4,tmp7_fp4;
    //set
    fp4_mul(&v0,&A->x0,&B->x0);//x0*y0
    fp4_mul(&v1,&A->x1,&B->x1);//x1*y1
    fp4_mul(&v2,&A->x2,&B->x2);//x2*y2
    
    fp4_add(&a0a1,&A->x0,&A->x1);//x0+x1
    fp4_add(&b0b1,&B->x0,&B->x1);//y0+y1
    fp4_mul(&tmp1_fp4,&a0a1,&b0b1);//(x0+x1)(y0+y1)
    
    fp4_add(&a1a2,&A->x1,&A->x2);//x1+x2
    fp4_add(&b1b2,&B->x1,&B->x2);//y1+y2
    fp4_mul(&tmp2_fp4,a1a2,b1b2);//(x1+x2)(y1+y2)
    
    fp4_add(&a0a2,&B->x0,&B->x2);//y2+y0
    fp4_add(&b0b2,&A->x0,&A->x2);//x2+x0
    fp4_mul(&tmp3_fp4,&a0a2,&b0b2);//(x2+x0)(y2+y0)
    
    //x0
    fp4_sub(&buf,&tmp2_fp4,&v1);
    fp4_sub(&buf,&buf,&v2);//(x1+x2)(y1+y2)-x1y1-x2y2
    fp4_mul_basis(&buf,&buf);
    fp4_add(&ANS->x0,&v0,&buf);
    
    //mada
    //x1
    fp4_sub(&buf,&tmp1_fp4,&tmp1_fp4);
    fp4_sub(&buf,&tmp5_fp4,&tmp2_fp4);
    fp4_mul_basis(&tmp4_fp4,&tmp3_fp4);
    fp4_add(&ANS->x1,&tmp4_fp4,&tmp5_fp4);
    //x2
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp1_fp4);
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp3_fp4);
    fp4_add(&ANS->x2,&tmp2_fp4,&tmp7_fp4);
}
*/
void fp12_mul_lazy(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4,tmp6_fp4,tmp7_fp4;
    //set
    fp4_mul_lazy(&tmp1_fp4,&A->x0,&B->x0);//tmp1
    fp4_mul_lazy(&tmp2_fp4,&A->x1,&B->x1);//tmp2
    fp4_mul_lazy(&tmp3_fp4,&A->x2,&B->x2);//tmp3
    
    fp4_add_lazy(&tmp5_fp4,&A->x0,&A->x1);
    fp4_add_lazy(&tmp4_fp4,&B->x0,&B->x1);
    fp4_mul_lazy(&tmp5_fp4,&tmp5_fp4,&tmp4_fp4);//tmp5

    fp4_add_lazy(&tmp6_fp4,&A->x1,&A->x2);
    fp4_add_lazy(&tmp4_fp4,&B->x1,&B->x2);
    fp4_mul_lazy(&tmp6_fp4,&tmp6_fp4,&tmp4_fp4);//tmp6
    
    fp4_add_lazy(&tmp7_fp4,&B->x0,&B->x2);
    fp4_add_lazy(&tmp4_fp4,&A->x0,&A->x2);
    fp4_mul_lazy(&tmp7_fp4,&tmp7_fp4,&tmp4_fp4);//tmp7
    //x0
    fp4_sub(&tmp6_fp4,&tmp6_fp4,&tmp2_fp4);
    fp4_sub(&tmp6_fp4,&tmp6_fp4,&tmp3_fp4);
    fp4_mul_basis(&tmp4_fp4,&tmp6_fp4);
    fp4_add(&ANS->x0,&tmp1_fp4,&tmp4_fp4);
    //x1
    fp4_sub(&tmp5_fp4,&tmp5_fp4,&tmp1_fp4);
    fp4_sub(&tmp5_fp4,&tmp5_fp4,&tmp2_fp4);
    fp4_mul_basis(&tmp4_fp4,&tmp3_fp4);
    fp4_add(&ANS->x1,&tmp4_fp4,&tmp5_fp4);
    //x2
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp1_fp4);
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp3_fp4);
    fp4_add(&ANS->x2,&tmp2_fp4,&tmp7_fp4);
}
void fp12_mul_lazy_montgomery(fp12_t *ANS,fp12_t *A,fp12_t *B){
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4,tmp6_fp4,tmp7_fp4;
    //set
    fp4_mul_lazy_montgomery(&tmp1_fp4,&A->x0,&B->x0);//tmp1
    fp4_mul_lazy_montgomery(&tmp2_fp4,&A->x1,&B->x1);//tmp2
    fp4_mul_lazy_montgomery(&tmp3_fp4,&A->x2,&B->x2);//tmp3
    
    fp4_add_lazy(&tmp5_fp4,&A->x0,&A->x1);
    fp4_add_lazy(&tmp4_fp4,&B->x0,&B->x1);
    fp4_mul_lazy_montgomery(&tmp5_fp4,&tmp5_fp4,&tmp4_fp4);//tmp5

    fp4_add_lazy(&tmp6_fp4,&A->x1,&A->x2);
    fp4_add_lazy(&tmp4_fp4,&B->x1,&B->x2);
    fp4_mul_lazy_montgomery(&tmp6_fp4,&tmp6_fp4,&tmp4_fp4);//tmp6
    
    fp4_add_lazy(&tmp7_fp4,&B->x0,&B->x2);
    fp4_add_lazy(&tmp4_fp4,&A->x0,&A->x2);
    fp4_mul_lazy_montgomery(&tmp7_fp4,&tmp7_fp4,&tmp4_fp4);//tmp7
    //x0
    fp4_sub(&tmp6_fp4,&tmp6_fp4,&tmp2_fp4);
    fp4_sub(&tmp6_fp4,&tmp6_fp4,&tmp3_fp4);
    fp4_mul_basis(&tmp4_fp4,&tmp6_fp4);
    fp4_add(&ANS->x0,&tmp1_fp4,&tmp4_fp4);
    //x1
    fp4_sub(&tmp5_fp4,&tmp5_fp4,&tmp1_fp4);
    fp4_sub(&tmp5_fp4,&tmp5_fp4,&tmp2_fp4);
    fp4_mul_basis(&tmp4_fp4,&tmp3_fp4);
    fp4_add(&ANS->x1,&tmp4_fp4,&tmp5_fp4);
    //x2
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp1_fp4);
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp3_fp4);
    fp4_add(&ANS->x2,&tmp2_fp4,&tmp7_fp4);
}
void fp12_mul_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp4_mul_ui(&ANS->x0,&A->x0,UI);
    fp4_mul_ui(&ANS->x1,&A->x1,UI);
    fp4_mul_ui(&ANS->x2,&A->x2,UI);
}

void fp12_mul_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B){
    fp4_mul_mpn(&ANS->x0,&A->x0,B);
    fp4_mul_mpn(&ANS->x1,&A->x1,B);
    fp4_mul_mpn(&ANS->x2,&A->x2,B);
}

void fp12_mul_basis(fp12_t *ANS,fp12_t *A){
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4;
    fp4_set(&tmp1_fp4,&A->x0);
    fp4_set(&tmp2_fp4,&A->x1);
    fp4_set(&tmp3_fp4,&A->x2);
	
    fp2_sub(&ANS->x0.x0,&tmp3_fp4.x0,&tmp3_fp4.x1);
    fp2_add(&ANS->x0.x1,&tmp3_fp4.x0,&tmp3_fp4.x1);
    fp2_set(&ANS->x1.x0,&tmp1_fp4.x0);
    fp2_set(&ANS->x1.x1,&tmp1_fp4.x1);
    fp2_set(&ANS->x2.x0,&tmp2_fp4.x0);
    fp2_set(&ANS->x2.x1,&tmp2_fp4.x1);
}

void fp12_mul_basis_lazy(fp12_t *ANS,fp12_t *A){
    static mp_limb_t tmp1[FPLIMB],tmp2[FPLIMB];
    static fp4_t tmp1_fp4,tmp2_fp4;
    fp4_set(&tmp1_fp4,&A->x0);
    fp4_set(&tmp2_fp4,&A->x1);
    mpn_copyd(tmp1,A->x2.x0.x0,FPLIMB);
    mpn_copyd(tmp2,A->x2.x1.x0,FPLIMB);
	
    fp2_sub_lazy(ANS->x0.x0.x0,FPLIMB,tmp1,FPLIMB,tmp2,FPLIMB);
    fp2_add_lazy(ANS->x0.x1.x0,FPLIMB,tmp1,FPLIMB,tmp2,FPLIMB);

    fp2_set(&ANS->x1.x0,&tmp1_fp4.x0);
    fp2_set(&ANS->x1.x1,&tmp1_fp4.x1);
    fp2_set(&ANS->x2.x0,&tmp2_fp4.x0);
    fp2_set(&ANS->x2.x1,&tmp2_fp4.x1);
}
void fp12_sqr(fp12_t *ANS,fp12_t *A){
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4;
    fp4_sqr(&tmp1_fp4,&A->x0);        //x0^2
    fp4_sqr(&tmp4_fp4,&A->x2);        //x2^2
    fp4_add(&tmp5_fp4,&A->x1,&A->x1);        //2x1
    fp4_mul(&tmp2_fp4,&tmp5_fp4,&A->x2);  //2x1x2
    fp4_mul(&tmp3_fp4,&A->x0,&tmp5_fp4);  //2x0x1
    fp4_add(&tmp5_fp4,&A->x0,&A->x1);        //x0+x1+x2
    fp4_add(&tmp5_fp4,&tmp5_fp4,&A->x2);
    
    //x0
    fp4_mul_basis(&ANS->x0,&tmp2_fp4);
    fp4_add(&ANS->x0,&ANS->x0,&tmp1_fp4);
    //x1
    fp4_mul_basis(&ANS->x1,&tmp4_fp4);
    fp4_add(&ANS->x1,&ANS->x1,&tmp3_fp4);
    //x2
    fp4_sqr(&ANS->x2,&tmp5_fp4);
    fp4_add(&tmp5_fp4,&tmp1_fp4,&tmp4_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp2_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp3_fp4);
    fp4_sub(&ANS->x2,&ANS->x2,&tmp5_fp4);
}
void fp12_sqr_lazy(fp12_t *ANS,fp12_t *A){
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4;
    fp4_sqr_lazy(&tmp1_fp4,&A->x0);        //x0^2
    fp4_sqr_lazy(&tmp4_fp4,&A->x2);        //x2^2
    fp4_add_lazy(&tmp5_fp4,&A->x1,&A->x1);        //2x1

    fp4_mul_lazy(&tmp2_fp4,&tmp5_fp4,&A->x2);  //2x1x2
    fp4_mul_lazy(&tmp3_fp4,&A->x0,&tmp5_fp4);  //2x0x1
    fp4_add(&tmp5_fp4,&A->x0,&A->x1);        //x0+x1+x2
    fp4_add(&tmp5_fp4,&tmp5_fp4,&A->x2);
    
    //x0
    fp4_mul_basis(&ANS->x0,&tmp2_fp4);
    fp4_add(&ANS->x0,&ANS->x0,&tmp1_fp4);

    //x1
    fp4_mul_basis(&ANS->x1,&tmp4_fp4);
    fp4_add(&ANS->x1,&ANS->x1,&tmp3_fp4);

    //x2
    fp4_sqr_lazy(&ANS->x2,&tmp5_fp4);
    fp4_add(&tmp5_fp4,&tmp1_fp4,&tmp4_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp2_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp3_fp4);
    fp4_sub(&ANS->x2,&ANS->x2,&tmp5_fp4);
}
void fp12_sqr_lazy_montgomery(fp12_t *ANS,fp12_t *A){
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4;
    fp4_sqr_lazy(&tmp1_fp4,&A->x0);        //x0^2
    fp4_sqr_lazy(&tmp4_fp4,&A->x2);        //x2^2
    fp4_add_lazy(&tmp5_fp4,&A->x1,&A->x1);        //2x1

    fp4_mul_lazy_montgomery(&tmp2_fp4,&tmp5_fp4,&A->x2);  //2x1x2
    fp4_mul_lazy_montgomery(&tmp3_fp4,&A->x0,&tmp5_fp4);  //2x0x1
    fp4_add(&tmp5_fp4,&A->x0,&A->x1);        //x0+x1+x2
    fp4_add(&tmp5_fp4,&tmp5_fp4,&A->x2);
    
    //x0
    fp4_mul_basis(&ANS->x0,&tmp2_fp4);
    fp4_add(&ANS->x0,&ANS->x0,&tmp1_fp4);

    //x1
    fp4_mul_basis(&ANS->x1,&tmp4_fp4);
    fp4_add(&ANS->x1,&ANS->x1,&tmp3_fp4);

    //x2
    fp4_sqr_lazy(&ANS->x2,&tmp5_fp4);
    fp4_add(&tmp5_fp4,&tmp1_fp4,&tmp4_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp2_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp3_fp4);
    fp4_sub(&ANS->x2,&ANS->x2,&tmp5_fp4);
}
void fp12_add(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp4_add(&ANS->x0,&A->x0,&B->x0);
    fp4_add(&ANS->x1,&A->x1,&B->x1);
    fp4_add(&ANS->x2,&A->x2,&B->x2);
}
void fp12_add_final(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp4_add_final(&ANS->x0,&A->x0,&B->x0);
    fp4_add_final(&ANS->x1,&A->x1,&B->x1);
    fp4_add_final(&ANS->x2,&A->x2,&B->x2);
}
void fp12_add_lazy(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp4_add_lazy(&ANS->x0,&A->x0,&B->x0);
    fp4_add_lazy(&ANS->x1,&A->x1,&B->x1);
    fp4_add_lazy(&ANS->x2,&A->x2,&B->x2);
}

void fp12_add_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp4_add_ui(&ANS->x0,&A->x0,UI);
    fp4_add_ui(&ANS->x1,&A->x1,0);
    fp4_add_ui(&ANS->x2,&A->x2,0);
}

void fp12_add_ui_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp4_add_ui_ui(&ANS->x0,&A->x0,UI);
    fp4_add_ui_ui(&ANS->x1,&A->x1,UI);
    fp4_add_ui_ui(&ANS->x2,&A->x2,UI);
}
void fp12_add_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B){
    fp4_add_mpn(&ANS->x0,&A->x0,B);
    fp4_add_mpn(&ANS->x1,&A->x1,B);
    fp4_add_mpn(&ANS->x2,&A->x2,B);
}

void fp12_sub(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp4_sub(&ANS->x0,&A->x0,&B->x0);
    fp4_sub(&ANS->x1,&A->x1,&B->x1);
    fp4_sub(&ANS->x2,&A->x2,&B->x2);
}
void fp12_sub_final(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp4_sub_final(&ANS->x0,&A->x0,&B->x0);
    fp4_sub_final(&ANS->x1,&A->x1,&B->x1);
    fp4_sub_final(&ANS->x2,&A->x2,&B->x2);
}
void fp12_sub_lazy(fp12_t *ANS,fp12_t *A,fp12_t *B){
    fp4_sub_lazy(&ANS->x0,&A->x0,&B->x0);
    fp4_sub_lazy(&ANS->x1,&A->x1,&B->x1);
    fp4_sub_lazy(&ANS->x2,&A->x2,&B->x2);
}

void fp12_sub_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp4_sub_ui(&ANS->x0,&A->x0,UI);
    fp4_sub_ui(&ANS->x1,&A->x1,0);
    fp4_sub_ui(&ANS->x2,&A->x2,0);
}

void fp12_sub_ui_ui(fp12_t *ANS,fp12_t *A,unsigned long int UI){
    fp4_sub_ui_ui(&ANS->x0,&A->x0,UI);
    fp4_sub_ui_ui(&ANS->x1,&A->x1,UI);
    fp4_sub_ui_ui(&ANS->x2,&A->x2,UI);
}
void fp12_sub_mpn(fp12_t *ANS,fp12_t *A,mp_limb_t *B){
    fp4_sub_mpn(&ANS->x0,&A->x0,B);
    fp4_sub_mpn(&ANS->x1,&A->x1,B);
    fp4_sub_mpn(&ANS->x2,&A->x2,B);
}

void fp12_inv(fp12_t *ANS,fp12_t *A){
	static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4,tmp6_fp4,tmp7_fp4,tmp8_fp4;

    fp4_sqr(&tmp1_fp4,&A->x0);
    fp4_sqr(&tmp2_fp4,&A->x1);
    fp4_sqr(&tmp3_fp4,&A->x2);
    
    fp4_mul(&tmp4_fp4,&A->x1,&A->x2);
    fp4_mul_basis(&tmp4_fp4,&tmp4_fp4);
    fp4_sub(&tmp6_fp4,&tmp1_fp4,&tmp4_fp4);
    
    fp4_mul(&tmp4_fp4,&A->x0,&A->x1);
    fp4_mul_basis(&tmp7_fp4,&tmp3_fp4);
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp4_fp4);
    
    fp4_mul(&tmp4_fp4,&A->x0,&A->x2);
    fp4_sub(&tmp8_fp4,&tmp2_fp4,&tmp4_fp4);
    
    fp4_mul(&tmp1_fp4,&tmp1_fp4,&A->x0);
    fp4_mul(&tmp3_fp4,&tmp3_fp4,&A->x2);
    fp4_mul_basis(&tmp3_fp4,&tmp3_fp4);
    
    fp4_add(&tmp5_fp4,&tmp4_fp4,&tmp4_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp4_fp4);
    fp4_sub(&tmp5_fp4,&tmp2_fp4,&tmp5_fp4);
    fp4_mul(&tmp5_fp4,&tmp5_fp4,&A->x1);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp3_fp4);
    fp4_mul_basis(&tmp5_fp4,&tmp5_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp1_fp4);
    
    fp4_inv(&tmp5_fp4,&tmp5_fp4);
    
    fp4_mul(&ANS->x0,&tmp6_fp4,&tmp5_fp4);
    fp4_mul(&ANS->x1,&tmp7_fp4,&tmp5_fp4);
    fp4_mul(&ANS->x2,&tmp8_fp4,&tmp5_fp4);
}
void fp12_inv_lazy(fp12_t *ANS,fp12_t *A){
	static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4,tmp6_fp4,tmp7_fp4,tmp8_fp4;

    fp4_sqr_lazy(&tmp1_fp4,&A->x0);
    fp4_sqr_lazy(&tmp2_fp4,&A->x1);
    fp4_sqr_lazy(&tmp3_fp4,&A->x2);
    
    fp4_mul_lazy(&tmp4_fp4,&A->x1,&A->x2);
    fp4_mul_basis(&tmp4_fp4,&tmp4_fp4);
    fp4_sub(&tmp6_fp4,&tmp1_fp4,&tmp4_fp4);//tmp6
    
    fp4_mul_lazy(&tmp4_fp4,&A->x0,&A->x1);
    fp4_mul_basis(&tmp7_fp4,&tmp3_fp4);
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp4_fp4);//tmp7
    
    fp4_mul_lazy(&tmp4_fp4,&A->x0,&A->x2);//tmp4
    fp4_sub(&tmp8_fp4,&tmp2_fp4,&tmp4_fp4);//tmp8
    
    fp4_mul_lazy(&tmp1_fp4,&tmp1_fp4,&A->x0);//tmp1
    fp4_mul_lazy(&tmp3_fp4,&tmp3_fp4,&A->x2);
    fp4_mul_basis(&tmp3_fp4,&tmp3_fp4);//tmp3
    
    fp4_add_lazy(&tmp5_fp4,&tmp4_fp4,&tmp4_fp4);
    fp4_add_lazy(&tmp5_fp4,&tmp5_fp4,&tmp4_fp4);
    fp4_sub_lazy(&tmp5_fp4,&tmp2_fp4,&tmp5_fp4);
    fp4_mul_lazy(&tmp5_fp4,&tmp5_fp4,&A->x1);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp3_fp4);
    fp4_mul_basis(&tmp5_fp4,&tmp5_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp1_fp4);//mod
    
    fp4_inv_lazy(&tmp5_fp4,&tmp5_fp4);
    
    fp4_mul_lazy(&ANS->x0,&tmp6_fp4,&tmp5_fp4);
    fp4_mul_lazy(&ANS->x1,&tmp7_fp4,&tmp5_fp4);
    fp4_mul_lazy(&ANS->x2,&tmp8_fp4,&tmp5_fp4);
}
void fp12_inv_lazy_montgomery(fp12_t *ANS,fp12_t *A){
	static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4,tmp4_fp4,tmp5_fp4,tmp6_fp4,tmp7_fp4,tmp8_fp4;

    fp4_sqr_lazy_montgomery(&tmp1_fp4,&A->x0);
    fp4_sqr_lazy_montgomery(&tmp2_fp4,&A->x1);
    fp4_sqr_lazy_montgomery(&tmp3_fp4,&A->x2);
    
    fp4_mul_lazy_montgomery(&tmp4_fp4,&A->x1,&A->x2);
    fp4_mul_basis(&tmp4_fp4,&tmp4_fp4);
    fp4_sub(&tmp6_fp4,&tmp1_fp4,&tmp4_fp4);//tmp6
    
    fp4_mul_lazy_montgomery(&tmp4_fp4,&A->x0,&A->x1);
    fp4_mul_basis(&tmp7_fp4,&tmp3_fp4);
    fp4_sub(&tmp7_fp4,&tmp7_fp4,&tmp4_fp4);//tmp7
    
    fp4_mul_lazy_montgomery(&tmp4_fp4,&A->x0,&A->x2);//tmp4
    fp4_sub(&tmp8_fp4,&tmp2_fp4,&tmp4_fp4);//tmp8
    
    fp4_mul_lazy_montgomery(&tmp1_fp4,&tmp1_fp4,&A->x0);//tmp1
    fp4_mul_lazy_montgomery(&tmp3_fp4,&tmp3_fp4,&A->x2);
    fp4_mul_basis(&tmp3_fp4,&tmp3_fp4);//tmp3
    
    fp4_add_lazy(&tmp5_fp4,&tmp4_fp4,&tmp4_fp4);
    fp4_add_lazy(&tmp5_fp4,&tmp5_fp4,&tmp4_fp4);
    fp4_sub_lazy(&tmp5_fp4,&tmp2_fp4,&tmp5_fp4);
    fp4_mul_lazy_montgomery(&tmp5_fp4,&tmp5_fp4,&A->x1);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp3_fp4);
    fp4_mul_basis(&tmp5_fp4,&tmp5_fp4);
    fp4_add(&tmp5_fp4,&tmp5_fp4,&tmp1_fp4);//mod
    
    fp4_inv_lazy_montgomery(&tmp5_fp4,&tmp5_fp4);
    
    fp4_mul_lazy_montgomery(&ANS->x0,&tmp6_fp4,&tmp5_fp4);
    fp4_mul_lazy_montgomery(&ANS->x1,&tmp7_fp4,&tmp5_fp4);
    fp4_mul_lazy_montgomery(&ANS->x2,&tmp8_fp4,&tmp5_fp4);
}
int  fp12_legendre(fp12_t *A){
    mpz_t exp;
    mpz_init(exp);
    fp12_t tmp;
    fp12_init(&tmp);
    
    mpz_pow_ui(exp,prime_z,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    fp12_pow(&tmp,A,exp);
    
    if(fp12_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

int  fp12_isCNR(fp12_t *A){
    fp12_t tmp;
    fp12_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime_z,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    fp12_pow(&tmp,A,exp);
    
    if(fp12_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}

void fp12_sqrt(fp12_t *ANS,fp12_t *A){
    fp12_t tmp1,tmp2;
    fp12_init(&tmp1);
    fp12_init(&tmp2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    fp4_set(&tmp1.x0,&A->x0);
    fp4_mul_mpn(&tmp1.x1,&A->x1,frobenius_constant[f_p4][1].x0.x0);
    fp4_mul_mpn(&tmp1.x2,&A->x2,frobenius_constant[f_p4][2].x0.x0);
    
    fp4_set(&tmp2.x0,&A->x0);
    fp4_mul_mpn(&tmp2.x1,&A->x1,frobenius_constant[f_p2][1].x0.x0);
    fp4_mul_mpn(&tmp2.x2,&A->x2,frobenius_constant[f_p2][2].x0.x0);
    
    fp12_mul(&tmp1,&tmp1,&tmp2);
    fp12_mul(&tmp1,&tmp1,A);
    fp12_set_ui(&tmp2,0);
    fp4_sqrt(&tmp2.x0,&tmp1.x0);
    fp4_inv(&tmp2.x0,&tmp2.x0);
    fp4_set(&tmp2.x0,&tmp2.x0);
    mpz_pow_ui(exp,prime_z,8);
    mpz_pow_ui(buf,prime_z,4);
    mpz_add(exp,exp,buf);
    mpz_add_ui(exp,exp,2);
    mpz_tdiv_q_ui(exp,exp,2);
    fp12_pow(&tmp1,A,exp);
    fp12_mul(&tmp1,&tmp1,&tmp2);
    fp12_set(ANS,&tmp1);
    
    mpz_clear(exp);
    mpz_clear(buf);
}

void fp12_pow(fp12_t *ANS,fp12_t *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    fp12_t tmp;
    fp12_init(&tmp);
    fp12_set(&tmp,A);
    
    for(i=1;i<length; i++){
        fp12_sqr(&tmp,&tmp);
        if(binary[i]=='1'){
            fp12_mul(&tmp,A,&tmp);
        }
    }
    
    fp12_set(ANS,&tmp);
}

int  fp12_cmp(fp12_t *A,fp12_t *B){
    if(fp4_cmp(&A->x0,&B->x0)==0 && fp4_cmp(&A->x1,&B->x1)==0 && fp4_cmp(&A->x2,&B->x2)==0){
        return 0;   
    }
    return 1;
}

int  fp12_cmp_ui(fp12_t *A,unsigned long int UI){
    if(fp4_cmp_ui(&A->x0,UI)==0 && fp4_cmp_ui(&A->x1,UI)==0 && fp4_cmp_ui(&A->x2,UI)==0){
        return 0;
    }
    return 1;
}

int  fp12_cmp_mpn(fp12_t *A,mp_limb_t *B){
    if(fp4_cmp_mpn(&A->x0,B)==0 && fp4_cmp_mpn(&A->x1,B)==0 && fp4_cmp_mpn(&A->x2,B)==0){
        return 0;
    }
    return 1;
}

int  fp12_cmp_zero(fp12_t *A){
    if(fp4_cmp_zero(&A->x0)==0 && fp4_cmp_zero(&A->x1)==0 && fp4_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}

int  fp12_cmp_one(fp12_t *A){
    if(fp4_cmp_one(&A->x0)==0 && fp4_cmp_zero(&A->x1)==0 && fp4_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}
