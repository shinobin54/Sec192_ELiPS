#include <ELiPS/efp4.h>
//efp4
void efp4_init(efp4_t *P){
    fp4_init(&P->x);
    fp4_init(&P->y);
    P->infinity=0;
}
void efp4_printf(char *str,efp4_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp4_printf("",&P->x);
        printf(",");
        fp4_printf("",&P->y);
        printf(")");
    }else{
        printf("0");
    }
}
void efp4_println(char *str,efp4_t *P){
    printf("%s",str);
    if(P->infinity==0){
        printf("(");
        fp4_printf("",&P->x);
        printf(",");
        fp4_printf("",&P->y);
        printf(")\n");
    }else{
        printf("0\n");
    }
}
void efp4_set(efp4_t *ANS,efp4_t *A){
    fp4_set(&ANS->x,&A->x);
    fp4_set(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

void efp4_set_ui(efp4_t *ANS,unsigned long int UI1,unsigned long int UI2){
    fp4_set_ui(&ANS->x,UI1);
    fp4_set_ui(&ANS->y,UI2);
    ANS->infinity=0;
}

void efp4_set_mpn(efp4_t *ANS,mp_limb_t *A){
    fp4_set_mpn(&ANS->x,A);
    fp4_set_mpn(&ANS->y,A);
    ANS->infinity=0;
}

void efp4_set_neg(efp4_t *ANS,efp4_t *A){
    fp4_set(&ANS->x,&A->x);
    fp4_set_neg(&ANS->y,&A->y);
    ANS->infinity=A->infinity;
}

int  efp4_cmp(efp4_t *A,efp4_t *B){
    if(fp4_cmp(&A->x,&B->x)==0 && fp4_cmp(&A->y,&B->y)==0){
        return 0;   
    }else if(A->infinity==1&&B->infinity==1){
	return 0;
    }else{
    return 1;
    }
}

void efp4_rational_point(efp4_t *P){
    fp4_t tmp1,tmp2;
    fp4_init(&tmp1);
    fp4_init(&tmp2);
    while(1){
        fp4_set_random(&P->x,state);
        fp4_sqr(&tmp1,&P->x);
        fp4_mul(&tmp2,&tmp1,&P->x);
        fp_add_mpn(&tmp2.x0.x0,&tmp2.x0.x0,curve_b);
        if(fp4_legendre(&tmp2)==1){
            fp4_sqrt(&P->y,&tmp2);
            break;
        }
    }
}
void efp4_ecd(efp4_t *ANS,efp4_t *P){
    static efp4_t tmp1_efp4;
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4;
    if(fp4_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    efp4_set(&tmp1_efp4,P);
    
    fp4_add(&tmp1_fp4,&tmp1_efp4.y,&tmp1_efp4.y);
    
    fp4_inv(&tmp1_fp4,&tmp1_fp4);
    fp4_sqr(&tmp2_fp4,&tmp1_efp4.x);
    fp4_add(&tmp3_fp4,&tmp2_fp4,&tmp2_fp4);
    fp4_add(&tmp2_fp4,&tmp2_fp4,&tmp3_fp4);
    fp4_mul(&tmp3_fp4,&tmp1_fp4,&tmp2_fp4);
    
    fp4_sqr(&tmp1_fp4,&tmp3_fp4);
    fp4_add(&tmp2_fp4,&tmp1_efp4.x,&tmp1_efp4.x);
    fp4_sub(&ANS->x,&tmp1_fp4,&tmp2_fp4);
    
    fp4_sub(&tmp1_fp4,&tmp1_efp4.x,&ANS->x);
    fp4_mul(&tmp2_fp4,&tmp3_fp4,&tmp1_fp4);
    fp4_sub(&ANS->y,&tmp2_fp4,&tmp1_efp4.y);
}

void efp4_ecd_lazy(efp4_t *ANS,efp4_t *P){
    static efp4_t tmp1_efp4;
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4;
    if(fp4_cmp_zero(&P->y)==0){
        ANS->infinity=1;
        return;
    }
    
    efp4_set(&tmp1_efp4,P);
    
    fp4_add(&tmp1_fp4,&tmp1_efp4.y,&tmp1_efp4.y);
    
    fp4_inv(&tmp1_fp4,&tmp1_fp4);
    fp4_sqr_lazy(&tmp2_fp4,&tmp1_efp4.x);
    fp4_add_lazy(&tmp3_fp4,&tmp2_fp4,&tmp2_fp4);
    fp4_add_lazy(&tmp2_fp4,&tmp2_fp4,&tmp3_fp4);
    fp4_mul_lazy(&tmp3_fp4,&tmp1_fp4,&tmp2_fp4);
    
    fp4_sqr_lazy(&tmp1_fp4,&tmp3_fp4);
    fp4_add(&tmp2_fp4,&tmp1_efp4.x,&tmp1_efp4.x);
    fp4_sub(&ANS->x,&tmp1_fp4,&tmp2_fp4);
    
    fp4_sub_lazy(&tmp1_fp4,&tmp1_efp4.x,&ANS->x);
    fp4_mul_lazy(&tmp2_fp4,&tmp3_fp4,&tmp1_fp4);
    fp4_sub(&ANS->y,&tmp2_fp4,&tmp1_efp4.y);
}

void efp4_eca(efp4_t *ANS,efp4_t *P1,efp4_t *P2){
    static efp4_t tmp1_efp4,tmp2_efp4;
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4;
    if(P1->infinity==1){
        efp4_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp4_set(ANS,P1);
        return;
    }else if(fp4_cmp(&P1->x,&P2->x)==0){
        if(fp4_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp4_ecd(ANS,P1);
            return;
        }
    }
    
    efp4_set(&tmp1_efp4,P1);
    efp4_set(&tmp2_efp4,P2);
    
    fp4_sub(&tmp1_fp4,&tmp2_efp4.x,&tmp1_efp4.x);
    fp4_inv(&tmp1_fp4,&tmp1_fp4);
    fp4_sub(&tmp2_fp4,&tmp2_efp4.y,&tmp1_efp4.y);
    fp4_mul(&tmp3_fp4,&tmp1_fp4,&tmp2_fp4);
    fp4_sqr(&tmp1_fp4,&tmp3_fp4);
    fp4_sub(&tmp2_fp4,&tmp1_fp4,&tmp1_efp4.x);
    fp4_sub(&ANS->x,&tmp2_fp4,&tmp2_efp4.x);
    fp4_sub(&tmp1_fp4,&tmp1_efp4.x,&ANS->x);
    fp4_mul(&tmp2_fp4,&tmp3_fp4,&tmp1_fp4);
    fp4_sub(&ANS->y,&tmp2_fp4,&tmp1_efp4.y);
}
void efp4_eca_lazy(efp4_t *ANS,efp4_t *P1,efp4_t *P2){
    static efp4_t tmp1_efp4,tmp2_efp4;
    static fp4_t tmp1_fp4,tmp2_fp4,tmp3_fp4;
    if(P1->infinity==1){
        efp4_set(ANS,P2);
        return;
    }else if(P2->infinity==1){
        efp4_set(ANS,P1);
        return;
    }else if(fp4_cmp(&P1->x,&P2->x)==0){
        if(fp4_cmp(&P1->y,&P2->y)!=0){
            ANS->infinity=1;
            return;
        }else{
            efp4_ecd_lazy(ANS,P1);
            return;
        }
    }
    
    efp4_set(&tmp1_efp4,P1);
    efp4_set(&tmp2_efp4,P2);
    
    fp4_sub(&tmp1_fp4,&tmp2_efp4.x,&tmp1_efp4.x);
    fp4_inv(&tmp1_fp4,&tmp1_fp4);
    fp4_sub_lazy(&tmp2_fp4,&tmp2_efp4.y,&tmp1_efp4.y);
    fp4_mul_lazy(&tmp3_fp4,&tmp1_fp4,&tmp2_fp4);
    fp4_sqr_lazy(&tmp1_fp4,&tmp3_fp4);
    fp4_sub_lazy(&tmp2_fp4,&tmp1_fp4,&tmp1_efp4.x);
    fp4_sub(&ANS->x,&tmp2_fp4,&tmp2_efp4.x);
    fp4_sub_lazy(&tmp1_fp4,&tmp1_efp4.x,&ANS->x);
    fp4_mul_lazy(&tmp2_fp4,&tmp3_fp4,&tmp1_fp4);
    fp4_sub(&ANS->y,&tmp2_fp4,&tmp1_efp4.y);
}
void efp4_scm(efp4_t *ANS,efp4_t *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        efp4_set(ANS,P);
        return;
    }
    
    efp4_t Tmp_P,Next_P;
    efp4_init(&Tmp_P);
    efp4_set(&Tmp_P,P);
    efp4_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    efp4_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        efp4_ecd(&Next_P,&Next_P);
        if(binary[i]=='1'){
            efp4_eca(&Next_P,&Next_P,&Tmp_P);
        }
    }
    efp4_set(ANS,&Next_P);
}
void efp4_scm_lazy(efp4_t *ANS,efp4_t *P,mpz_t scalar){
    if(mpz_cmp_ui(scalar,0)==0){
        ANS->infinity=1;
        return;
    }else if(mpz_cmp_ui(scalar,1)==0){
        efp4_set(ANS,P);
        return;
    }
    
    efp4_t Tmp_P,Next_P;
    efp4_init(&Tmp_P);
    efp4_set(&Tmp_P,P);
    efp4_init(&Next_P);
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    
    efp4_set(&Next_P,&Tmp_P);
    for(i=1;i<length; i++){
        efp4_ecd_lazy(&Next_P,&Next_P);
        if(binary[i]=='1'){
            efp4_eca_lazy(&Next_P,&Next_P,&Tmp_P);
        }
    }
    efp4_set(ANS,&Next_P);
}