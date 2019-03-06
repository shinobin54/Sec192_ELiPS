#include <ELiPS/Fp.h>

void Fp_init(Fp *A){
    mpn_zero(A->x0,FPLIMB);
}

void Fp_printf(Fp *A,char *str){
    gmp_printf("%s%Nu",str,A->x0,FPLIMB);
}

void Fp_set(Fp *ANS,Fp *A){
    mpn_copyd(ANS->x0,A->x0,FPLIMB);
}

void Fp_set_ui(Fp *ANS,unsigned long int UI){
    mpn_set_ui(ANS->x0,FPLIMB,UI);
}
void Fp_set_mpn(Fp *ANS,mp_limb_t *A){
    mpn_copyd(ANS->x0,A,FPLIMB);
}

void Fp_set_neg(Fp *ANS,Fp *A){
    mpn_sub_n(ANS->x0,prime,A->x0,FPLIMB);
    mpn_mod(ANS,ANS->x0,FPLIMB);
}
void Fp_lshift(Fp *ANS,Fp *A,unsigned long int UI){
    static mp_limb_t buf[FPLIMB];
    mpn_lshift(buf,A->x0,FPLIMB,UI);
    mpn_mod(ANS,buf,FPLIMB);
}
void Fp_set_random(Fp *ANS,gmp_randstate_t state){
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_urandomm(tmp,state,prime_z);
    mpn_set_mpz(ANS->x0,tmp);
    /*
    mpn_random(buf,FPLIMB);
    mpn_mod(ANS,buf,FPLIMB);
    */
    mpz_clear(tmp);
}

void Fp_MR(mp_limb_t *ANS,mp_limb_t *T,mp_size_t T_size){
    static mp_limb_t s[FPLIMB];
    static mp_limb_t buf[FPLIMB],bufL[FPLIMB2];

    mpn_and_n(buf,T,R1,FPLIMB);
    mpn_mul_n(bufL,buf,Ni,FPLIMB);
    mpn_and_n(s,bufL,R1,FPLIMB);

    mpn_mul_n(bufL,s,prime,FPLIMB);
    mpn_add(bufL,bufL,FPLIMB2,T,T_size);
    mpn_rshift_ext(bufL,bufL,FPLIMB2,m);
    mpn_copyd(buf,bufL,FPLIMB);

    if(mpn_cmp(buf,prime,FPLIMB)>0)mpn_sub_n(ANS,buf,prime,FPLIMB);
    else mpn_copyd(ANS,buf,FPLIMB);
}

void Lazy_add(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size,mp_limb_t *B,mp_size_t B_size){
    if(mpn_zero_p(A,A_size)==1){
	mpn_copyd(ANS,B,B_size);
    }else if(mpn_zero_p(B,B_size)==1){
	mpn_copyd(ANS,A,A_size);
    }else{
	mpn_add_n(ANS,A,B,ANS_size);
    }
}
void Lazy_sub(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size,mp_limb_t *B,mp_size_t B_size){
    static mp_limb_t buf[FPLIMB],bufL[FPLIMB2];
    //Only:A_size>=B_size
    if(mpn_zero_p(B,B_size)==1){
	//mpn_zero(ANS,ANS_size);
	mpn_copyd(ANS,A,A_size);
    }else if(mpn_zero_p(A,A_size)==1){
	    mpn_zero(ANS,ANS_size);
	    Lazy_mod(buf,B,B_size);
    	    mpn_sub_n(ANS,prime,buf,FPLIMB);
    }else{
	if(A_size>B_size){
	    if(mpn_chk_limb(A,FPLIMB,FPLIMB2)==0){
	        mpn_sub(ANS,A,A_size,B,B_size);
            }else{
		if(mpn_cmp(A,B,FPLIMB)<0){
		    mpn_sub_n(buf,B,A,FPLIMB);
		    Lazy_mod(buf,buf,FPLIMB);
		    mpn_zero(ANS,ANS_size);
    		    mpn_sub_n(ANS,prime,buf,FPLIMB);
                }else{
		    mpn_zero(ANS,ANS_size);
		    mpn_sub_n(ANS,A,B,FPLIMB);
	        }
	    }
	}else{
	    if(mpn_cmp(A,B,ANS_size)<0){
    		mpn_sub_n(ANS,B,A,ANS_size);
		Lazy_mod(buf,ANS,ANS_size);
		mpn_zero(ANS,ANS_size);
    		mpn_sub_n(ANS,prime,buf,FPLIMB);
            }else{
		mpn_sub_n(ANS,A,B,ANS_size);
	    }
	}
    }
}

void Lazy_sub_mod(Fp *ANS,mp_limb_t *A,mp_limb_t *B){
    //Only:A_size=B_size=FPLIMB2
    static mp_limb_t buf[FPLIMB2];
    if(mpn_zero_p(B,FPLIMB2)==1){
	mpn_mod(ANS,A,FPLIMB2);
    }else{
	if(mpn_cmp(A,B,FPLIMB2)<0){
    	    mpn_sub_n(buf,B,A,FPLIMB2);
	    mpn_mod(ANS,buf,FPLIMB2);
    	    mpn_sub_n(ANS->x0,prime,ANS->x0,FPLIMB);
        }else{
	    mpn_sub_n(buf,A,B,FPLIMB2);
	    mpn_mod(ANS,buf,FPLIMB2);
	}
    }
}
void Lazy_mul(mp_limb_t *ANS,mp_limb_t *A,mp_limb_t *B){
    if(mpn_zero_p(A,FPLIMB)==1||mpn_zero_p(B,FPLIMB)==1){
	mpn_set_ui(ANS,FPLIMB2,0);
    }else if(mpn_cmp_ui(A,FPLIMB,1)==0){
	mpn_zero(ANS,FPLIMB2);
	mpn_copyd(ANS,B,FPLIMB);
    }else if(mpn_cmp_ui(B,FPLIMB,1)==0){
	mpn_zero(ANS,FPLIMB2);
	mpn_copyd(ANS,A,FPLIMB);
    }else{
        mpn_mul_n(ANS,A,B,FPLIMB);
    }

}
void Lazy_sqr(mp_limb_t *ANS,mp_limb_t *A){
    if(mpn_zero_p(A,FPLIMB)==1){
	mpn_set_ui(ANS,FPLIMB2,1);
    }else if(mpn_cmp_ui(A,FPLIMB,0)==0){
	mpn_set_ui(ANS,FPLIMB2,0);
    }else{
        mpn_sqr(ANS,A,FPLIMB);
    }

}

void Lazy_mod(mp_limb_t *ans,mp_limb_t *a,mp_size_t size_a){
    mp_limb_t dumy[size_a];
    mpn_tdiv_qr(dumy,ans,0,a,size_a,prime,FPLIMB);
}

void Fp_mul(Fp *ANS,Fp *A,Fp *B){
    static mp_limb_t tmp_mul[FPLIMB2];
    if(mpn_zero_p(A->x0,FPLIMB)==1||mpn_zero_p(B->x0,FPLIMB)==1){
	mpn_set_ui(ANS->x0,FPLIMB,0);
    }else if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
	Fp_set(ANS,B);
    }else if(mpn_cmp_ui(B->x0,FPLIMB,1)==0){
        Fp_set(ANS,A);
    }else{
        mpn_mul_n(tmp_mul,A->x0,B->x0,FPLIMB);
        mpn_mod(ANS,tmp_mul,FPLIMB2);
    }
}
void Fp_mul_montgomery(mp_limb_t *ANS,mp_size_t ANS_size,mp_limb_t *A,mp_size_t A_size){
    static mp_limb_t tmp_mul[FPLIMB2];
    if(mpn_zero_p(A,A_size)==1){
	mpn_set_ui(ANS,ANS_size,0);
    }else if(mpn_cmp_ui(A,A_size,1)==0){
	mpn_zero(ANS,ANS_size);
	mpn_copyd(ANS,RR,FPLIMB);
    }else{
        mpn_mul_n(ANS,A,RR,FPLIMB);
    }
}
void Fp_mul_final(Fp *ANS,Fp *A,Fp *B){
	static mp_limb_t tmp_mul[FPLIMB2];
    if(mpn_zero_p(A->x0,FPLIMB)==1||mpn_zero_p(B->x0,FPLIMB)==1){
	mpn_set_ui(ANS->x0,FPLIMB,0);
    }else if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
	mpn_mod(ANS,B->x0,FPLIMB);
    }else if(mpn_cmp_ui(B->x0,FPLIMB,1)==0){
	mpn_mod(ANS,A->x0,FPLIMB);
    }else{
        mpn_mul_n(tmp_mul,A->x0,B->x0,FPLIMB);
        mpn_mod(ANS,tmp_mul,FPLIMB2);
    }
}
void Fp_mul_ui(Fp *ANS,Fp *A,unsigned long int UI){
    static mp_limb_t tmp_mul[FPLIMB2];
    mpn_mul_ui(tmp_mul,A->x0,FPLIMB,UI);
    mpn_mod(ANS,tmp_mul,FPLIMB2);
}
void Fp_mul_mpn(Fp *ANS,Fp *A,mp_limb_t *B){
    static mp_limb_t tmp_mul[FPLIMB2];
    if(mpn_zero_p(A->x0,FPLIMB)==1||mpn_zero_p(B,FPLIMB)==1){
	mpn_set_ui(ANS->x0,FPLIMB,0);
    }else if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
	mpn_copyd(ANS->x0,B,FPLIMB);
    }else if(mpn_cmp_ui(B,FPLIMB,1)==0){
        Fp_set(ANS,A);
    }else{
        mpn_mul_n(tmp_mul,A->x0,B,FPLIMB);
        mpn_mod(ANS,tmp_mul,FPLIMB2);
    }
}
void Fp_add(Fp *ANS,Fp *A,Fp *B){
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	Fp_set(ANS,B);
    }else if(mpn_zero_p(B->x0,FPLIMB)==1){
	Fp_set(ANS,A);
    }else{
	mpn_add_n(buf,A->x0,B->x0,FPLIMB);
	if(mpn_cmp(buf,prime,FPLIMB)>=0)mpn_sub_n(ANS->x0,buf,prime,FPLIMB);
	else mpn_copyd(ANS->x0,buf,FPLIMB);
    }
}

void Fp_add_final(Fp *ANS,Fp *A,Fp *B){
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	mpn_mod(ANS,B->x0,FPLIMB);
    }else if(mpn_zero_p(B->x0,FPLIMB)==1){
	mpn_mod(ANS,A->x0,FPLIMB);
    }else{
	mpn_add_n(buf,A->x0,B->x0,FPLIMB);
	mpn_mod(ANS,buf,FPLIMB);
    }

}
void Fp_add_ui(Fp *ANS,Fp *A,unsigned long int UI){
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	Fp_set_ui(ANS,UI);
    }else if(UI==0){
	Fp_set(ANS,A);
    }else{
	mpn_add_ui(buf,A->x0,FPLIMB,UI);
	if(mpn_cmp(buf,prime,FPLIMB)>0)mpn_sub_n(ANS->x0,buf,prime,FPLIMB);
	else mpn_copyd(ANS->x0,buf,FPLIMB);
    }
}
void Fp_add_mpn(Fp *ANS,Fp *A,mp_limb_t *B){
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(A->x0,FPLIMB)==1){
	mpn_copyd(ANS->x0,B,FPLIMB);
    }else if(mpn_zero_p(B,FPLIMB)==1){
	Fp_set(ANS,A);
    }else{
	mpn_add_n(buf,A->x0,B,FPLIMB);
	if(mpn_cmp(buf,prime,FPLIMB)>0)mpn_sub_n(ANS->x0,buf,prime,FPLIMB);
	else mpn_copyd(ANS->x0,buf,FPLIMB);
    }
}
void Fp_sub(Fp *ANS,Fp *A,Fp *B){
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(B->x0,FPLIMB)==1){
	Fp_set(ANS,A);
    }else{
        if(mpn_cmp(A->x0,B->x0,FPLIMB)<0){
    	    mpn_sub_n(buf,B->x0,A->x0,FPLIMB);
    	    mpn_sub_n(ANS->x0,prime,buf,FPLIMB);
        }else{
	    mpn_sub_n(ANS->x0,A->x0,B->x0,FPLIMB);
	}
    }
}

void Fp_sub_final(Fp *ANS,Fp *A,Fp *B){
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(B->x0,FPLIMB)==1){
	mpn_mod(ANS,A->x0,FPLIMB);
    }else{
        if(mpn_cmp(A->x0,B->x0,FPLIMB)<0){
    	    mpn_sub_n(buf,B->x0,A->x0,FPLIMB);
	    mpn_mod(ANS,buf,FPLIMB);
    	    mpn_sub_n(ANS->x0,prime,ANS->x0,FPLIMB);
        }else{
	    mpn_sub_n(buf,A->x0,B->x0,FPLIMB);
	    mpn_mod(ANS,buf,FPLIMB);
	}
    }
}
void Fp_sub_ui(Fp *ANS,Fp *A,unsigned long int UI){
	static mp_limb_t buf[FPLIMB];
    if(UI==0){
	Fp_set(ANS,A);
    }else{
	mpn_sub_ui(ANS->x0,A->x0,FPLIMB,UI);
    }
}
void Fp_sub_mpn(Fp *ANS,Fp *A,mp_limb_t *B){
    static mp_limb_t buf[FPLIMB];
    if(mpn_zero_p(B,FPLIMB)==1){
	Fp_set(ANS,A);
    }else{
        if(mpn_cmp(A->x0,B,FPLIMB)<0){
    	    mpn_sub_n(buf,B,A->x0,FPLIMB);
    	    mpn_sub_n(buf,prime,buf,FPLIMB);
	    mpn_mod(ANS,buf,FPLIMB);
        }else{
	    mpn_sub_n(buf,A->x0,B,FPLIMB);
	    mpn_mod(ANS,buf,FPLIMB);
	}
    }
}
void Fp_inv(Fp *ANS,Fp *A){
    static mp_limb_t prime_tmp[FPLIMB],gp[FPLIMB],sp[FPLIMB],tmp[FPLIMB],buf[FPLIMB];
    static mp_size_t buf_size;

    //mpn_init(gp,FPLIMB);
    mpn_init(sp,FPLIMB);
    //mpn_init(tmp,FPLIMB);
    //mpn_init(prime_tmp,FPLIMB);

    mpn_copyd(prime_tmp,prime,FPLIMB);
	
    mpn_add_n(buf,A->x0,prime,FPLIMB);
    mpn_gcdext(gp,sp,&buf_size,buf,FPLIMB,prime_tmp,FPLIMB);
	
    if( buf_size < 0 ){
        mpn_sub_n(tmp, prime,sp,FPLIMB);
    }else{	
	mpn_copyd(tmp,sp,FPLIMB);
    }
    
    mpn_mod(ANS,tmp,FPLIMB);
}
int  Fp_legendre(Fp *A){

    int i;
    mpz_t tmp1,tmp2;
    Fp tmp1_Fp;
    mpz_init(tmp1);
    mpz_init(tmp2);
    Fp_init(&tmp1_Fp);
    
    mpz_sub_ui(tmp1,prime_z,1);
    mpz_tdiv_q_ui(tmp2,tmp1,2);
    Fp_pow(&tmp1_Fp,A,tmp2);
	
    if(mpn_cmp_char(tmp1_Fp.x0,"1")==0)		i=1;
    else if(mpn_cmp_char(tmp1_Fp.x0,"0")==0)	i=0;
    else					i=-1;
    
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    
    return i;
}

int  Fp_isCNR(Fp *A){
    Fp tmp;
    Fp_init(&tmp);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_sub_ui(exp,prime_z,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp_pow(&tmp,A,exp);
    
    if(Fp_cmp_one(&tmp)==0){
        mpz_clear(exp);
        return 1;
    }else{
        mpz_clear(exp);
        return -1;
    }
}
void Fp_sqrt(Fp *ANS,Fp *A){
    Fp x,y,t,k,n,tmp;
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&t);
    Fp_init(&k);
    Fp_init(&n);
    Fp_init(&tmp);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp_set_random(&n,state);
    
    while(Fp_legendre(&n)!=-1){
        Fp_set_random(&n,state);
    }
    mpz_sub_ui(q,prime_z,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp_pow(&x,A,exp);
    Fp_mul(&tmp,&x,&x);
    Fp_mul(&k,&tmp,A);
    Fp_mul(&x,&x,A);
    while(mpn_cmp_ui(k.x0,FPLIMB,1)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp_pow(&tmp,&k,exp);
        while(mpn_cmp_ui(tmp.x0,FPLIMB,1)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp_pow(&tmp,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp_pow(&t,&y,result);
        Fp_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp_mul(&x,&x,&t);
        Fp_mul(&k,&k,&y);
    }
    Fp_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
}
void Fp_pow(Fp *ANS,Fp *A,mpz_t scalar){
    int i,length;
    length=(int)mpz_sizeinbase(scalar,2);
    char binary[length+1];
    mpz_get_str(binary,2,scalar);
    Fp tmp;
    Fp_init(&tmp);
    
    Fp_set(&tmp,A);
    
    for(i=1;i<length; i++){
        Fp_mul(&tmp,&tmp,&tmp);
        if(binary[i]=='1'){
            Fp_mul(&tmp,A,&tmp);
        }
    }
    Fp_set(ANS,&tmp);
    
}
void Fp_pow_mpn(Fp *ans,Fp *a,mp_limb_t *r,mp_size_t n){
    Fp Temp,tmp;
    mp_limb_t bit[FPLIMB],bit_copy[FPLIMB];
    size_t bit_size;
    mp_size_t size;
    int i,cnt;
    size=FPLIMB;
	
    Fp_init(&Temp);
    Fp_init(&tmp);
    mpn_init(bit,size);
    mpn_init(bit_copy,size);
	
    if(mpn_zero_p(r,size)==1){
    	Fp_set_ui(ans,1);
    	return;
}
    bit_size=mpn_sizeinbase(r,n,2);
    bit_size--;
	
    mpn_set_ui(bit,1,size);
    mpn_lshift_ext(bit,bit,size,bit_size);
	
    //SCM
    mpn_copyd(Temp.x0,a->x0,size);
    while(bit_size>0){
	mpn_copyd(tmp.x0,Temp.x0,size);
	Fp_mul(&Temp,&tmp,&tmp);
	//bit
	mpn_rshift(bit,bit,size,1);
	mpn_and_n(bit_copy,r,bit,size);
	//bit=1 -> 
	if(mpn_zero_p(bit_copy,size)==0){
		mpn_copyd(tmp.x0,Temp.x0,size);
		Fp_mul(&Temp,&tmp,a);
	}
	mpn_init(bit_copy,size);
	bit_size--;
	
    }
    mpn_copyd(ans->x0,Temp.x0,size);
}

int  Fp_cmp(Fp *A,Fp *B){
    if(mpn_cmp(A->x0,B->x0,FPLIMB)==0){
        return 0;   
    }
    return 1;
}

int  Fp_cmp_ui(Fp *A,unsigned long int UI){
    if(mpn_cmp_ui(A->x0,FPLIMB,UI)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_mpn(Fp *A,mp_limb_t *B){
    if(mpn_cmp(A->x0,B,FPLIMB)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_zero(Fp *A){
    if(mpn_cmp_ui(A->x0,FPLIMB,0)==0){
        return 0;
    }
    return 1;
}

int  Fp_cmp_one(Fp *A){
    if(mpn_cmp_ui(A->x0,FPLIMB,1)==0){
        return 0;
    }
    return 1;
}
/*----------------------------------------------------------------------------*/