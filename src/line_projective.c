#include <ELiPS/line_projective.h>
/*mine*/
void FpJ12_init(FpJ12 *A){
    Fp12_init(&A->x);
    Fp12_set_ui(&A->z,1);
}
void FpJ12_set_ui(FpJ12 *A,unsigned int UI){
  Fp12_set_ui(&A->x,UI);
  Fp12_set_ui(&A->z,1);
}
void FpJ12_to_Fp12(Fp12 *ANS,FpJ12 *A){
  static Fp12 tmp1_Fp12;
  Fp12_inv(&tmp1_Fp12,&A->z);
  Fp12_mul(ANS,&A->x,&tmp1_Fp12);
}
void FpJ12_mul(FpJ12 *ANS,FpJ12 *P,FpJ12 *Q){
    Fp12_mul(&ANS->x,&P->x,&Q->x);
    Fp12_mul(&ANS->z,&P->z,&Q->z);
}
void FpJ12_sqr(FpJ12 *ANS,FpJ12 *A){
    Fp12_sqr(&ANS->x,&A->x);
    Fp12_sqr(&ANS->z,&A->z);
}

void ff_ltq_projective_projective(FpJ12 *f,EFpJ12 *T,EFp12 *Q,EFpJ12 *P){
	static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
	//XpZt-Xt
	Fp12_mul(&tmp1_Fp12,&P->x,&T->z);
	Fp12_sub(&f->z,&tmp1_Fp12,&T->x);

	//(Yq-Yp)(XpZt-Xt)
	Fp12_sub(&tmp1_Fp12,&Q->y,&P->y);
	Fp12_mul(&tmp1_Fp12,&tmp1_Fp12,&f->z);

	//YpZt-Yt
	Fp12_mul(&tmp2_Fp12,&P->y,&T->z);
	Fp12_sub(&tmp2_Fp12,&tmp2_Fp12,&T->y);
	//(Xq-Xp)(YpZt-Yt)
	Fp12_sub(&tmp3_Fp12,&Q->x,&P->x);
	Fp12_mul(&tmp2_Fp12,&tmp3_Fp12,&tmp2_Fp12);

	Fp12_sub(&f->x,&tmp1_Fp12,&tmp2_Fp12);
}
void ff_ltt_projective_projective(FpJ12 *f,EFpJ12 *T,EFp12 *Q){
	static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
	//tmp1=2YtZt,f->z=2YtZt^2
	Fp12_mul(&tmp1_Fp12,&T->y,&T->z);
	Fp12_mul_ui(&tmp1_Fp12,&tmp1_Fp12,2);
	Fp12_mul(&f->z,&tmp1_Fp12,&T->z);

	//(YqZt-Yt)*2YtZt
	Fp12_mul(&tmp2_Fp12,&Q->y,&T->z);
	Fp12_sub(&tmp2_Fp12,&tmp2_Fp12,&T->y);
	Fp12_mul(&tmp1_Fp12,&tmp2_Fp12,&tmp1_Fp12);

	//XqZt-Xt
	Fp12_mul(&tmp2_Fp12,&Q->x,&T->z);
	Fp12_sub(&tmp2_Fp12,&tmp2_Fp12,&T->x);
	//(XqZt-Xt)*3Xt^2
	Fp12_sqr(&tmp3_Fp12,&T->x);
	Fp12_mul_ui(&tmp3_Fp12,&tmp3_Fp12,3);
	Fp12_mul(&tmp2_Fp12,&tmp2_Fp12,&tmp3_Fp12);

	Fp12_sub(&f->x,&tmp1_Fp12,&tmp2_Fp12);
}
void ff_ltq_projective_plain(FpJ12 *f,EFp12 *T,EFp12 *Q,EFp12 *P){
	static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
	//Xp-Xt
	Fp12_sub(&f->z,&P->x,&T->x);

	//(Yq-Yp)(Xp-Xt)
	Fp12_sub(&tmp1_Fp12,&Q->y,&P->y);
	Fp12_mul(&tmp1_Fp12,&tmp1_Fp12,&f->z);

	//Yp-Yt
	Fp12_sub(&tmp2_Fp12,&P->y,&T->y);
	//(Xq-Xp)(Yp-Yt)
	Fp12_sub(&tmp3_Fp12,&Q->x,&P->x);
	Fp12_mul(&tmp2_Fp12,&tmp3_Fp12,&tmp2_Fp12);

	Fp12_sub(&f->x,&tmp1_Fp12,&tmp2_Fp12);

  //getchar();
  //Fp12_println("Lx=",&f->x);
  //Fp12_println("Lz=",&f->z);
}
void ff_ltt_projective_plain(FpJ12 *f,EFp12 *T,EFp12 *Q){
	static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
	//tmp1=2YtZt,f->z=2Yt
	Fp12_mul_ui(&f->z,&T->y,2);

	//(Yq-Yt)*2Yt
	Fp12_sub(&tmp1_Fp12,&Q->y,&T->y);
	Fp12_mul(&tmp1_Fp12,&tmp1_Fp12,&f->z);

	//Xq-Xt
	Fp12_sub(&tmp2_Fp12,&Q->x,&T->x);
	//(Xq-Xt)*3Xt^2
	Fp12_sqr(&tmp3_Fp12,&T->x);
	Fp12_mul_ui(&tmp3_Fp12,&tmp3_Fp12,3);
	Fp12_mul(&tmp2_Fp12,&tmp2_Fp12,&tmp3_Fp12);

	Fp12_sub(&f->x,&tmp1_Fp12,&tmp2_Fp12);

  //getchar();
  //Fp12_println("Lx=",&f->x);
  //Fp12_println("Lz=",&f->z);

}
/*
void ff_ltq_projective(FpJ2 *f,EFpJ2 *T,EFp2 *Q,EFpJ2 *P){
	static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
	//XpZt-YtZt
	Fp2_mul(&tmp1_Fp2,&P->x,&T->z);
	Fp2_sub(&f->z,&tmp1_Fp2,&T->y);

	//(Yq-Yp)(XpZt-Yt)
	Fp2_sub(&tmp1_Fp2,&Q->y,&P->y);
	Fp2_mul(&tmp1_Fp2,&tmp1_Fp2,&f->z);

	//YpZt-Yt
	Fp2_mul(&tmp2_Fp2,&P->y,&T->z);
	Fp2_sub(&tmp2_Fp2,&tmp2_Fp2,&T->y);
	//(Xq-Xp)(YpZt-Yt)
	Fp2_sub(&tmp3_Fp2,&Q->x,&P->x);
	Fp2_mul(&tmp2_Fp2,&tmp3_Fp2,&tmp2_Fp2);

	Fp2_sub(&f->x,&tmp1_Fp2,&tmp2_Fp2);
}
void ff_ltt_projective(FpJ2 *f,EFpJ2 *T,EFp2 *Q){
	static Fp2 tmp1_Fp2,tmp2_Fp2,tmp3_Fp2;
	//2YtZt^2
	Fp2_mul(&tmp1_Fp2,&T->y,&T->z);
	Fp2_mul_ui(&tmp1_Fp2,&tmp1_Fp2,2);
	Fp2_mul(&f->z,&tmp1_Fp2,&T->z);

	//(YqZt-Yt)*2YtZt
	Fp2_mul(&tmp2_Fp2,&Q->y,&T->z);
	Fp2_sub(&tmp2_Fp2,&tmp2_Fp2,&T->y);
	Fp2_mul(&tmp1_Fp2,&tmp2_Fp2,&tmp1_Fp2);

	//XqZt-Xt
	Fp2_mul(&tmp2_Fp2,&Q->x,&T->z);
	Fp2_sub(&tmp2_Fp2,&tmp2_Fp2,&T->x);
	//(XqZt-Xt)*3Xt^2
	Fp2_sqr(&tmp3_Fp2,&T->x);
	Fp2_mul_ui(&tmp3_Fp2,&tmp3_Fp2,3);
	Fp2_mul(&tmp2_Fp2,&tmp2_Fp2,&tmp3_Fp2);

	Fp2_sub(&f->x,&tmp1_Fp2,&tmp2_Fp2);
}
*/
/*
void ff_ltq_projective(FpJ12 *f,EFpJ12 *T,EFp12 *Q,EFpJ12 *P){
	static Fp12 tmp1_Fp12,tmp2_Fp12,tmp3_Fp12;
	//XpZt-YtZp
	Fp12_mul(&tmp1_Fp12,&P->z,&T->y);
	Fp12_mul(&tmp2_Fp12,&P->x,&T->z);
	Fp12_sub(&tmp2_Fp12,&tmp1_Fp12,&tmp2_Fp12);
	Fp12_mul(&f->z,&tmp2_Fp12,&P->z);

	//(ZpYq-Yp)(XpZt-YtZp)
	Fp12_mul(&tmp3_Fp12,&P->z,&Q->y);
	Fp12_sub(&tmp3_Fp12,&tmp3_Fp12,&P->y);
	Fp12_mul(&tmp3_Fp12,&tmp3_Fp12,&tmp2_Fp12);

	//YpZt-ZpYt (=YpZt-tmp1)
	Fp12_mul(&tmp2_Fp12,&P->y,&T->z);
	Fp12_sub(&tmp2_Fp12,&tmp2_Fp12,&tmp1_Fp12);
	//(Xq-Xp)(YpZt-Yt)
	Fp12_mul(&tmp1_Fp12,&P->z,&Q->x);
	Fp12_sub(&tmp1_Fp12,&tmp1_Fp12,&P->x);
	Fp12_mul(&tmp2_Fp12,&tmp1_Fp12,&tmp2_Fp12);

	Fp12_sub(&f->x,&tmp3_Fp12,&tmp2_Fp12);
}
*/
void EFp12_lineTP(Fp12 *rop,EFp12 *q,EFp12 *p,EFp12 *t){/*L_t,p(q)*/
  if(Fp12_cmp(&(p->x),&(t->x))==0){/*Xp=Xtの時の例外処理*/
    Fp12_sub(rop,&(q->x),&(t->x));
  }else if(q->infinity==1){
    //printf("\nTP::q.infinity==1\n");
    //getchar();
  }else if(p->infinity==1){
    //printf("\nTP::p.infinity==1\n");
    //getchar();
  }else if(t->infinity==1){
    //printf("\nTP::t.infinity==1\n");
    //getchar();
    /*add*/
    Fp12_set_ui(rop,1);
  }else{
    EFp12 temp;
    EFp12_init(&temp);
    Fp12 slope;
    Fp12_init(&slope);

    //傾き計算
    Fp12_sub(&(temp.y), &(p->y), &(t->y));
    Fp12_sub(&(temp.x), &(p->x), &(t->x));
    Fp12_inv(&temp.x,&temp.x);
    Fp12_mul(&slope, &(temp.y), &(temp.x));
    //temp.x=slope*(xq-xp)
    Fp12_sub(&(temp.x), &(q->x), &(p->x));
    Fp12_mul(&(temp.x), &slope, &(temp.x));
    //temp.y=(yq-yp)
    Fp12_sub(&(temp.y), &(q->y), &(p->y));
    Fp12_sub(rop, &(temp.y), &(temp.x));

    getchar();
    Fp12_println("L=",rop);
  }
}
void EFp12_lineTT(Fp12 *rop,EFp12 *q,EFp12 *t){/*l_T,T(Q)*/
  if(Fp12_cmp_ui(&(t->y),0)==0){
    /*
    printf("\n=====check line ttf cmp====\n");
    getchar();
    */
    Fp12_sub(rop,&(q->x),&(t->x));
  }else if(q->infinity==1){
    //printf("\nTT::q.infinity==1\n");
    //getchar();
  }else if(t->infinity==1){
    //printf("\nTT::t.infinity==1\n");
    //getchar();
    /*add*/
    Fp12_set_ui(rop,1);
  }else{
    EFp12 temp;
    EFp12_init(&temp);
    Fp12 slope;
    Fp12_init(&slope);

    //傾き計算
    //temp.x=3*xt^2+a
    Fp12_sqr(&(temp.x), &(t->x));
    Fp12_mul_ui(&(temp.x), &(temp.x), 3);
    //Fp12_add(&temp.x,&temp.x,&Fp2_a);
    //tenp.y=2*yt

    Fp12_mul_ui(&(temp.y), &(t->y), 2);
    //slope=temp.x/temp.y
    //Fp12_div(&slope, &(temp.x), &(temp.y));
    Fp12_inv(&temp.y,&temp.y);
    Fp12_mul(&slope, &(temp.x), &(temp.y));

    //temp.x=slope*(xq-xt)
    Fp12_sub(&(temp.x), &(q->x), &(t->x));
    Fp12_mul(&(temp.x), &slope, &(temp.x));

    //temp.y=yq-yt
    Fp12_sub(&(temp.y), &(q->y), &(t->y));
    Fp12_sub(rop,&(temp.y),&(temp.x));

    getchar();
    Fp12_println("L=",rop);
  }
}
void lineTT_projective_tst(Fp12 *rop,EFp2 *q,EFp2 *t){

}
