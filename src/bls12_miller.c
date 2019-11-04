#include <ELiPS/bls12_miller.h>
//bls12
void BLS12_Miller_algo_for_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
  EFp12 Test;
  EFp12_init(&Test);
  EFp2 T;
  EFp2_init(&T);
  EFp2 mapped_Q;
  EFp2_init(&mapped_Q);
  EFp mapped_P;
  EFp_init(&mapped_P);
  Fp12 f;
  Fp12_init(&f);
  Fp L;
  Fp_init(&L);
  mpz_t loop;
  mpz_init(loop);
  mpz_sub_ui(loop,trace_z,1);
  int i,length;
  length=(int)mpz_sizeinbase(loop,2);
  char binary[length];
  mpz_get_str(binary,2,loop);

  EFp12_to_EFp(&mapped_P,P);
  EFp12_to_EFp2(&mapped_Q,Q);
  Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
  EFp2_set(&T,&mapped_Q);     //set T
  Fp_set_ui(&f.x0.x0.x0,1);   //set f

  //miller
  for(i=1; binary[i]!='\0'; i++){
    ff_ltt(&f,&T,&mapped_P,&L);
    if(binary[i]=='1'){
      f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
    }
  }

  Fp12_set(ANS,&f);


  mpz_clear(loop);
}

void BLS12_Miller_algo_for_opt_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
  EFp12 Buf;
  EFp12_init(&Buf);
  EFp2 T;
  EFp2_init(&T);
  EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
  EFp2_init(&mapped_Q);
  EFp2_init(&mapped_Q_neg);
  EFp2_init(&mapped_Q1);
  EFp2_init(&mapped_Q2_neg);
  EFp mapped_P;
  EFp_init(&mapped_P);
  Fp12 f;
  Fp12_init(&f);
  Fp L;
  Fp_init(&L);
  int i;

  //set
  EFp12_to_EFp(&mapped_P,P);//set mapped_P
  EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
  Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
  EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
  Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
  EFp2_set(&T,&mapped_Q);     //set T
  Fp_set_ui(&f.x0.x0.x0,1);
  //miller
  for(i=BLS12_X_length-1; i>=0; i--){
    switch(BLS12_X_binary[i]){
      case 0:
      ff_ltt(&f,&T,&mapped_P,&L);
      break;
      case 1:
      ff_ltt(&f,&T,&mapped_P,&L);
      f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
      break;
      case -1:
      ff_ltt(&f,&T,&mapped_P,&L);
      f_ltq(&f,&T,&mapped_Q_neg,&mapped_P,&L);
      break;
      default:
      break;
    }
  }


  Fp12_set(ANS,&f);
}
void BLS12_Miller_algo_for_opt_ate_lazy(Fp12 *ANS,EFp12 *P,EFp12 *Q){
  EFp12 Buf;
  EFp12_init(&Buf);
  EFp2 T;
  EFp2_init(&T);
  EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
  EFp2_init(&mapped_Q);
  EFp2_init(&mapped_Q_neg);
  EFp2_init(&mapped_Q1);
  EFp2_init(&mapped_Q2_neg);
  EFp mapped_P;
  EFp_init(&mapped_P);
  Fp12 f;
  Fp12_init(&f);
  Fp L;
  Fp_init(&L);
  int i;

  //set
  EFp12_to_EFp(&mapped_P,P);//set mapped_P
  EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
  Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
  EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
  Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
  EFp2_set(&T,&mapped_Q);     //set T
  Fp_set_ui(&f.x0.x0.x0,1);
  //miller
  for(i=BLS12_X_length-1; i>=0; i--){
    switch(BLS12_X_binary[i]){
      case 0:
      ff_ltt_lazy(&f,&T,&mapped_P,&L);
      break;
      case 1:
      ff_ltt_lazy(&f,&T,&mapped_P,&L);
      f_ltq_lazy(&f,&T,&mapped_Q,&mapped_P,&L);
      break;
      case -1:
      ff_ltt_lazy(&f,&T,&mapped_P,&L);
      f_ltq_lazy(&f,&T,&mapped_Q_neg,&mapped_P,&L);
      break;
      default:
      break;
    }
  }

  Fp12_set(ANS,&f);
}
void BLS12_Miller_algo_for_opt_ate_lazy_montgomery(Fp12 *ANS,EFp12 *P,EFp12 *Q){
  EFp12 Buf;
  EFp12_init(&Buf);
  EFp2 T;
  EFp2_init(&T);
  EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
  EFp2_init(&mapped_Q);
  EFp2_init(&mapped_Q_neg);
  EFp2_init(&mapped_Q1);
  EFp2_init(&mapped_Q2_neg);
  EFp mapped_P;
  EFp_init(&mapped_P);
  Fp12 f;
  Fp12_init(&f);
  Fp L;
  Fp_init(&L);
  int i;

  //set
  EFp12_to_EFp(&mapped_P,P);//set mapped_P
  EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
  EFp_to_montgomery(&mapped_P,&mapped_P);
  EFp2_to_montgomery(&mapped_Q,&mapped_Q);

  Pseudo_8_sparse_mapping_montgomery(&mapped_P,&mapped_Q,&L);
  EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
  Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
  EFp2_set(&T,&mapped_Q);     //set T
  Fp_set_ui(&f.x0.x0.x0,1);
  //miller
  for(i=BLS12_X_length-1; i>=0; i--){
    switch(BLS12_X_binary[i]){
      case 0:
      ff_ltt_lazy_montgomery(&f,&T,&mapped_P,&L);
      break;
      case 1:
      ff_ltt_lazy_montgomery(&f,&T,&mapped_P,&L);
      f_ltq_lazy_montgomery(&f,&T,&mapped_Q,&mapped_P,&L);
      break;
      case -1:
      ff_ltt_lazy_montgomery(&f,&T,&mapped_P,&L);
      f_ltq_lazy_montgomery(&f,&T,&mapped_Q_neg,&mapped_P,&L);
      break;
      default:
      break;
    }
  }

  Fp12_mod_montgomery(ANS,&f);
}

/*                              mine                                    */
/*
void BLS12_Miller_algo_for_plain_ate_homogeneou(Fp12 *ANS,EFp12 *P,EFp12 *Q){
EFp12 Test;
EFp12_init(&Test);
EFp12 T;
EFp12_init(&T);
Fp12 f;
Fp12_init(&f);

mpz_t loop;
mpz_init(loop);
mpz_sub_ui(loop,trace_z,1);
int i,length;
length=(int)mpz_sizeinbase(loop,2);
char binary[length];
mpz_get_str(binary,2,loop);


//miller
for(i=1; binary[i]!='\0'; i++){
ff_ltt(&f,&T,&mapped_P,&L);
if(binary[i]=='1'){
f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
}
}

Fp12_set(ANS,&f);


mpz_clear(loop);
}
*/
/*
void Miller_opt_projective(Fp12 *ANS,EFp12 *P,EFp12 *Q){
printf("start MILLER\n");
FpJ12 f,miller_tmp1;
EFpJ12 JT,JP;
FpJ12_init(&f);
FpJ12_init(&miller_tmp1);
EFpJ12_init(&JT);
EFpJ12_init(&JP);
//f<-1,T<-P
Fp_set_ui(&f.x.x0.x0.x0,1);
EFp12_to_EFpJ12(&JP,P);
EFp12_to_EFpJ12(&JT,Q);


mpz_t loop;
mpz_init(loop);
mpz_sub_ui(loop,trace_z,1);
int i,length;
length=(int)mpz_sizeinbase(loop,2);
char binary[length];
mpz_get_str(binary,2,loop);


for(i=BLS12_X_length-1; i>=0; i--){
FpJ12_sqr(&f,&f);
ff_ltt_projective_Fp12(&miller_tmp1,&JT,Q);
FpJ12_mul(&f,&f,&miller_tmp1);
EFp12_ECD_Jacobian(&JT,&JT);
if(JT.infinity==1) printf("error::JT=0(ECD)\n");
if(binary[i]=='1'){
//void EFp12_lineTP(Fp12 *rop,EFp12 *q,EFp12 *p,EFp12 *t)
ff_ltq_projective_Fp12(&miller_tmp1,&JT,Q,&JP);
FpJ12_mul(&f,&f,&miller_tmp1);
EFp12_ECA_Jacobian(&JT,&JT,&JP);
if(JT.infinity==1) printf("error::JT=0(ECD)\n");
}
}

if(JT.infinity==1) printf("error::JT=0(end)\n");
Fp12_inv(&f.z,&f.z);
Fp12_mul(ANS,&f.x,&f.z);

Fp12_clear(&f);
Fp12J_clear(&miller_tmp1);
EFpJ12_clear(&JT);
EFpJ12_clear(&JP);

}
*/

void Miller_debag(Fp12 *ANS,EFp12 *P,EFp12 *Q){
  EFp12 T;
  EFp12_init(&T);
  Fp12 f,Fp12_temp2;
  Fp12_init(&f);
  Fp12_init(&Fp12_temp2);
  mpz_t loop;
  mpz_init(loop);
  mpz_sub_ui(loop,trace_z,1);
  int i,length;
  length=(int)mpz_sizeinbase(loop,2);
  char binary[length];
  mpz_get_str(binary,2,loop);

  EFp12_set(&T,Q);     //set T
  Fp_set_ui(&f.x0.x0.x0,1);   //set f

  //miller
  for(i=1; binary[i]!='\0'; i++){
    Fp12_sqr(&f,&f);

    getchar();
    Fp12_println("test(f^2 DBL)=",&f);

    //void EFp2_lineTT(Fp2 *rop,EFp2 *q,EFp2 *t)
    EFp12_lineTT(&Fp12_temp2,P,&T);

    Fp12_mul(&f,&f,&Fp12_temp2);
    EFp12_ECD(&T,&T);

    getchar();
    Fp12_println("test(DBL)=",&f);
    if(binary[i]=='1'){
      EFp12_lineTP(&Fp12_temp2,P,Q,&T);
      Fp12_mul(&f,&f,&Fp12_temp2);
      EFp12_ECA(&T,&T,Q);

      getchar();
      Fp12_println("test(ADD)=",&f);
    }
  }
  Fp12_set(ANS,&f);


  mpz_clear(loop);
}
void Miller_opt_projective_plain(Fp12 *ANS,EFp12 *P,EFp12 *Q){
  EFp12 T;
  EFp12_init(&T);
  FpJ12 f,miller_tmp1;
  FpJ12_init(&f);
  FpJ12_init(&miller_tmp1);

  Fp12 test;
  Fp12_init(&test);

  mpz_t loop;
  mpz_init(loop);
  mpz_sub_ui(loop,trace_z,1);
  int i,length;
  length=(int)mpz_sizeinbase(loop,2);
  char binary[length];
  mpz_get_str(binary,2,loop);

  EFp12_set(&T,Q);     //set T
  FpJ12_set_ui(&f,1);   //set f

  //miller
  for(i=1; binary[i]!='\0'; i++){
    FpJ12_sqr(&f,&f);

    ff_ltt_projective_plain(&miller_tmp1,&T,P);
    FpJ12_mul(&f,&f,&miller_tmp1);
    EFp12_ECD(&T,&T);
/*
    getchar();
    FpJ12_to_Fp12(&test,&f);
    Fp12_println("test(DBL)=",&test);
*/
    if(binary[i]=='1'){
      ff_ltq_projective_plain(&miller_tmp1,&T,P,Q);
      FpJ12_mul(&f,&f,&miller_tmp1);
      EFp12_ECA(&T,&T,Q);
/*
      getchar();
      FpJ12_to_Fp12(&test,&f);
      Fp12_println("test(ADD)=",&test);
*/
    }
  }
  //FpJ12_to_Fp12(ANS,&f);
  Fp12_set(ANS,&f.x);
  printf("test start\n");
  Fp12 uni_test;
  Fp12_init(&uni_test);
  BLS12_Final_exp_plain(&uni_test,&f.z);
  Fp12_println("unitest=",&uni_test);

  Fp12_printf("ANS=",ANS);
  mpz_clear(loop);
}
void Miller_opt_projective_projective(Fp12 *ANS,EFp12 *P,EFp12 *Q){
  EFpJ12 JT,JQ;
  EFpJ12_init(&JT);
  EFpJ12_init(&JQ);
  FpJ12 f,miller_tmp1;
  FpJ12_init(&f);
  FpJ12_init(&miller_tmp1);

  Fp12 test;
  Fp12_init(&test);

  mpz_t loop;
  mpz_init(loop);
  mpz_sub_ui(loop,trace_z,1);
  int i,length;
  length=(int)mpz_sizeinbase(loop,2);
  char binary[length];
  mpz_get_str(binary,2,loop);

  EFp12_to_EFpJ12(&JT,Q);     //set T
  EFp12_to_EFpJ12(&JQ,Q);
  FpJ12_set_ui(&f,1);   //set f

  //miller
  for(i=1; binary[i]!='\0'; i++){
    FpJ12_sqr(&f,&f);

    getchar();
    FpJ12_to_Fp12(&test,&f);
    Fp12_println("test(f^2DBL)=",&test);

    ff_ltt_projective_projective(&miller_tmp1,&JT,P);
    FpJ12_mul(&f,&f,&miller_tmp1);
    EFp12_ECD_Jacobian(&JT,&JT);

    getchar();
    FpJ12_to_Fp12(&test,&f);
    Fp12_println("test(DBL)=",&test);

    if(binary[i]=='1'){
      ff_ltq_projective_projective(&miller_tmp1,&JT,P,&JQ);
      FpJ12_mul(&f,&f,&miller_tmp1);
      EFp12_ECA_Jacobian(&JT,&JT,&JQ);

      getchar();
      FpJ12_to_Fp12(&test,&f);
      Fp12_println("test(ADD)=",&test);
    }
  }
  FpJ12_to_Fp12(ANS,&f);
  Fp12_printf("ANS=",&f);
  mpz_clear(loop);
}
/*
void Miller_opt_projective(Fp12 *ANS,EFp12 *P,EFp12 *Q){
  FpJ12 f,miller_tmp1;
  int i;
  EFpJ12 JT,JP;
  FpJ12_init(&f);
  FpJ12_init(&miller_tmp1);
  EFpJ12_init(&JT);
  EFpJ12_init(&JP);
  //f<-1,T<-P
  mpz_t loop;
  mpz_init(loop);
  mpz_sub_ui(loop,trace_z,1);

  Fp_set_ui(&f.x.x0.x0.x0,1);
  EFp12_to_EFpJ12(&JP,P);
  EFp12_to_EFpJ12(&JT,Q);


  for(i=(int)(mpz_sizeinbase(loop,2))-2;i>=0;i--){
    //getchar();
    FpJ12_sqr(&f,&f);
    ff_ltt_projective_Fp12(&miller_tmp1,&JT,P);
    FpJ12_mul(&f,&f,&miller_tmp1);
    EFp12_ECD_Jacobian(&JT,&JT);
    //Fp12_printf("\n\nf(ltt)x=",&f.x);
    //Fp12_printf("\nf(ltt)z=",&f.z);

    if(mpz_tstbit(loop,i)==1){
      //void EFp12_lineTP(Fp12 *rop,EFp12 *q,EFp12 *p,EFp12 *t)
      ff_ltq_projective_Fp12(&miller_tmp1,&JT,Q,&JP);
      FpJ12_mul(&f,&f,&miller_tmp1);
      EFp12_ECA_Jacobian(&JT,&JT,&JP);

      //Fp12_printf("\n\nf(ltq)x=",&f.x);
      //Fp12_printf("\nf(ltq)z=",&f.z);
    }
  }
  Fp12_inv(&f.z,&f.z);
  Fp12_mul(ANS,&f.x,&f.z);
}
*/
/*
void BLS12_Miller_algo_for_projective(Fp12 *ANS,EFp12 *P,EFp12 *Q){
EFp12 Test;
EFp12_init(&Test);
EFp12 T;
EFp12_init(&T);
EFpJ12 JP,JT;

Fp12 f;
Fp12_init(&f);
mpz_t loop;
mpz_init(loop);
mpz_sub_ui(loop,trace_z,1);
int i,length;
length=(int)mpz_sizeinbase(loop,2);
char binary[length];
mpz_get_str(binary,2,loop);
EFp12_set(&T,Q);     //set T
Fp_set_ui(&f.x0.x0.x0,1);   //set f

//miller
for(i=1; binary[i]!='\0'; i++){
ff_ltt_projective(&f,&T,P);
if(binary[i]=='1'){
f_ltq_projective(&f,&T,Q,P);
}
}

Fp12_set(ANS,&f);


mpz_clear(loop);
}
*/
