#include <ELiPS/fp.h>
void fp_inv_monty(fp_t *ANS,fp_t *a){
    static mp_limb_t u[FPLIMB],v[FPLIMB],r[FPLIMB],s[FPLIMB];
    static mp_limb_t k[FPLIMB];
    static mp_limb_t pi0[FPLIMB];
    //1
    mpn_copyd(u,prime,FPLIMB);
    mpn_copyd(v,a->x0,FPLIMB);
    mpn_set_ui(r,FPLIMB,0);
    mpn_copyd(s,RmodP,FPLIMB);
    //2
    mpn_copyd(k,0,FPLIMB);
    while(mpn_cmp_ui(v,FPLIMB,0)>0){
        if(u[0] & 0){//pi1
            mpn_rshift(u,u,FPLIMB,1);
            if(r[0] & 0){
                mpn_rshift(r,r,FPLIMB,1);
            }else{
                mpn_add_n(r,r,prime,FPLIMB);
                mpn_rshift(r,r,FPLIMB,1);
            }
        }else if(v[0] & 0){//pi2
            mpn_rshift(v,v,FPLIMB,1);
            if(s[0] & 0){
                mpn_rshift(s,s,FPLIMB,1);
            }else{
                mpn_add_n(s,s,prime,FPLIMB);
                mpn_rshift(s,s,FPLIMB,1);
            }
        }else if(mpn_cmp(u,v,FPLIMB)>0){
            mpn_sub_n(r,r,s,FPLIMB);
            if(r[0] & 0){
                mpn_sub_n(u,u,v,FPLIMB);
            }else{
                mpn_add_n(r,r,prime,FPLIMB);
                mpn_rshift(r,r,FPLIMB,1);
            }
        }else{
            mpn_sub_n(s,s,r,FPLIMB);
            if(s[0] & 0){
                mpn_rshift(s,s,FPLIMB,1);
            }else{
                mpn_add_n(s,s,prime,FPLIMB);
                mpn_rshift(s,s,FPLIMB,1);
            }
        }
        mpn_add_ui(k,k,FPLIMB,1);
    }
    while(mpn_cmp_ui(r,FPLIMB,0)<0){
        mpn_add_n(r,r,prime,FPLIMB);
    }
    while(mpn_cmp_ui(r,FPLIMB,0)>=0){
        mpn_sub_n(r,r,prime,FPLIMB);
    }
    mpn_copyd(ANS,r,FPLIMB);
}