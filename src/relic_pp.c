void pp_dbl_k12_projc_basic(fp12_t l, ep2_t r, ep2_t q, ep_t p) {
	fp2_t t0, t1, t2, t3, t4, t5, t6;
	int one = 1, zero = 0;

	fp2_null(t0);
	fp2_null(t1);
	fp2_null(t2);
	fp2_null(t3);
	fp2_null(t4);
	fp2_null(t5);
	fp2_null(t6);

	TRY {
		fp2_new(t0);
		fp2_new(t1);
		fp2_new(t2);
		fp2_new(t3);
		fp2_new(t4);
		fp2_new(t5);
		fp2_new(t6);

		if (ep2_curve_is_twist() == EP_MTYPE) {
			one ^= 1;
			zero ^= 1;
		}

		if (ep_curve_opt_b() == RLC_TWO) {
			/* t0 = x1^2. */
			fp2_sqr(t0, q->x);
			/* t1 = B = y1^2. */
			fp2_sqr(t1, q->y);
			/* t2 = C = z1^2. */
			fp2_sqr(t2, q->z);
			/* t4 = A = (x1 * y1)/2. */
			fp2_mul(t4, q->x, q->y);
			fp_hlv(t4[0], t4[0]);
			fp_hlv(t4[1], t4[1]);
			/* t3 = E = 3b'C = 3C * (1 - i). */
			fp2_dbl(t3, t2);
			fp2_add(t2, t2, t3);
			fp_add(t3[0], t2[0], t2[1]);
			fp_sub(t3[1], t2[1], t2[0]);
			/* t2 = F = 3E. */
			fp2_dbl(t2, t3);
			fp2_add(t2, t3, t2);
			/* x3 = A * (B - F). */
			fp2_sub(r->x, t1, t2);
			fp2_mul(r->x, r->x, t4);
			/* t2 = G = (B + F)/2. */
			fp2_add(t2, t1, t2);
			fp_hlv(t2[0], t2[0]);
			fp_hlv(t2[1], t2[1]);
			/* y3 = G^2 - 3E^2. */
			fp2_sqr(t2, t2);
			fp2_sqr(t4, t3);
			/* t5 = y1 * z1. */
			fp2_mul(t5, q->y, q->z);
			fp2_dbl(r->y, t4);
			fp2_add(r->y, r->y, t4);
			fp2_sub(r->y, t2, r->y);

			/* t2 = H = 2 * y1 * z1. */
			fp2_dbl(t2, t5);
			/* z3 = B * H. */
			fp2_mul(r->z, t1, t2);

			/* l11 = E - B. */
			fp2_sub(l[one][one], t3, t1);

			/* l10 = (3 * xp) * t0. */
			fp_mul(l[one][zero][0], p->x, t0[0]);
			fp_mul(l[one][zero][1], p->x, t0[1]);

			/* l00 = H * (-yp). */
			fp_mul(l[zero][zero][0], t2[0], p->y);
			fp_mul(l[zero][zero][1], t2[1], p->y);
		} else {
			/* A = x1^2. */
			fp2_sqr(t0, q->x);
			/* B = y1^2. */
			fp2_sqr(t1, q->y);
			/* C = z1^2. */
			fp2_sqr(t2, q->z);
			/* D = 3bC, general b. */
			fp2_dbl(t3, t2);
			fp2_add(t3, t3, t2);
			ep2_curve_get_b(t4);
			fp2_mul(t3, t3, t4);
			/* E = (x1 + y1)^2 - A - B. */
			fp2_add(t4, q->x, q->y);
			fp2_sqr(t4, t4);
			fp2_sub(t4, t4, t0);
			fp2_sub(t4, t4, t1);

			/* F = (y1 + z1)^2 - B - C. */
			fp2_add(t5, q->y, q->z);
			fp2_sqr(t5, t5);
			fp2_sub(t5, t5, t1);
			fp2_sub(t5, t5, t2);

			/* G = 3D. */
			fp2_dbl(t6, t3);
			fp2_add(t6, t6, t3);

			/* x3 = E * (B - G). */
			fp2_sub(r->x, t1, t6);
			fp2_mul(r->x, r->x, t4);

			/* y3 = (B + G)^2 -12D^2. */
			fp2_add(t6, t6, t1);
			fp2_sqr(t6, t6);
			fp2_sqr(t2, t3);
			fp2_dbl(r->y, t2);
			fp2_dbl(t2, r->y);
			fp2_dbl(r->y, t2);
			fp2_add(r->y, r->y, t2);
			fp2_sub(r->y, t6, r->y);

			/* z3 = 4B * F. */
			fp2_dbl(r->z, t1);
			fp2_dbl(r->z, r->z);
			fp2_mul(r->z, r->z, t5);

			/* l11 = D - B. */
			fp2_sub(l[one][one], t3, t1);

			/* l10 = (3 * xp) * A. */
			fp_mul(l[one][zero][0], p->x, t0[0]);
			fp_mul(l[one][zero][1], p->x, t0[1]);

			/* l00 = F * (-yp). */
			fp_mul(l[zero][zero][0], t5[0], p->y);
			fp_mul(l[zero][zero][1], t5[1], p->y);
		}
		r->norm = 0;
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp2_free(t0);
		fp2_free(t1);
		fp2_free(t2);
		fp2_free(t3);
		fp2_free(t4);
		fp2_free(t5);
		fp2_free(t6);
	}
}

void pp_add_k12_projc_basic(fp12_t l, ep2_t r, ep2_t q, ep_t p) {
	fp2_t t0, t1, t2, t3, t4;
	int one = 1, zero = 0;

	fp2_null(t0);
	fp2_null(t1);
	fp2_null(t2);
	fp2_null(t3);
	fp2_null(t4);

	TRY {
		fp2_new(t0);
		fp2_new(t1);
		fp2_new(t2);
		fp2_new(t3);
		fp2_new(t4);

		/* B = t0 = x1 - x2 * z1. */
		fp2_mul(t0, r->z, q->x);
		fp2_sub(t0, r->x, t0);
		/* A = t1 = y1 - y2 * z1. */
		fp2_mul(t1, r->z, q->y);
		fp2_sub(t1, r->y, t1);

		/* D = B^2. */
		fp2_sqr(t2, t0);
		/* G = x1 * D. */
		fp2_mul(r->x, r->x, t2);
		/* E = B^3. */
		fp2_mul(t2, t2, t0);
		/* C = A^2. */
		fp2_sqr(t3, t1);
		/* F = E + z1 * C. */
		fp2_mul(t3, t3, r->z);
		fp2_add(t3, t2, t3);

		if (ep2_curve_is_twist() == EP_MTYPE) {
			one ^= 1;
			zero ^= 1;
		}

		/* l10 = - (A * xp). */
		fp_mul(l[one][zero][0], t1[0], p->x);
		fp_mul(l[one][zero][1], t1[1], p->x);
		fp2_neg(l[one][zero], l[one][zero]);

		/* t4 = B * x2. */
		fp2_mul(t4, q->x, t1);

		/* H = E + F - 2 * G. */
		fp2_sub(t3, t3, r->x);
		fp2_sub(t3, t3, r->x);
		/* y3 = A * (G - H) - y1 * E. */
		fp2_sub(r->x, r->x, t3);
		fp2_mul(t1, t1, r->x);
		fp2_mul(r->y, t2, r->y);
		fp2_sub(r->y, t1, r->y);
		/* x3 = B * H. */
		fp2_mul(r->x, t0, t3);
		/* z3 = z1 * E. */
		fp2_mul(r->z, r->z, t2);

		/* l11 = J = B * x2 - A * y2. */
		fp2_mul(t2, q->y, t0);
		fp2_sub(l[one][one], t4, t2);

		/* l00 = B * yp. */
		fp_mul(l[zero][zero][0], t0[0], p->y);
		fp_mul(l[zero][zero][1], t0[1], p->y);

		r->norm = 0;
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		fp2_free(t0);
		fp2_free(t1);
		fp2_free(t2);
		fp2_free(t3);
		fp2_free(t4);
	}
}

void pp_map_oatep_k12(fp12_t r, ep_t p, ep2_t q) {
	ep_t _p[1];
	ep2_t t[1], _q[1];
	bn_t a;

	ep_null(_p[0]);
	ep2_null(_q[0]);
	ep2_null(t[0]);
	bn_null(a);

	TRY {
		ep_new(_p[0]);
		ep2_new(_q[0]);
		ep2_new(t[0]);
		bn_new(a);

		fp_prime_get_par(a);
		fp12_set_dig(r, 1);

		ep_norm(_p[0], p);
		ep2_norm(_q[0], q);

		if (!ep_is_infty(_p[0]) && !ep2_is_infty(_q[0])) {
			switch (ep_curve_is_pairf()) {
				case EP_BN:
					bn_mul_dig(a, a, 6);
					bn_add_dig(a, a, 2);
					/* r = f_{|a|,Q}(P). */
					pp_mil_k12(r, t, _q, _p, 1, a);
					if (bn_sign(a) == RLC_NEG) {
						/* f_{-a,Q}(P) = 1/f_{a,Q}(P). */
						fp12_inv_cyc(r, r);
						ep2_neg(t[0], t[0]);
					}
					pp_fin_k12_oatep(r, t[0], _q[0], _p[0]);
					pp_exp_k12(r, r);
					break;
				case EP_B12:
					/* r = f_{|a|,Q}(P). */
					pp_mil_k12(r, t, _q, _p, 1, a);
					if (bn_sign(a) == RLC_NEG) {
						fp12_inv_cyc(r, r);
						ep2_neg(t[0], t[0]);
					}
					pp_exp_k12(r, r);
					break;
			}
		}
	}
	CATCH_ANY {
		THROW(ERR_CAUGHT);
	}
	FINALLY {
		ep_free(_p[0]);
		ep2_free(_q[0]);
		ep2_free(t[0]);
		bn_free(a);
	}
}
