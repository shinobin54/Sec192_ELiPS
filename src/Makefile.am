
lib_LTLIBRARIES = libelips.la

libelips_la_LIBADD = bls12/bls12_generate_curve.lo\
bls12/bls12_generate_parameter.lo\
bls12/bls12_init.lo

SUBDIRS = bls12

AM_CPPFLAGS = -I../include

LDADD = bls12/libelipsbls12.la

libelips_la_SOURCES = define.c\
scalar.c\
mpn.c\
fp.c\
fp2.c\
fp4.c\
efp.c\
efp2.c\
efp4.c\
time.c\
../include/ELiPS/define.h\
../include/ELiPS/scalar.h\
../include/ELiPS/mpn.h\
../include/ELiPS/fp.h\
../include/ELiPS/fp2.h\
../include/ELiPS/fp4.h\
../include/ELiPS/efp.h\
../include/ELiPS/efp2.h\
../include/ELiPS/efp4.h\
../include/ELiPS/time.h
