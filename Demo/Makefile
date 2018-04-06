# -*- Makefile -*-
#
# Standalone Makefile for Demo
# Demo/Makefile.sa
# In the 'nicaea' package, this file is renamed to 'Demo/Makefile'
# Martin Kilbinger 2008 - 2014

CC         = gcc

FFTW       = /usr/local
GSL        = $(HOME)

WARNINGS   = -Wall -Wuninitialized -pedantic
OPTIM      = -O3
FPIC	   = -fPIC
ARCHFLAGS  = -std=gnu9x

ifeq ("$(shell uname)", "Darwin")
   ARCH         = MacOSX
endif

ifeq ($(ARCH),MacOSX)
   EXT_DYN  = dylib
   SHARED  =  -macosx_version_min 10.7 # -dynamiclib -single_module
else
   EXT_DYN  = so
   SHARED  = -shared
endif

# Shared
#EXT = $(EXT_DYN)

# Static
EXT = a



cflags     = -I$(FFTW)/include -I$(GSL)/include $(OPTIM) $(WARNINGS) $(ARCHFLAGS) $(FPIC)
lflags     = -L$(FFTW)/lib -L$(GSL)/lib -lm -lgsl -lgslcblas -lfftw3

#objs      := $(patsubst %.h,%.o,$(wildcard *.h))
objs       = cosmo.o lensing.o reduced_fit.o decomp_eb.o sn1a.o nofz.o lensing_3rd.o \
             constants.o fr_constants.o emu.o fr_emu.o hubble.o \
	     cmb_bao.o halomodel.o hod.o errorlist.o io.o maths.o mvdens.o config.o par.o
exec       = lensingdemo halomodeldemo decomp_eb_demo sn1ademo cmb_bao_demo \
	     third_order_demo cosebi_demo


all: links_check libnicaea.$(EXT)

links_check:
	@ if ! test -e cosmo.c ; \
	then make links; \
	fi
	make $(exec)

#	then echo "Type 'make links' then 'make' to compile demo programs"; \

lensingdemo: lensingdemo.o libnicaea.$(EXT)
	$(CC) -o $@ $^ $(lflags) -lnicaea -L.

halomodeldemo: halomodeldemo.o libnicaea.$(EXT)
	$(CC) -o $@ $^ $(lflags) -lnicaea -L.

decomp_eb_demo: decomp_eb_demo.o libnicaea.$(EXT)
	$(CC) -o $@ $^ $(lflags) -lnicaea -L.

cosebi_demo: cosebi_demo.o libnicaea.$(EXT)
	$(CC) -o $@ $^ $(lflags) -lnicaea -L.

sn1ademo: sn1ademo.o libnicaea.$(EXT)
	$(CC) -o $@ $^ $(lflags) -lnicaea -L.

cmb_bao_demo: cmb_bao_demo.o libnicaea.$(EXT)
	$(CC) -o $@ $^ $(lflags) -lnicaea -L.

third_order_demo: third_order_demo.o libnicaea.$(EXT)
	$(CC) -o $@ $^ $(lflags) -lnicaea -L.

libnicaea.$(EXT_DYN): $(objs)
	$(LD) $(SHARED) -o $ $@ $^ $(lflags)

libnicaea.a: $(objs)
	$(AR) -r $@ $^


%.o:	%.c
	$(CC) $< -c $(cflags) -o $@

.h:



PHONY: links clean linksclean

links:
	ln -sf ../Cosmo/src/cosmo.c
	ln -sf ../Cosmo/src/lensing.c
	ln -sf ../Cosmo/src/reduced_fit.c
	ln -sf ../Cosmo/src/decomp_eb.c
	ln -sf ../Cosmo/src/sn1a.c
	ln -sf ../Cosmo/src/nofz.c
	ln -sf ../Cosmo/src/cmb_bao.c
	ln -sf ../Cosmo/src/lensing_3rd.c
	ln -sf ../Coyote/src/constants.c
	ln -sf ../Coyote/src/fr_constants.c
	ln -sf ../Coyote/src/emu.c
	ln -sf ../Coyote/src/fr_emu.c
	ln -sf ../Coyote/src/hubble.c
	ln -sf ../halomodel/src/halomodel.c
	ln -sf ../halomodel/src/hod.c
	ln -sf ../tools/src/errorlist.c
	ln -sf ../tools/src/io.c
	ln -sf ../tools/src/maths.c
	ln -sf ../tools/src/config.c
	ln -sf ../tools/src/mvdens.c
	ln -sf ../tools/src/par.c

	ln -sf ../Cosmo/include/cosmo.h
	ln -sf ../Cosmo/include/lensing.h
	ln -sf ../Cosmo/include/reduced_fit.h
	ln -sf ../Cosmo/include/decomp_eb.h
	ln -sf ../Cosmo/include/sn1a.h
	ln -sf ../Cosmo/include/nofz.h
	ln -sf ../Cosmo/include/cmb_bao.h
	ln -sf ../Cosmo/include/lensing_3rd.h
	ln -sf ../Coyote/include/constants.h
	ln -sf ../Coyote/include/fr_constants.h
	ln -sf ../Coyote/include/coyote.h
	ln -sf ../halomodel/include/halomodel.h
	ln -sf ../halomodel/include/hod.h
	ln -sf ../tools/include/errorlist.h
	ln -sf ../tools/include/io.h
	ln -sf ../tools/include/maths.h
	ln -sf ../tools/include/maths_base.h
	ln -sf ../tools/include/config.h
	ln -sf ../tools/include/mvdens.h
	ln -sf ../tools/include/par.h

	ln -sf ../par_files/cosmo.par
	ln -sf ../par_files/cosmoDP.par
	ln -sf ../par_files/cosmo_3rd.par
	ln -sf ../par_files/cosmo_SN.par
	ln -sf ../par_files/cosmo_lens.par
	ln -sf ../par_files/halomodel.par
	ln -sf ../par_files/nofz.par
	ln -sf ../par_files/nofz_1
	ln -sf ../par_files/nofz_2
	ln -sf ../par_files/nofz_hm.dat
	ln -sf ../par_files/nofz_single

linksclean:
	@ rm -f cosmo.{c,h} lensing.{c,h} reduced_fit.{c,h} decomp_eb.{c,h} \
	sn1a.{c,h} nofz.{c,h} cmb_bao.{c,h} lensing_3rd.{c,h} constants.{c,h} \
	fr_constants.{c,h} emu.c fr_emu.cd hubble.c coyote.h \
	halomodel.{c,h} hod.{c,h} errorlist.{c,h} io.{c,h} \
	maths.{c,h} config.{c,h} mvdens.{c,h} maths_base.h par.{c,h} \
	cosmo.par cosmo_lens.par halomodel.par cosmo_SN.par cosmo_DP.par

clean:
	rm -rf *~ *.o lib*.a *.oct lib*.so lib*.dylib $(exec)

