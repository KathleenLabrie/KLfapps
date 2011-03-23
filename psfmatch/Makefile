################################################################
#
# Makefile for psfmatch
#
# PSFmatch matches the PSF of a good seeing image to the PSF 
# of a worse seeing image.  PSFmatch was written by CJP and
# modified by KL (see modifications list in CHANGES).
#
# Coding languages: Fortran, C
# Libraries required: IRAF's IMFORT, PGPLOT, CJP, KLstr, KLimgio [Fortran]
#			 KLfio [Fortran], KLcfitsio, CFITSIO, KLutil [C]
# Compilers: g77, gcc
# OS: RedHat Linux, Mac OS X
#
################################################################

MAKE = make
GZIP = gzip -S .gz

#### Environment Variable Definitions ####
OS = macosx
PREFIX = /Users/klabrie/prgf
LIBDIR = $(PREFIX)/lib

C_PREFIX = /Users/klabrie/prgc
C_LIBDIR = $(C_PREFIX)/lib
C_INCLDIR = $(C_PREFIX)/include
IRAFDIR = /Local/iraf
XLIBDIR = /usr/X11R6/lib
#PGPLOT_DIR = /pel0/pritchet/pgplot

EPREFIX = $(PREFIX)
LPREFIX = $(EPREFIX)
IPREFIX = $(EPREFIX)
SPREFIX = $(EPREFIX)/packages/

#### Compilers Variable Definitions ####
FC = g77
FFLAGS = -g -O2 -I$(PGPLOT_DIR)
#CC = gcc
#CFLAGS = -g -O2 -I$(CINCLDIR) -I$(PGPLOT_DIR)

################ DO NOT EDIT BELOW THIS LINE ###################

#### Variable Definitions ####
PKGNAME = psfmatch
BASIC_FILES = Makefile Changes TODO
DOCS = $(PKGNAME).README
VERSION = 1
RELEASE = $(VERSION).1.4

#### Libraries ####
IMFORTLIBS = $(IRAFDIR)/irafbin/bin.$(OS)/libimfort.a \
             $(IRAFDIR)/irafbin/bin.$(OS)/libsys.a \
	      $(IRAFDIR)/irafbin/bin.$(OS)/libvops.a \
	      $(IRAFDIR)/iraf/unix/bin.$(OS)/libos.a
PGLIBS = $(PGPLOT_DIR)/libpgplot.a -L$(XLIBDIR) -lX11
FLIBS = -L$(LIBDIR) -lCJP -lKLstr -lKLimgio -lKLfio

CLIBS = -L$(C_LIBDIR) -lKLutil -lKLcfitsio -lcfitsio -lm 

#### Object Definitions ####
PSFOBJ =  psfmatch.o rmnsd_sampled.o subarray.o getkern.o
FIXOBJ = fixkernel.o rmnsd_sampled.o subarray.o getkern.o
ARROBJ = array2list.o
OBJECTS = $(PSFOBJ) $(FIXOBJ) $(ARROBJ)
SRC = $(OBJECTS:.o=.f)
BINS = psfmatch fixkernel array2list

#### Rule Definitions ####
#%.o : %.c
#	$(CC) -c $(CFLAGS) $<

%.o : %.f
	$(FC) -c $(FFLAGS) $<

main : array2list fixkernel $(PSFOBJ)
	$(FC) $(FFLAGS) -o psfmatch $(PSFOBJ) $(FLIBS) $(CLIBS) \
		$(IMFORTLIBS) $(PGLIBS)
all : main

array2list : $(ARROBJ)
	$(FC) $(FFLAGS) -o array2list $(ARROBJ) $(FLIBS)

fixkernel : $(FIXOBJ)
	$(FC) $(FFLAGS) -o fixkernel $(FIXOBJ) $(FLIBS) $(CLIBS) \
		$(IMFORTLIBS)

## Installing ##
.PHONY: install isntall-src
install :
	install -d -m 0755 $(EPREFIX)/bin
	@for bin in $(BINS) ; do \
	  install -m 0755 $$bin $(EPREFIX)/bin/; \
	done

install-src :
	install -d -m 0755 $(SPREFIX)/$(PKGNAME)
	@for file in $(SRC) ; do \
	  install -m 0644 $$file $(SPREFIX)/$(PKGNAME)/; \
	done
	@for file in $(BASIC_FILES) ; do \
	  install -m 0644 $$file $(SPREFIX)/$(PKGNAME)/; \
	done
	@for file in $(DOCS) ; do \
	  install -m 0644 $$file $(SPREFIX)/$(PKGNAME)/; \
	done

## Packaging ##
.PHONY: package-src
package-src :
	tar cvf $(PKGNAME)-$(RELEASE).tar $(BASIC_FILES) $(DOCS) \
	$(SRC)
	$(GZIP) $(PKGNAME)-$(RELEASE).tar

.PHONY: clean
clean :
	\rm -f $(BINS) $(OBJECTS)

### End of Makefile ###
