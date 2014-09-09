# Generic instructions:  make, then as root, make install.
# To uninstall, as root, make uninstall .
# To remove executable and objects from build directory, make clean .
# You can edit this file to change the choice of compiler CC, 
# install directories
# FLCT_BINDIR and FLCT_MANDIR, compiler options COPTS, or location of
# fftw3 library LIBFFTW3, or location of FFTW include files INCLUDEFFTW3 . 
# If no root permission, will need to change these locations. 
#CC = gcc4
#CC = gcc
CC = gcc
#CC = tcc
# linux:
#FLCT_MANDIR = /usr/local/man/man1
# macosx:
#FLCT_MANDIR = /usr/share/man/man1
FLCT_MANDIR = /usr/local/man/man1
FLCT_BINDIR = /usr/local/bin
LIBFFTW3 = /usr/local/lib
LOPTS = -lfftw3 -lm
INCLUDEFFTW3 = /usr/local/include
# next version of options for a reasonably optimized executable
COPTS = -O3 -Wall -fomit-frame-pointer
# next version of options for debugging with gdb or valgrind
#COPTS = -g
flct:	flct.c 
	$(CC) $(COPTS) flct.c -I$(INCLUDEFFTW3) -L$(LIBFFTW3) $(LOPTS) -o flct
warp:	warp.c
	$(CC) $(COPTS) warp.c -I$(INCLUDEFFTW3) -L$(LIBFFTW3) $(LOPTS) -o warp
install: flct
	cp ./flct $(FLCT_BINDIR)
	cp ../doc/flct.1 $(FLCT_MANDIR)
warpinstall: warp
	cp ./warp $(FLCT_BINDIR)
	cp ../doc/warp.1 $(FLCT_MANDIR)
uninstall: 
	rm -f $(FLCT_BINDIR)/flct
	rm -f $(FLCT_MANDIR)/flct.1
warpuninstall: 
	rm -f $(FLCT_BINDIR)/warp
	rm -f $(FLCT_MANDIR)/warp.1
clean:
	rm -f *.o flct flct.exe warp warp.exe
