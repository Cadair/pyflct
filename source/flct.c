/*

 FLCT: http://solarmuri.ssl.berkeley.edu/overview/publicdownloads/software.html
 Copyright (C) 2007-2009 Regents of the University of California

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
/* # include <complex.h> */
# include <fftw3.h> 
/* 
 * April 10, 2009
 * Version 1.01
 * Eliminated the g1tmp and g2tmp arrays and computed the g1 and g2 arrays
 * directly from f1, f2, f1bar, f2bar, and gaussian mask (gaussdata) arrays.  
 * Got a big speedup (factor of 9 on a mac using case in tests subdirectory).
 *
 * August 12, 2008
 * This version of flct will now be called 1.0 and supercedes all 
 * test_* versions.
 * July 16, 2008 - Have added the ability to interpolate points via cubic
 * convolution at points that are skipped, potentially
 * saving a great deal of computing effort.
 * Makefile now has 'make install', which will install executable and man page.
 * Found and fixed a bug in interpcc2d in which incorrect edge values were
 * assigned.
 *
 * July 10, 2008 - Included the subtraction of a local mean before each
 * subimage is multiplied by a gaussian.  This, in conjunction with filtering
 * (e.g. -k 0.25) results in much better behavior for noisy images.
 *
 * May 31, 2008 - Added ability to skip every N pixels, with p,q as x,y
 * pixel offsets, in response to suggestion by Karin Muglach.  vm mask is
 * set to zero when calculation is skipped, and vx, vy are set to missing value.
 *
 * Fixed bug in which fabs(x) was previously computed as sqrt(x^2).
 * 
 * April 20, 2007 - Version 1.0 - Supercedes version 12 aka test_12
 * Changes from Version test_12 to Version 1.0:
 * Eliminate any C++ style comments to be compatible with older C compilers
 * Change name of executable from vel_ccor to flct
 * Added a Makefile to compile the executable.
 * To compile now, just type "make".  Before doing that,
 * Make sure that LIBFFTW3 is set to the
 * directory containing libfftw3.a, and make sure INCLUDEFFTW3 is set to the
 * directory containing fftw3.h .
 *
 * Add capability of returning single vx, vy values for whole array by setting
 * sigma = 0. (for Deborah Haber's application)
 *
 *
 * May 12, 2006 - Version 12 
 * Changes from version 10 to version 12:
 * We have changed the algorithm for locating the peak of the cc fn to
 * sub-pixel resolution.  In version 10 and before, the maximum of the 
 * function, interpolated to .1 or .02 pixel (hires mode) was returned.
 * This resulted in a "quantization" effect for velocities computed from
 * very small shifts, corresponding to small time steps.  Now, we have
 * adopted the concept of Taylor-expanding the function in the neighborhood
 * of the peak to find the sub-pixel shift values.  The idea for doing this
 * came from examining the "curve fitting" part of Chae's LCT code, but we
 * have implemented this in a somewhat different way that avoids the need
 * to do any explicit linear algebra solves.  
 *
 * The default mode of operation of the code simply uses the un-interpolated
 * cross-correlation function, in the neighborhood of the peak, to estimate
 * the sub-pixel shift.  The "hires" (-h) flag, if present, turns on cubic
 * convolution interpolation at the 0.1 pixel resolution level, and then
 * uses the Taylor-expansion on the interpolated data to find the peak.  This
 * has the result of a smoother derived velocity field, at the possible expense
 * of some numerical accuracy.
 *
 * May 10, 2005
 * Changes from version 8 to version 10:  
 * The binary I/O is now all done in large endian format.  If the code
 * is run on a small endian machine (e.g. linux or windows), 
 * the data is byteswapped before it
 * is read in or written out.  This eliminates the problem of incompatible
 * binary data format between e.g. SUNs or Macs and linux/windows machines.  
 * This should all be transparent to the user, as long as the new companion 
 * IDL i/o routines
 * are used to prepare the input image arrays and to read the output from
 * vel_ccor.  This version also reads/writes a "signature" integer at the front
 * of the file, so that if an input file that is not meant for vel_ccor is
 * mistakenly entered on the input line, it exits gracefully. *NOTE THIS MEANS
 * THE I/O FILES ARE NOT INTERCHANGEABLE BETWEEN VERSIONS test_8 AND test_10*.
 * 
 * The new IDL companion procedures for reading/writing data for vel_ccor
 * are now called vcimage1in.pro, vcimage2in.pro, vcimage3in.pro (input
 * procedures for 1, 2, or 3 2-d arrays), and vcimage1out.pro, vcimage2out.pro, 
 * and vcimage3out.pro (output procedures for 1, 2, or 3 2-d arrays).
 * The syntax for these procedures is the same as for the old versions
 * described below, except now there is no /bs keyword (no need for it).
 * 
 * The threshold value can now be forced into "absolute" mode by appending
 * an 'a' to the numerical threshold value in the calling argument.
 *
 * The calling arguments for vel_ccor ARE changed from version 8.  The 1st
 * five required arguments are the same.  The optional arguments may now
 * appear in any order, and the syntax is -q (for quiet mode), -h (for
 * hires (high resolution) mode, and -t thresh to specify the threshold value
 * for computing the velocity.  The full syntax is now:
 *
 * vel_ccor infile outfile deltat deltas sigma [-q -h -t thresh]
 * 
 * Running vel_ccor with no arguments now results in some terse documentation,
 * in addition to the expected syntax.  
 *
 * Installation instructions identical to those described below.
 *
 * April 21, 2005 - notes on version test_8:
 *
 * This version of vel_ccor (in C) that seems to have
 * most of the functionality and speed we need.  Even though much improved
 * from earlier versions, I am sure much more could be done to speed this up.
 * This version also introduces the use of a "threshold" below which the LCT
 * velocity is not computed.
 *
 * Almost everywhere, the 2-d arrays are stored and accessed as 1-d arrays
 * using pointer arithmetic to convert i,j indices to a single index.  In
 * general the single index has the form ptr+i*ny+j, where i,j are the x,y
 * indices and ptr is the pointer to the beginning of the array.  Here ny
 * is the "y" limit (the 2nd limit) to the array.  Note that in C, arrays
 * are stored backwards (transposed) from the order in IDL and fortran, so that
 * while the images are being processed, the roles of i,j are reversed from 
 * what we'd think
 * of in IDL.  But when the output velocity arrays are read back into IDL,
 * everything comes back in the expected order.
 *
 * To install from source:  First install fftw version 3 (you can get it from
 * www.fftw.org) by downloading the tarball, unpacking it, and compiling the
 * source code.  If you don't need special compilers or other special options,
 * it should build with just 
 *
 * ./configure
 * make 
 * make install (as root).
 *
 * Then compile this (the vel_ccor) source code.  For me, this is just the
 * single line
 * gcc -O3 <name of this file> -lfftw3 -lm -o vel_ccor
 *
 * Depending on how fftw3 was installed, you may need to add an option 
 * -I <directory containing fftw3.h> and/or
 * -L <directory containing libfftw3.a>
 * depending on whether gcc finds these things on its own or not.  You may or
 * may not need to include the -lm option (but I do).
 * 
 * To build a statically compiled version, the compilation command is
 * gcc -static -O3 <name of this file> -lfftw3 -lm -o vel_ccor
 *
 * FFTW3 and this version of vel_ccor build in linux, windows xp, 
 * and solaris.  In windows I used
 * the mingw/msys environment.  The msys version needs to be at least as recent
 * as 1.10.  The procedure to compile vel_ccor in windows and solaris is
 * identical to the procedure for linux described above.
 *
 * vel_ccor currently runs only from the command line.
 *
 * The command line syntax for running vel_ccor is:
 *
 * vel_ccor infile outfile deltat deltas sigma thresh hires quiet
 *
 * infile and outfile are the input and output files; infile will contain
 * 2 images from which LCT flows are determined; outfile will contain the
 * arrays vx, vy, and a velocity mask vm derived from the 2 images.  
 * deltat is the time between
 * the 2 images, deltas is the unit of pixel spacing (ie length), and 
 * sigma is the width * (in pixels) of the gaussian mask used to 
 * modulate the 2 images.  
 *
 * The parameters thresh, hires and quiet are optional, and default values are 
 * to turn these options off (=0).  If thresh is non-zero, then if the
 * average pixel value between the 2 images is less than thresh, no velocity
 * is computed, and the velocity mask value for that pixel
 * is set to 0.  If the value of
 * thresh is between 0 and 1, thresh is assumed to be given in fractional units
 * of the maximum absolute value of the 2 images.  If thresh > 1, then thresh
 * is assumed to be given in absolute units.  If thresh is set to the
 * noise level of the images, significant increases in run time can be
 * achieved, if one is willing to settle for having no computed velocities at
 * those points.
 *
 * If hires is on (!=0) then .02 pixel precision is used
 * in determining flow velocities, otherwise 0.1 pixel precision is assumed.
 * If quiet keyword (!=0) all printing to the terminal or console is suppressed.
 * Turning hires on can substantially increase run time.
 *
 * To write the input file for vel_ccor containing the pair of images, 
 * use the 'images2out' procedure in IDL.
 * To read the output file from vel_ccor containing the vx, vy arrays into IDL, 
 * use the 'images2in' procedure in IDL.  If you also want to read the velocity
 * mask array vm, use the images3in procedure in IDL.
 *
 * Authors:
 * George H. Fisher, Space Sciences Lab # 7450, University of California,
 * 7 Gauss Way, Berkeley CA 94720-7450 email: fisher at ssl dot berkeley dot edu
 *
 * Brian T. Welsch, Space Sciences Lab # 7450, University of California,
 * 7 Gauss Way, Berkeley CA 94720-7450 email: fisher at ssl dot berkeley dot edu
 *
*/

/* global declarations */

/* i4 and f4 are supposed to be definitions that give rise to 4 byte integers
 * and 4 byte floats */


typedef int i4;
typedef float f4;

/* function prototypes: */

int main (int argc, char *argv[]);
i4 read2images (char *fname, i4 * nx, i4 * ny, double **arr, double **barr,
		i4 transp);
i4 where (char *cond, i4 xsize, i4 ** index, i4 * length_index);
i4 cross_cor (i4 init, i4 hires, i4 expand, double *arr, double *barr,
	   double **absccor, i4 nx, i4 ny, double *shiftx, double *shifty, 
           i4 filterflag, double kr);
i4 writeimage (char *fname, double *arr, i4 nx, i4 ny, i4 transp);
i4 write2images (char *fname, double *arr, double *barr, i4 nx, i4 ny,
		 i4 transp);
i4 write3images (char *fname, double *arr, double *barr, double *carr,
		 i4 nx, i4 ny, i4 transp);
i4 shift2d (double *arr, i4 nx, i4 ny, i4 ishift, i4 jshift);
i4 maxloc (double *arr, i4 xsize);
i4 minloc (double *arr, i4 xsize);
i4 iminloc (i4 * arr, i4 xsize);
i4 imaxloc (i4 * arr, i4 xsize);
double r (double t);
i4 interpcc2d (double *fdata, double xmiss, i4 nx, i4 ny, 
double *xwant, i4 nxinterp, double *ywant, i4 nyinterp, double **finterp);
i4 gaussfilt(double *filter, double *kx, double *ky, i4 nx, i4 ny, double kr);
i4 make_freq(double *k, i4 ndim);
i4 filter_image(double *arr, double *barr, double *outarr, double *outbarr,
        i4 nx, i4 ny, double kr);
i4 is_big_endian ();
i4 byteswap (unsigned char *arr, i4 arrsize, i4 nbpw);

/* end function prototypes */

int main (int argc, char *argv[])
{

/* BEGIN MAIN PROGRAM */

  char *version ="1.01   ";
  char infile[100], outfile[100], deltats[100], deltass[100], sigmas[100],
    threshs[100], ks[100],skips[100];
  char *aloc = NULL;
  char *ploc = NULL;
  char *qloc = NULL;
  char *intloc = NULL;
  i4 iarg, quiet, verbose, hires, expand, sigmaeq0, filter;
  i4 nx, nxorig, ny, nyorig, i, j, ii, jj, ic, jc, ixcount=0, iycount=0;
  i4 ier, icc, nsize=0, nt, ndmin, init, ibe;
  i4 iloc1, iloc2, imin0, jmax0, jmin0, imax0, imin, jmin, imax, jmax,
    isize, jsize;
  i4 poffset=0, qoffset=0, skip=0, skipon=0, skipxy,noskipx,noskipy,noskipxy;
  i4 xoffset=0,yoffset=0, interpolate=0;
  double *f1, *f2, *vx, *vy, *vm;
  double *vxnoi=NULL, *vynoi=NULL, *vmnoi=NULL, *xwant, *ywant, *vxint, *vyint;
  double *gaussdata, *f1temp, *f2temp,
    *g1, *g2, *absccor;
  double shiftx, shifty, argx, argy, f1max, f2max, fmax, fabsbar, vmask;
  double deltat, deltas, sigma, thresh, kr=0., sigminv=-999., 
     deltinv, vxx=0., vyy=0., f1bar=0., f2bar=0.;
/*      double tol=1e-4; */
  double tol = 1e-2;
  double xmiss = -999999.;
  i4 hardworkneeded, belowthresh;
  i4 transp = 1; /* This flag is nonzero to transpose input/output arrays */

  /* CODE TO READ IN ARGUMENTS AND OPTIONS: */

  /* check to see if number args is in right range - 
   * if not, print syntax & quit */

  if ((argc < 6) || (argc > 14))
    {
      printf("flct: Version %s Copyright: 2007-2009 University of California\n",
          version);
      printf("Authors: G.H. Fisher, B.T. Welsch, UCB Space Sciences Lab\n\n");
      printf
        ("Syntax: %s ifile ofile deltat deltas sigma -t thr -k kr -s N[pP][qQ][i] -h -q\n\n"
            ,argv[0]);
      printf("ifile - contains 2 images for local correlation tracking\n");
      printf("   (to create ifile use the IDL procedure vcimage2out.pro)\n");
      printf("ofile - contains output vx, vy, vm (x,y velocity and mask\n");
      printf("   arrays).  Mask array 0 where velocity not computed.\n");
      printf("   To read ofile use the IDL procedure vcimage3in.pro)\n");
      printf("deltat - amount of time between images\n");
      printf("deltas - units of length of the side of a single pixel\n");
      printf("   (note velocity is computed in units of deltas/deltat)\n");
      printf("sigma - images modulated by gaussian of width sigma pixels\n");
      printf("(set sigma=0. to just get the overall shift between 2 images)\n");
      printf("Optional parameters:\n");
      printf("-t thr - if avg abs value image pixel < thr, skip calc\n");
      printf("   (if thr between 0 & 1, thr is assumed to be in\n");
      printf("   relative units of the max. abs. pixel value.  To force\n");
      printf("   thr to be in absolute units, append an 'a' to it.)\n");
      printf("-s NpPqQ - skip calc except every N pixels in x and y\n");
      printf("   P is # of offset pixels in x, Q is pixel offset in y\n");
      printf("   e.g. -s 10p5q5 means every 10 pixels, with 5 pixel offsets\n");
      printf("   appending 'i' to string will interpolate at skipped points\n");
      printf("   e.g. -s 5i means skip every 5 points, then interpolate\n");
      printf("-k kr - filter subimages w/gaussian w/ roll-off wavenumber kr\n");
      printf("    - kr is in units of max of kx, ky (typically 0 < kr < 1) \n");
      printf("-h  - 'hires' flag - Use cub. conv. interp. to find shift\n");
      printf("-q  - flag to suppress printing all non-error messages\n");
      exit (1);
    }

  /* GET THE 5 REQUIRED ARGUMENTS */

  /* get input file name */

  strncpy (infile, argv[1], 99);

  /* get output file name */

  strncpy (outfile, argv[2], 99);

  /* get deltat */

  strncpy (deltats, argv[3], 99);
  deltat = atof (deltats);
  deltinv = 1. / deltat;

  /* get deltas */

  strncpy (deltass, argv[4], 99);
  deltas = atof (deltass);

  /* get sigma */

  strncpy (sigmas, argv[5], 99);
  sigma = atof (sigmas);
  if(sigma > 0.) 
  {
     sigminv = 1. / sigma;
     sigmaeq0=0;
  }
  else
  {
     sigmaeq0 = 1;
  }

  /*  GET OPTIONAL ARGUMENTS quiet, hires, expand, thresh, kr, skip : */

  hires = -1;
  quiet = 0;
  expand = 1;
  verbose = 1;
  filter = 0;
  thresh = (double) 0.;

  if(argc >= 7)
  {
    for (iarg=6; iarg < argc; iarg++)
    {
       if(!strncmp("-h",argv[iarg],2))
       {
          hires=0;
       }
       if(!strncmp("-q",argv[iarg],2))
       {
          quiet=1;
       }
       if(!strncmp("-t",argv[iarg],2) && ((iarg+1) < argc))
       {
          /* read in threshhold value */

          strncpy (threshs, argv[iarg+1], 99);

          /* Now check to see if there's an 'a' at end of threshold value
           * string.  If so, threshold value is treated as absolute
           * instead of relative even if it's between 0 and 1. */

          aloc = strchr (threshs, 'a');
          /* if threshold string ends in 'a', aloc will not be NULL.
           * Now, remove the 'a' from the string before doing atof(). */
          if (aloc != NULL) *aloc = '\0';
          thresh = (double) atof (threshs);
       }
       if(!strncmp("-s",argv[iarg],2) && ((iarg+1) < argc))
       {
          /* read in skip information */

          skip=0;
          poffset=0;
          qoffset=0;
          skipon=0;
          interpolate=0;
          strncpy (skips, argv[iarg+1], 99);

          /* Check to see if there's an 'i' in the skip string.
           * If so, set interpolate flag to 1. */
          intloc = strchr(skips, 'i');
          
          /* if skip string has an 'i', intloc will not be NULL */

          if(intloc != NULL)
          {
             interpolate=1;
             *intloc='\0';
          }

          /* Now check to see if there's a 'q' in the skip 
           * string.  If so, find qoffset (stuff after 'q') */

          qloc = strchr (skips, 'q');

          /* if skip string has a 'q', qloc will not be NULL. */

          if (qloc != NULL) 
          {
            qoffset = (i4) atoi (qloc+1);
            *qloc='\0';
          }

          /* Now check to see if there's a 'p' in the skip 
           * string.  If so, find poffset (stuff after 'p') */

          ploc = strchr (skips, 'p');

          /* if skip string has a 'p', ploc will not be NULL. */
          if (ploc != NULL) 
          {
            poffset = (i4) atoi (ploc+1);
            *ploc='\0';
          }

          /* Finally extract skip integer */

          skip= (i4) atoi (skips);
 
          /* Now check that values of poffset,qoffset, and skip make sense */
 
          if( skip <= 0 )
          {
             printf("flct: skip =%d is zero or negative, fatal\n",
               skip);
             exit(1);
          }

          if((abs(poffset) >= skip) || (abs(qoffset) >= skip))

          {
             printf("flct: p=%d,q=%d, abs(p) or abs(q) >= skip=%d\n",
                   poffset,qoffset,skip);
             exit(1);
          }
       }

       skipon=skip+abs(qoffset)+abs(poffset);

       if(!strncmp("-k",argv[iarg],2) && ((iarg+1) < argc))
       {
          /* read in filter roll-off value, relative to max. kx, ky */

          strncpy (ks, argv[iarg+1], 99);
          kr = (double) atof (ks);
          filter=1;
       }

     }
  }

  if (quiet) verbose = 0;

 /* DONE FINDING ARGUMENTS AND OPTIONS */

  /* determine if this is a large endian or small endian platform */
  ibe = is_big_endian ();
  if (verbose)
    {
      printf("flct: Version %s Copyright: 2007,2008 University of California\n",
          version);
      if (ibe)
	{
	  printf ("flct: large endian machine; i/o not byteswapped\n");
	}
      else
	{
	  printf ("flct: small endian machine; i/o will be byteswapped\n");
	}
    }

  /* print out arguments and options */

  if (verbose) printf ("flct: infile = %s\n", infile);
  if (verbose) printf ("flct: outfile = %s\n", outfile);
  if (verbose) printf ("flct: deltat = %g\n", deltat);
  if (verbose) printf ("flct: deltas = %g\n", deltas);
  if (verbose) printf ("flct: sigma = %g\n", sigma);
/* Comment out this warning:
  if (verbose && (sigma < 5.) && (sigma > 0.)) 
       printf ("flct: WARNING: sigma < 5 not recommended\n");
*/
  if (verbose) 
       printf ("flct: threshold image value for LCT is %g\n", thresh);
  if (verbose && aloc) 
       printf ("flct: threshold forced to be in abs. units\n");
  if (verbose && skipon) 
       printf ("flct: skip = %d pixels with p=%d, q=%d\n",skip,poffset,qoffset);
  if (verbose && skipon && ((poffset < 0) || (qoffset < 0))) 
       printf ("flct: p=%d, q=%d: negative p,q will be reset to skip-|p,q|\n",
       poffset,qoffset);
  if (verbose && interpolate) 
       printf ("flct: skipped pixels interpolated with cubic convolution\n");
  if (verbose && filter)
       printf ("flct: filter rolloff value for input images is %g\n", kr);
  if (verbose && (hires == 0)) printf ("flct: hires option turned on\n");

  /* For negative poffset or qoffset, set to skip-poffset or skip-qoffset */
 
  if(poffset < 0) poffset=skip-abs(poffset);
  if(qoffset < 0) qoffset=skip-abs(qoffset);
 
  /* Debug 
  printf("flct: poffset = %d, qoffset = %d\n",poffset,qoffset);
  */
  

  /*
   * read nx, ny, and return references to nx and ny to main prog. *
   * NOTE -- roles of nx, ny are reversed from IDL!!!! In the C version,
   * must work in transposed space.  Therefore transp is set to 1.
   */

  ier = read2images (infile, &nx, &ny, &f1, &f2, transp);

  /* nx, ny may get set to 1 if sigma=0. so copy original values */
  nxorig=nx;
  nyorig=ny;

  if (verbose)
    printf ("flct: from input file, nx = %d, ny = %d\n", nx, ny);
  if ((skip >= nx) || (skip >= ny))
  {
     printf("flct: skip = %d is too big compared to nx or ny, fatal\n",skip);
     exit(1);
  }

  /* Figure out size of sliding box in which FFTs will be done */

  if(!sigmaeq0)
  {
     /* This stuff done only if sigma > 0 */
     nt = (i4) sigma *sqrt (log (1. / tol));
     ndmin = (nx < ny) ? (((nx / 3) / 2)) * 2 : (((ny / 3) / 2)) * 2;
     nsize = ((nt / 2) * 2 < ndmin) ? (nt / 2) * 2 : ndmin;
     if (verbose) printf ("flct: nominal sliding box size = %d\n", 
        2 * nsize);
     if(nsize <= 0)
     {
        printf("flct: error - illegal box size, exiting\n");
        exit(1);
     }
  }
  if(sigmaeq0) 
  {
     /* sigma = 0 means we'll only compute one point */
     nx=1;
     ny=1;
  }

  /* figure out if threshold is in absolute or fractional units
   * and if fractional, convert to absolute units.  If thresh is between
   * zero and 1 (not inclusive) then it's assumed to be fractional,
   * (unless the threshold string ends with 'a') and must be scaled.
   * if aloc == NULL, there's no 'a' in the threshold string. */

  if(!sigmaeq0)
  {
    if ((thresh > 0.) && (thresh < 1.) && (aloc == NULL))
      {
        f1temp = (double *) malloc (sizeof (double) * nx * ny);
        f2temp = (double *) malloc (sizeof (double) * nx * ny);

        for (i = 0; i < nx * ny; i++)

  	{

	  /* compute abs value of f1,f2 arrays as f1temp,
	   * f2temp arrays */

	  *(f1temp + i) = (double) fabs (*(f1 + i));
	  *(f2temp + i) = (double) fabs (*(f2 + i));
	}

      /* now find maximum absolute value of both images */

        iloc1 = maxloc (f1temp, nx * ny);
        iloc2 = maxloc (f2temp, nx * ny);
        f1max = *(f1temp + iloc1);
        f2max = *(f2temp + iloc2);
        fmax = (f1max > f2max) ? f1max : f2max;

      /* now convert relative thresh to absolute threshhold */

        thresh *= fmax;
        if (verbose) 
           printf ("flct: relative threshold in abs. units = %g\n", thresh);

        free (f1temp);
        free (f2temp);
    }
  }

  /* debug: output the two input images to file "f1f2.dat" */
/*
      write2images("f1f2.dat",f1,f2,nx,ny,transp);
*/

  /* Create velocity arrays vx,vy and the velocity mask array vm */

  vx = (double *) malloc (sizeof (double) * nx * ny);
  vy = (double *) malloc (sizeof (double) * nx * ny);

  /* the vm array (velocity mask array, not to be confused with the
   * gaussian mask array that is used to modulate the images) 
   * will later be set to 1.0 where the velocity is computed,
   * and to 0.0 where the velocity is not computed -- because the image
   * pixels are below the noise threshold value "thresh" input from
   * the command line. */

  vm = (double *) malloc (sizeof (double) * nx * ny);

  /* Now create master gaussian image mask: */

  gaussdata = (double *) malloc (sizeof (double) * (2 * nxorig) * (2 * nyorig));

  if(!sigmaeq0) /* this case for sigma > 0 */
  {
    for (i = 0; i < 2 * nxorig; i++)
    {
        argx = sigminv * (double) (i - nxorig);
        for (j = 0; j < 2 * nyorig; j++)
        {
	  argy = sigminv * (double) (j - nyorig);
	  *(gaussdata + i * (2 * ny) + j) = exp (-argx * argx - argy * argy);
        }
    }
  }
  else /* this case for sigma = 0. ie set gaussian to 1.0 */
  {
    for (i = 0; i < 2 * nxorig; i++)
    {
        for (j = 0; j < 2 * nyorig; j++)
        {
	  *(gaussdata + i * (2 * nyorig) + j) = (double) 1.;
        }
    }
  }

      /* Debug -  output the gaussian image mask data: */
      /*
      writeimage("gaussdata.dat",gaussdata,2*nxorig,2*nyorig,transp);
      */

  /* Now do the master loop over i,j for computing velocity field: */

  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
	{
          if((nx ==1) && (ny == 1))
            {
              init=2; /* special case: must initialize AND destroy plans */
            }
	  else if ((i == 0) && (j == 0) && ((nx+ny) > 2))
	    {
	      /* 1st time through, set init to 1 so that
	       * fftw FFT plans are initialized */

	      init = 1;
	    }
	  else if ((i == (nx - 1)) && (j == (ny - 1)) && ((nx+ny) > 2))
	    {
	      /* last time through, set init to -1 so that
	       * fftw static variables are freed and plans
	       * destroyed */

	      init = -1;
	    }
	  else
	    {
	      /* the rest of the time just chunk along */

	      init = 0;
	    }

	  /* Now, figure out if image value is below
	   * threshold: */

	  /* the data is considered below theshhold if the
	   * absolute value of average of the pixel value from the 2 images 
	   * is below threshold */

	  fabsbar = 0.5 * (fabs (*(f1 + i * ny + j) + *(f2 + i * ny + j)));
	  belowthresh = (fabsbar < thresh);

	  /* Or alternatively:
	     belowthresh = ((fabs1 < thresh) || (fabs2 < thresh));
	   */

	  /* all the hard work of doing the cross-correlation
	   * needs to be done if the avg data is above the
	   * threshold OR if init != 0 */

          /* added skip logic here */

          if(skipon)
          {
             if(transp)
             {
                xoffset=qoffset;
                yoffset=poffset;
             }
             else
             {
                xoffset=poffset;
                yoffset=qoffset;
             }
             noskipx = !((i-xoffset) % skip);
             noskipy = !((j-yoffset) % skip);
             noskipxy=noskipx*noskipy;
          }
          else
          {
             noskipxy=1;
          }

	  hardworkneeded = (((!belowthresh) && (noskipxy)) || (init != 0));

	  if (hardworkneeded)
	    {

	      /* the hard work for this particular pixel starts
	       * now */


	      /* Now find where the gaussian modulated image 
	       * is
	       * chopped off by the sliding box.  The sliding
	       * box is centered at i,j, unless the edges of
	       * the box would go outside the array -- 
	       * then the
	       * sliding box just sits at edges and/or corners
	       * of the array   */

              if(!sigmaeq0) /* for sigma > 0 */
              {
	        imin0 = (0 > (i - (nsize - 1))) ? 0 : i - (nsize - 1);
                imax0 = ((nx - 1) < (i + nsize)) ? nx - 1 : i + nsize;
                imin = (imax0 == nx - 1) ? nx - 1 - (2 * nsize - 1) : imin0;
                imax = (imin0 == 0) ? 0 + (2 * nsize - 1) : imax0;

                jmin0 = (0 > (j - (nsize - 1))) ? 0 : j - (nsize - 1);
                jmax0 = ((ny - 1) < (j + nsize)) ? ny - 1 : j + nsize;
                jmin = (jmax0 == ny - 1) ? ny - 1 - (2 * nsize - 1) : jmin0;
                jmax = (jmin0 == 0) ? 0 + (2 * nsize - 1) : jmax0;

                isize = imax - imin + 1;
                jsize = jmax - jmin + 1;

                 /* If the following tests for isize, jsize fail,
	         this is very bad:  exit */

	        if (isize != 2 * nsize)
                {
		  printf ("flct: exiting, bad isize = %d\n", isize);
		  exit (1);
                }
	        if (jsize != 2 * nsize)
                {
		  printf ("flct: exiting, bad jsize = %d\n", jsize);
		  exit (1);
                }
              }
              else /* if sigma = 0. just set isize=nxorig, jsize=nyorig */
              {
                 isize=nxorig;
                 jsize=nyorig;
                 imin=0;
                 jmin=0;
              }
              /* debug:
              printf("isize = %d, jsize = %d,\n",isize,jsize);
              */

              /* Compute sub-image means of f1 and f2: */

              f1bar=0.;
              f2bar=0.;
              for (ii = 0; ii < isize; ii++)
                { 
                   for (jj = 0; jj < jsize; jj++)
                      {
                         f1bar=f1bar+ *(f1 + (ii+imin)*nyorig + (jj+jmin));
                         f2bar=f2bar+ *(f2 + (ii+imin)*nyorig + (jj+jmin));
                      }
                }

              f1bar=f1bar/((double)isize*jsize);
              f2bar=f2bar/((double)isize*jsize);

	      g1 = (double *) malloc (sizeof (double) * isize * jsize);
	      g2 = (double *) malloc (sizeof (double) * isize * jsize);

/* comment out the g1bar calculation:  
   This is not being used, but code retained just in case it is resurrected. */

/*
              g1bar=0.;
              g2bar=0.;
              for (ii = 0; ii < isize; ii++)
                {
                  for (jj = 0; jj < jsize; jj++)
                    {
                       g1bar=g1bar+ *(g1tmp + (ii+imin)*nyorig
                        +(jj+jmin));
                       g2bar=g2bar+ *(g2tmp + (ii+imin)*nyorig
                        +(jj+jmin));
                    }
                }
              g1bar=g1bar/((double)isize*jsize);
              g2bar=g2bar/((double)isize*jsize);
*/

	      /* Now fill the reduced size arrays (sub-images) with the 
	       * appropriate values from the 
	       * full-sized arrays: */

	      for (ii = 0; ii < isize; ii++)
		{
		  for (jj = 0; jj < jsize; jj++)
		    {
		      *(g1 + ii * jsize + jj) = 
                          *(gaussdata + (nxorig-i+(ii+imin))*2*nyorig
                           +nyorig-j+(jj+jmin)) *
                          (*(f1 + (ii + imin) * nyorig + (jj + jmin))-f1bar) ;

		      *(g2 + ii * jsize + jj) = 
                          *(gaussdata + (nxorig-i+(ii+imin))*2*nyorig
                           +nyorig-j+(jj+jmin)) *
                          (*(f2 + (ii + imin) * nyorig + (jj + jmin))-f2bar) ;
		    }
		}

	      /* Call to cross_cor is used to find the 
	       * relative
	       * shift of image g2 relative to g1: */

	      icc = cross_cor (init, hires, expand, g1, g2, &absccor,
			       isize, jsize, &shiftx, &shifty, filter, kr);

/*                              debug:  output of absccor */

/*
                              writeimage("absccor.dat",absccor,isize,jsize,
                                 transp);
*/


	      /* Now free up all the temporary arrays created
	       * during the loop */

	      free (g1);
	      free (g2);
	      free (absccor);

	      /* Now convert shifts to units of velocity using
	       * deltas and deltat */

	      /* Note: if (transp), then the meaning of 
	       * velocities
	       * has to be switched between x and y */

	      if (transp)
		{
		  vxx = shifty * deltinv * deltas;
		  vyy = shiftx * deltinv * deltas;
		}
	      else
		{
		  vxx = shiftx * deltinv * deltas;
		  vyy = shifty * deltinv * deltas;
		}

	      /* all the hard work for this pixel is now done */

	    }

	  /* default value for vmask is 1. */

	  vmask = 1.;

	  if ((belowthresh || !noskipxy) && !sigmaeq0)

	    /* If data below threshold, set vxx, vyy to xmiss 
	     * and vmask to 0, meaning vel not computed. */
            /* If sigma=0 just ignore the threshold and compute anyway */

	    {
	      vxx = xmiss;
	      vyy = xmiss;
	      vmask = 0.;
	    }


	  if ((j == 0) && (verbose))

	    {
	      printf ("flct: progress  i = %d out of %d\r", i, nx - 1);
	      fflush (stdout);
	    }

	  *(vx + i * ny + j) = vxx;
	  *(vy + i * ny + j) = vyy;
	  *(vm + i * ny + j) = vmask;
          if(verbose && sigmaeq0)
          {
              printf("\nflct: vx = %g vy = %g \n",vxx,vyy);
              fflush(stdout);
          }

	}
    }
/* Debug
    printf("\n"); 
*/

  /*  If interpolation to non-computed pixels is desired, this next code is
      where that is done */

  if(skipon && interpolate)
  {
  /* First step is to figure out the number of computed points in x,y 
  which we'll call ixcount and iycount */
  if(transp)
    {
       xoffset=qoffset;
       yoffset=poffset;
    }
  else
    {
       xoffset=poffset;
       yoffset=qoffset;
    }
     ixcount=0;
     for(i=0;i<nx;i++)
     {
        noskipx = !((i-xoffset) % skip);
        if(noskipx) ixcount++;
     }

     iycount=0;
     for(j=0;j<ny;j++)
     {
        noskipy = !((j-yoffset) % skip);
        if(noskipy) iycount++;
     }

  /* Now that we know ixcount,iycount create arrays big enough to hold the
  computed points.  Call these arrays vxnoi, vynoi, plus the mask, vmnoi */

  vxnoi=(void *)malloc(ixcount*iycount*sizeof(double));
  vynoi=(void *)malloc(ixcount*iycount*sizeof(double));
  vmnoi=(void *)malloc(ixcount*iycount*sizeof(double));
  }

  if(skipon && interpolate)
  /* Next step is to fill in the arrays vxnoi, vynoi with the computed
     values from the vy, vy arrays */
  {
  if(transp)
    {
       xoffset=qoffset;
       yoffset=poffset;
    }
  else
    {
       xoffset=poffset;
       yoffset=qoffset;
    }

     ic=-1;
     for(i=0;i<nx;i++)
     {
        noskipx = !((i-xoffset) % skip);
        if(noskipx) ic++;
        jc=-1;
        for(j=0;j<ny;j++)
        {
           noskipy = !((j-yoffset) % skip);
           if(noskipy) jc++;
           noskipxy=noskipx*noskipy; 
           if(noskipxy) *(vxnoi+ic*iycount+jc)=*(vx+i*ny+j);
           if(noskipxy) *(vynoi+ic*iycount+jc)=*(vy+i*ny+j);
           if(noskipxy) *(vmnoi+ic*iycount+jc)=*(vm+i*ny+j);
        }
     }
  }
  /* DEBUG:
  write3images ("noiarrays.dat", vxnoi,vynoi,vmnoi,ixcount,iycount,transp);
  */

  /* Next step is to compute xwant, ywant arrays, to get values of x,y at
     which interpolation is desired.  Note that sometimes these arrays
     will be outside the range of vxnoi, vynoi, especially if non-zero
     values of p,q are used. */

  if(skipon && interpolate)
  {
    xwant=(void *)malloc(nx*sizeof(double));
    ywant=(void *)malloc(ny*sizeof(double));

    for (i=0;i<nx;i++)
    {
    if(transp)
      {
         xoffset=qoffset;
      }
    else
      {
         xoffset=poffset;
      }
      *(xwant+i)=((double)(i-xoffset))/(double)(skip);
/* debug: 
      printf("for i= %d, xwant[i] = %g\n",i,*(xwant+i)); 
*/
    }

    for (j=0;j<ny;j++)
    {
    if(transp)
      {
         yoffset=poffset;
      }
    else
      {
         yoffset=qoffset;
      }
      *(ywant+j)=((double)(j-yoffset))/(double)(skip);
/* debug:
      printf("for j= %d, ywant[j] = %g\n",j,*(ywant+j)); 
*/

    }

/* Next step is to use interpcc2d function to interpolate vxnoi,vynoi arrays
   to vxint,vyint arrays.  Note vxint, vyint are allocated within interpcc2d
   and should be freed when you are done with them. */

      interpcc2d(vxnoi,xmiss,ixcount,iycount,xwant,nx,ywant,ny,&vxint);
      interpcc2d(vynoi,xmiss,ixcount,iycount,xwant,nx,ywant,ny,&vyint);

/* Debug: 
      write3images ("intarrays.dat", vxint,vyint,vyint,nx,ny,transp);
*/

/* Now put in the interpolated values of vxint, vyint back into the vx,vy arrays
   and set the vm array to 0.5 for interpolated values */

      for(i=0;i<nx;i++)
      {
        noskipx = !((i-xoffset) % skip);
        for(j=0;j<ny;j++)
         {
           if(transp)
           {
             xoffset=qoffset;
             yoffset=poffset;
           }
           else
           {
             xoffset=poffset;
             yoffset=qoffset;
           }
           noskipy = !((j-yoffset) % skip);
           noskipxy=noskipx*noskipy; 
           skipxy=!noskipxy;
           if(skipxy) *(vx+i*ny+j)=*(vxint+i*ny+j);
           if(skipxy) *(vy+i*ny+j)=*(vyint+i*ny+j);
           if(skipxy) *(vm+i*ny+j)=0.5;
         }
      }
     /* Now need to free all the arrays that have been created for
        the interpolation operation */
    free(vxnoi);
    free(vynoi);
    free(vmnoi);
    free(xwant);
    free(ywant);
    free(vxint);
    free(vyint);
    /* Should finally be done with the interpolation over skipped data */
  }

  /* Finally, reset any values of vx,vy that are equal to xmiss to zero, and
     make sure corresponding value of vm is also zero. */

  for(i=0;i<nx;i++)
   {
      for(j=0;j<ny;j++)
       {
          if(*(vx+i*ny+j) == xmiss) *(vm+i*ny+j)=0.;
          if(*(vx+i*ny+j) == xmiss) *(vx+i*ny+j)=0.;
          if(*(vy+i*ny+j) == xmiss) *(vy+i*ny+j)=0.;
       }
   }

  /* Outer loops over i,j finally done! */
  /* Output the vx, vy arrays to the output file 'outfile': */

  write3images (outfile, vx, vy, vm, nx, ny, transp);

  /* free the gaussian mask array, the original images, and the
   * velocity arrays */

  free (gaussdata);
  free (f1);
  free (f2);
  free (vx);
  free (vy);
  free (vm);

  if (verbose)
    printf ("\nflct: finished\n");

  /* we're done! */
  return 0;
  /*  END MAIN PROGAM */
}

i4 read2images (char *fname, i4 * nx, i4 * ny, double **arr, double **barr,
	     i4 transp)
/* Function to read array dims, create space for 2 double arrays, read them in.
 * Note the double pointer to the double precision * arrays in the calling 
 * argument.  Note also these are referenced as pointers and returned to the
 * calling function. */
{
  FILE *f1;			/* pointer to input file */
  i4 newsize;			/* size of the new double prec. array to be 
				   read in =nx*ny */
  i4 i, ier, ibe, ise, vcid;
  f4 *farr, *fbarr;
  ibe = is_big_endian ();
  ise = 0;
  if (ibe == 0) ise = 1;	/* set small endian flag if not big  endian */
  f1 = fopen (fname, "rb");	/* open file for binary unformatted read */
  ier = 0;
  if (f1 == NULL)		/* error exit if file can't be opened */
    {
      printf ("read2images: cannot open file %s\n", fname);
      exit (1);
    }

  /* check that files begins with special vel_ccor flag: */
  fread (&vcid, sizeof (i4), 1, f1);
  if (ise) byteswap ((void *) &vcid, 1, sizeof (i4));
  if (vcid != 2136967593)
    {
      printf ("read2images: input file is not a vel_ccor i/o file\n");
      exit (1);
    }

  /* order of nx, ny read in must be switched if transp is true */

  if (transp)
    {
      fread (ny, sizeof (i4), 1, f1);
      fread (nx, sizeof (i4), 1, f1);
    }
  else
    {
      fread (nx, sizeof (i4), 1, f1);
      fread (ny, sizeof (i4), 1, f1);
    }
  if (ise)			/* byteswap nx,ny if small endian */
    {
      byteswap ((void *) nx, 1, sizeof (i4));
      byteswap ((void *) ny, 1, sizeof (i4));
    }
/*
	printf("\n\nnx,ny read in from file arr = %d,%d\n",*nx,*ny);
*/
  newsize = (*nx) * (*ny) * sizeof (double);	/* size of new double array */

  /* now create enough space in memory to hold the array arr */

  *arr = malloc (newsize);

  /* allocate space for the temporary, f4 array farr */

  farr = malloc (sizeof (f4) * (*nx) * (*ny));

  if (!*arr)
    {				/* check for error in memory allocation */
      printf ("read2images: memory request for arr failed\n");
      exit (1);
    }
/*
	printf("%d bytes of memory now allocated for arr \n",newsize);
*/

  /* now read in the arr array */

  fread (farr, sizeof (f4), (*nx) * (*ny), f1);
  if (ise) byteswap ((void *) farr, (*nx) * (*ny), sizeof (f4));
  /*byteswap if needed */

  /* now create space for temp. f4 array fbarr, and returned
   * array barr: */

  fbarr = malloc (sizeof (f4) * (*nx) * (*ny));
  *barr = malloc (newsize);

  if (!*barr)
    {				/* check for error in memory allocation */
      printf ("read2images: memory request for barr failed\n");
      exit (1);
    }

  /* now read in the fbarr array */

  fread (fbarr, sizeof (f4), (*nx) * (*ny), f1);
  /*byteswap if needed */
  if (ise) byteswap ((void *) fbarr, (*nx) * (*ny), sizeof (f4));

  /* now transfer data from temp. arrays to arr and barr: */

  for (i = 0; i < (*nx) * (*ny); i++)
    {
      *(*arr + i) = (double) *(farr + i);
      *(*barr + i) = (double) *(fbarr + i);
    }

  /* now free the temp. arrays and close the files */

  free (farr);
  free (fbarr);
  fclose (f1);
  ier = 1;
  return ier;
}

i4 where (char *cond, i4 xsize, i4 ** index, i4 * length_index)
/* This function serves a similar purpose to the "where" function in IDL.
 * Given the array *cond (char) of length xsize 
 * containing a pre-computed condition, the function finds **index, a double
 * pointer to an array of indices which reference those values
 * where *cond is
 * non-zero.  The integer *length_index returns the length of *index. 
 */
/* Note - this function no longer used in vel_ccor */
{
  i4 ier;	/* function return value - not thought through yet */
  i4 i, ii;	/* counter variables */
  i4 *indtmp;	/* temporary local array of indices of *x */
  ier = 0;	/*return value of function */
  *length_index = 0;	/* initialize length of *index array to 0 */

/*	printf("\nxsize = %d",xsize); */

  indtmp = (i4 *) malloc (sizeof (i4) * xsize);	/* create temp. ind. array */

  /* Ready to start */

  ii = 0;		/* set initial counter of temp index array to 0 */
  for (i = 0; i < xsize; i++)	/* start incrementing the *cond array: */
    {
      if (*(cond + i))
	{
	  /* if the condition is true, record i into temp index */
	  *(indtmp + ii) = (i4) i;
	  ii++;		/* and then increment ii */
	}
      /* otherwise just keep incrementing i and doing nothing */
    }

/*	printf("\nii= %d\n", ii) ;
	fflush (stdout) ; */

  /* Now create right amount of space for *index: */

  *index = (i4 *) malloc (sizeof (i4) * ii);

  /* Now copy index values from temp array into *index array */

  memcpy ((void *) *index, (void *) indtmp, ii * sizeof (i4));

  /* Now set the length of the *index array */

  *length_index = (i4) ii;

  /* Now free memory from temp. index array */

  free (indtmp);
  return ier;			/* always 0 at the moment */
}

i4 cross_cor (i4 init, i4 hires, i4 expand, double *arr, double *barr,
	   double **absccor, i4 nx, i4 ny, double *shiftx, double *shifty, 
           i4 filterflag, double kr)
{
/*  To use C99 complex arithmetic, uncomment next line */
/* #include <complex.h> */

/* #include <fftw3.h> */
  i4 i, j, ixx, iyy, maxind, ixmax, iymax, ishft, maxfine, absccmax;
  i4 nxinterp, nyinterp, nfgppergp;
  double normfac, rangex, rangey, shiftx0, shifty0, shiftxx, shiftyy;
  double *xwant, *ywant, *peakarea;
  double shiftsubx, shiftsuby, fx, fy, fxx, fyy, fxy;
  double xmiss=0.;

  /* following variables must be saved between calls; declared static */

  static double *ina, *inb, *ccor;
  static double *filter, *kx, *ky;
  static fftw_complex *outa, *outb, *ccorconj;
  static fftw_plan pa, pb, pback;

  /* absccor is a double pointer containing abs. value of cc function */

  /* debug:
  printf("cross_cor: nx = %d, ny = %d\n",nx,ny);
  */

  *absccor = malloc (sizeof (double) * nx * ny);

  /* Set up interpolation grid depending on whether hires set or not */

  if (hires == 1)
    {
      nxinterp = 101;
      nyinterp = 101;
      nfgppergp = 50;
    }
  else
    {
      nxinterp = 21;
      nyinterp = 21;
      nfgppergp = 10;
    }
/*	printf("initialization stuff done in cross_cor\n"); */
  if ((init == 1) || (init == 2))
    {
      /* First time through: */
      /* Initialization of FFT variables and FFTW plans. */
      /* NOTE -- empirically had to add 1 to "y" dimensions of outa,
       * outb, and ccorconj to
       * avoid a memory leak and a seg fault at fftw_free */

      /* should check to see if still a problem */

      outa = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
					   nx * ((ny / 2) + 2));
      outb = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * 
              nx * ((ny / 2) + 2));	/* valgrind sometimes complains */
      ccorconj = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
					       nx * ((ny / 2) + 2));

      ina = (double *) fftw_malloc (sizeof (double) * nx * ny);
      inb = (double *) fftw_malloc (sizeof (double) * nx * ny);
      ccor = (double *) fftw_malloc (sizeof (double) * nx * ny);
      filter = (double *) fftw_malloc (sizeof (double) * nx * ny);
      kx=(double *) fftw_malloc (sizeof(double)*nx);
      ky=(double *) fftw_malloc (sizeof(double)*ny);
      if(filterflag)
      {
         make_freq(kx,nx);
         make_freq(ky,ny);
         gaussfilt(filter,kx,ky,nx,ny,kr);
      }

      for (i = 0; i < nx * ny; i++)
	{
	  *(ina + i) = (double) 0.;
	  *(inb + i) = (double) 0.;
	}
      for (i = 0; i < nx * ((ny / 2 + 1)); i++);
      {
/*                      *(ccorconj+i)=0.; */

	/* above line works if complex.h is invoked, but we instead
	 * follow the advice on complex arithmetic described in
	 * the fftw3 tutorial */

	ccorconj[i][0] = 0.;
	ccorconj[i][1] = 0.;
      }
      pa = fftw_plan_dft_r2c_2d (nx, ny, ina, outa, FFTW_MEASURE);
      pb = fftw_plan_dft_r2c_2d (nx, ny, inb, outb, FFTW_MEASURE);
      pback = fftw_plan_dft_c2r_2d (nx, ny, ccorconj, ccor, FFTW_MEASURE);
    }

    for (i = 0; i < nx * ny; i++)

    {
/*		printf("1st loop: i = %d, *(arr+i)= %g, *(barr+i) = %g\n",
				i,*(arr+i),*(barr+i)); */

      /* copy from input doubles to fftw variables */

      *(ina + i) = (double) (*(arr + i));
      *(inb + i) = (double) (*(barr + i));
    }

  /* actually do the forward FFTs: */

  fftw_execute (pa);
  fftw_execute (pb);

  /* calculate normalization factor */

  normfac = (1. / ((double) nx * ny));
  normfac *= normfac; /* square of above line */
 
/* Now apply the gaussian filter to the FFT'd data in frequency space */

    if(filterflag)
    {
      for (i=0;i<nx;i++)
      {
         for (j=0;j<(ny/2)+1;j++)
         {
            outa[i*((ny/2)+1)+j][0]=outa[i*((ny/2)+1)+j][0]*filter[i*ny+j];
            outa[i*((ny/2)+1)+j][1]=outa[i*((ny/2)+1)+j][1]*filter[i*ny+j];
            outb[i*((ny/2)+1)+j][0]=outb[i*((ny/2)+1)+j][0]*filter[i*ny+j];
            outb[i*((ny/2)+1)+j][1]=outb[i*((ny/2)+1)+j][1]*filter[i*ny+j];
         }
      }
    }

  /* Now calculate product of conj(outa) * outb */

  for (i = 0; i < nx * ((ny/2) + 1); i++)
    {

/*          *(ccorconj+i)=(conj(*(outa+i))*(*(outb+i)))*normfac; */

      /* the above commented out line can be used if complex.h
       * is invoked, but this is only accepted by C99 compliant
       * compilers.  For now, we used the suggested approach of
       * of the fftw tutorial.  This appears to be required if
       * one uses old e.g. 2.95 or 2.96 versions of gcc or 
       * mingw/msys in windows systems */

      ccorconj[i][0] = (outa[i][0] * outb[i][0] + outa[i][1] * outb[i][1])
	* normfac;
      ccorconj[i][1] = (outa[i][0] * outb[i][1] - outa[i][1] * outb[i][0])
	* normfac;
    }

  /* now do the inverse transform to get cc function */

  fftw_execute (pback);

  /* now calculate the absolute value of cc function */

  for (i = 0; i < nx * ny; i++)
    {
      *(*absccor + i) = (double) fabs(*(ccor+i));
    }

  if ((init == -1) || (init == 2))
    {
      /* Last time through: free all the plans and static variables */

      fftw_free (outa);
      fftw_free (outb);
      fftw_free (ccorconj);
      fftw_free (ccor);
      fftw_free (filter);
      fftw_free (kx);
      fftw_free (ky);
      fftw_free (ina);
      fftw_free (inb);
      fftw_destroy_plan (pa);
      fftw_destroy_plan (pback);
      fftw_destroy_plan (pb);
    }

/* Now shift the absccor array by nx/2, ny/2 to avoid aliasing problems */

  ishft = shift2d (*absccor, nx, ny, nx / 2, ny / 2);

  /* Now find maximum of the shifted cross-correlation function to 1 pixel
     accuracy:  */

  absccmax=1;
  maxind = maxloc (*absccor, nx * ny);
  if( *(*absccor+maxind) == (double)0.) 
  {
     absccmax=0;
  }
  if(absccmax == 1)
  {
     ixmax = maxind / ny;
     iymax = maxind % ny;
  }
  else
  {
     ixmax = nx/2;
     iymax = ny/2;
  }
  shiftx0 = ixmax;
  shifty0 = iymax;
  shiftsubx=0.;
  shiftsuby=0.;

  if(( expand == 1) && (hires == -1) && (ixmax > 0) && (ixmax < (nx-1))
     && (iymax > 0) && (iymax < (ny-1)) && (absccmax == 1))
  {
     fx=0.5* ( *(*absccor+(ixmax+1)*ny+iymax) - 
         *(*absccor+(ixmax-1)*ny+iymax) );
     fy=0.5* ( *(*absccor+ixmax*ny+iymax+1) - *(*absccor+ixmax*ny+iymax-1) );
     fxx = ( *(*absccor+(ixmax+1)*ny+iymax)+ *(*absccor+(ixmax-1)*ny+iymax)
        -2.*( *(*absccor+ixmax*ny+iymax))  );
     fyy = ( *(*absccor+ixmax*ny+iymax+1) + *(*absccor+ixmax*ny+iymax-1)
        -2.*( *(*absccor+ixmax*ny+iymax)) );
     fxy = 0.25*( *(*absccor+(ixmax+1)*ny+iymax+1) + 
            *(*absccor+(ixmax-1)*ny+iymax-1) -
            *(*absccor+(ixmax+1)*ny+iymax-1) - 
            *(*absccor+(ixmax-1)*ny+iymax+1) );
/* In following expressions for subshifts, shift is in units of pixel length */
     shiftsubx=(fyy*fx-fy*fxy)/(fxy*fxy-fxx*fyy);
     shiftsuby=(fxx*fy-fx*fxy)/(fxy*fxy-fxx*fyy);
  }

  shiftxx=shiftx0 + shiftsubx;
  shiftyy=shifty0 + shiftsuby;
/*
       printf("shiftx0-nx/2 = %g\n",(shiftx0-(double)(nx/2)));
       printf("shifty0-ny/2 = %g\n",(shifty0-(double)(ny/2)));
 
*/

/* Now create x, y arrays to define desired interpolation grid: */
  if(hires != -1)
  {

     rangex = (double) (nxinterp - 1) / nfgppergp;
     rangey = (double) (nyinterp - 1) / nfgppergp;

     xwant = (double *) malloc (sizeof (double) * nxinterp);
     ywant = (double *) malloc (sizeof (double) * nyinterp);

     for (i = 0; i < nxinterp; i++)
       {
         *(xwant + i) = ((((double) i) * rangex) / ((double) (nxinterp - 1)))
	   - 0.5 * rangex + shiftx0;
/*                 printf("xwant[%d] = %g\n",i,*(xwant+i)); */
       }
     for (j = 0; j < nyinterp; j++)
       {
         *(ywant + j) = ((((double) j) * rangey) / ((double) (nyinterp - 1)))
   	- 0.5 * rangey + shifty0;
/*                 printf("ywant[%d] = %g\n",j,*(ywant+j)); */
       }
   
  /* Now, do the interpolation of the region of the peak of the cc fn */

     interpcc2d (*absccor, xmiss, nx, ny, xwant, nxinterp, ywant, 
          nyinterp, &peakarea);

  /* Following writeimage stmt is available if you need to examine the
   * peakarea array for debugging - note transpose of x,y for IDL  read */

/*
      transp=1;
      writeimage("peakarea.dat",peakarea,nxinterp,nyinterp,transp);
*/

  /* Now find the peak of the interpolated function */

     maxfine = maxloc (peakarea, nxinterp * nyinterp);
     ixx = maxfine / nyinterp;
     iyy = maxfine % nyinterp;
/* Here is where to compute sub-pixel shifts in peakarea if they're wanted */
     shiftsubx=0.;
     shiftsuby=0.;
     if((expand == 1) && (ixx > 0) && (ixx < (nxinterp-1)) && (iyy > 0)
        && (iyy < (nyinterp-1)))
     {
        fx=0.5* ( *(peakarea+(ixx+1)*nyinterp+iyy) - 
              *(peakarea+(ixx-1)*nyinterp+iyy) );
        fy=0.5* ( *(peakarea+ixx*nyinterp+iyy+1) - 
              *(peakarea+ixx*nyinterp+iyy-1) );
        fxx = ( *(peakarea+(ixx+1)*nyinterp+iyy)+ 
              *(peakarea+(ixx-1)*nyinterp+iyy)
              -2.*( *(peakarea+ixx*nyinterp+iyy))  );
        fyy = ( *(peakarea+ixx*nyinterp+iyy+1) + *(peakarea+ixx*nyinterp+iyy-1)
           -2.*( *(peakarea+ixx*nyinterp+iyy)) );
        fxy = 0.25*( *(peakarea+(ixx+1)*nyinterp+iyy+1) +
               *(peakarea+(ixx-1)*nyinterp+iyy-1) -
               *(peakarea+(ixx+1)*nyinterp+iyy-1) -
               *(peakarea+(ixx-1)*nyinterp+iyy+1) );
/* In following expressions for subshifts, must mpy by unit of pixel length */
        shiftsubx=((fyy*fx-fy*fxy)/(fxy*fxy-fxx*fyy))*
          rangex/((double) (nxinterp -1));
        shiftsuby=((fxx*fy-fx*fxy)/(fxy*fxy-fxx*fyy))*
          rangey/((double) (nyinterp -1));
     }
     shiftxx = *(xwant + ixx) + shiftsubx;
     shiftyy = *(ywant + iyy) + shiftsuby;
  /* Free the variables created during interpolation */
     free (xwant);
     free (ywant);
     free (peakarea);
  }

/* Now, assign values to shiftx, shifty to return to calling program */

  *shiftx = shiftxx - (double) (nx / 2);
  *shifty = shiftyy - (double) (ny / 2);

/* Following expressions used if only 1 pixel accuracy needed from absccor 
 *
	*shiftx=((double)ixmax)-(double)(nx/2);
	*shifty=((double)iymax)-(double)(ny/2);
*/

  return 0;
}

i4 make_freq(double *k, i4 ndim)
{
/* k is assumed already allocated in main program, with dimension ndim */
i4 n21,i,inext;
n21=(ndim/2)-1;
for (i=0;i<n21+1;i++)
{
	k[i]=(double)i;
}

inext=n21+1;
if((ndim/2)*2 != ndim)
{
	k[inext]=(double)(ndim/2);
	inext++;
        k[inext]=-(double)(ndim/2);
        inext++;
}

else
{
	k[inext]=(double)(ndim/2);
        inext++;
}

for (i=inext;i<ndim;i++)
{
	k[i]=-(double)(n21-(i-inext));
}
/* debug */

/*
for (i=0;i<ndim;i++)
{
	printf("i = %d, k = %g\n",i,k[i]);
}
*/

/* end debug */

return 0;
}

i4 gaussfilt(double *filter, double *kx, double *ky, i4 nx, i4 ny, double kr)
{
/* Assumes kx of size nx, ky of size ny, and filter of size (nx,ny) */
double kxmax,kymax,kxroll,kyroll,smxinv,smyinv,argx,argy;
i4 i,j;
kxmax=(double)kx[nx/2];
kymax=(double)ky[ny/2];
kxroll=kr*kxmax;
kyroll=kr*kymax;
smxinv=(double)1./kxroll;
smyinv=(double)1./kyroll;
for (i=0;i<nx;i++)
{
	argx=kx[i]*smxinv;
	for(j=0;j<ny;j++)
	{
                argy=ky[j]*smyinv;
		filter[i*ny+j]=exp( -(argx*argx + argy*argy) );
	}
}
return 0;
}

i4 filter_image(double *arr, double *barr, double *outarr, double *outbarr,
        i4 nx, i4 ny, double kr)

/* Takes images arr, barr and filters them by gaussians in k-space of width
kr*kmax, where 0 < kr < 1, and kmax is the maximum wavenumber in x and y,
considered separately.  The input arrays are arr, barr, and the output
arrays are outarr, and outbarr.  They are assumed already allocated in
caller.  

This function is not used in this particular version of flct (filtering
is done within cross_cor), but is
included as it may be useful in the future.

*/

{

/*  To use C99 complex arithmetic, uncomment next line */
/* #include <complex.h> */

/* #include <fftw3.h> */

  i4 i,j;
  double *ina, *inb;
  double *filter, *kx, *ky;
  double normfac;
  fftw_complex *outa, *outb;
  fftw_plan pa, pb, pbacka, pbackb;
  outa = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
     nx * ((ny / 2) + 2));
  outb = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * 
     nx * ((ny / 2) + 2));	/* valgrind sometimes complains */
  ina = (double *) fftw_malloc (sizeof (double) * nx * ny);
  inb = (double *) fftw_malloc (sizeof (double) * nx * ny);
  filter = (double *) fftw_malloc (sizeof (double) * nx * ny);
  kx=(double *) fftw_malloc (sizeof(double)*nx);
  ky=(double *) fftw_malloc (sizeof(double)*ny);
  make_freq(kx,nx);
  make_freq(ky,ny);
  gaussfilt(filter,kx,ky,nx,ny,kr);
      for (i = 0; i < nx * ny; i++)
	{
	  *(ina + i) = (double) 0.;
	  *(inb + i) = (double) 0.;
	}
      for (i = 0; i < nx * ((ny / 2 + 1)); i++);
      {
/*                      *(outa+i)=0.; */
/*                      *(outb+i)=0.; */

	/* above lines work if complex.h is invoked, but we instead
	 * follow the advice on complex arithmetic described in
	 * the fftw3 tutorial */

	outa[i][0] = 0.;
	outa[i][1] = 0.;
	outb[i][0] = 0.;
	outb[i][1] = 0.;
      }
      /* set up plans for FFTs: */
      pa = fftw_plan_dft_r2c_2d (nx, ny, ina, outa, FFTW_MEASURE);
      pb = fftw_plan_dft_r2c_2d (nx, ny, inb, outb, FFTW_MEASURE);
      pbacka = fftw_plan_dft_c2r_2d (nx, ny, outa, ina, FFTW_MEASURE);
      pbackb = fftw_plan_dft_c2r_2d (nx, ny, outb, inb, FFTW_MEASURE);

    for (i = 0; i < nx * ny; i++)

    {
/*		printf("1st loop: i = %d, *(arr+i)= %g, *(barr+i) = %g\n",
				i,*(arr+i),*(barr+i)); */

      /* copy from input doubles to fftw variables */

      *(ina + i) = (double) (*(arr + i));
      *(inb + i) = (double) (*(barr + i));
    }

  /* actually do the forward FFTs: */

  fftw_execute (pa);
  fftw_execute (pb);
 /* calculate normalization factor */

  normfac = (1. / ((double) nx * ny));

/* Now apply the gaussian filter to the FFT'd data in frequency space */

    for (i=0;i<nx;i++)
    {
      for (j=0;j<(ny/2)+1;j++)
      {
         outa[i*((ny/2)+1)+j][0]=outa[i*((ny/2)+1)+j][0]*filter[i*ny+j]*normfac;
         outa[i*((ny/2)+1)+j][1]=outa[i*((ny/2)+1)+j][1]*filter[i*ny+j]*normfac;
         outb[i*((ny/2)+1)+j][0]=outb[i*((ny/2)+1)+j][0]*filter[i*ny+j]*normfac;
         outb[i*((ny/2)+1)+j][1]=outb[i*((ny/2)+1)+j][1]*filter[i*ny+j]*normfac;
      }
    }

/* now do the inverse transform to get filtered images: */

    fftw_execute (pbacka);
    fftw_execute (pbackb);

  for (i = 0; i < nx * ny; i++)
    {
      *(outarr+i)=(double) (*(ina+i));
      *(outbarr+i)=(double) (*(inb+i));
    }

/* Free the plans and locally created arrays */

      fftw_free (outa);
      fftw_free (outb);
      fftw_free (filter);
      fftw_free (kx);
      fftw_free (ky);
      fftw_free (ina);
      fftw_free (inb);
      fftw_destroy_plan (pa);
      fftw_destroy_plan (pbacka);
      fftw_destroy_plan (pb);
      fftw_destroy_plan (pbackb);

return 0;

}

i4 writeimage (char *fname, double *arr, i4 nx, i4 ny, i4 transp)
{

/* Function to write array dimensions, and then write out array */

  FILE *f1;
  i4 newsize;
  i4 i, ier, ibe, ise, vcid, vcidtmp;
  f4 *farr;
  i4 nxtmp, nytmp;
  vcid = 2136967593;
  vcidtmp = vcid;
  nxtmp = nx;
  nytmp = ny;
  ibe = is_big_endian ();
  ise = 0;
  if (ibe == 0)
    ise = 1;	/* set flag for byteswapping if small endian */

  /* get the file fname open for binary write */

  f1 = fopen (fname, "wb");
  ier = 0;
  if (f1 == NULL)
    {
      printf ("writeimage: cannot open file %s\n", fname);
      exit (1);
    }

  if (ise)
    {
      byteswap ((void *) &vcidtmp, 1, sizeof (i4));
      byteswap ((void *) &nxtmp, 1, sizeof (i4));
      byteswap ((void *) &nytmp, 1, sizeof (i4));
    }

  fwrite (&vcidtmp, sizeof (i4), 1, f1);  /* write vel_ccor id integer 1st */

  /* shift role of nx, ny if transp is set (for use of data w/ IDL */

  if (transp)
    {
      fwrite (&nytmp, sizeof (i4), 1, f1);
      fwrite (&nxtmp, sizeof (i4), 1, f1);
    }
  else
    {
      fwrite (&nxtmp, sizeof (i4), 1, f1);
      fwrite (&nytmp, sizeof (i4), 1, f1);
    }


/*
      printf("\n\nnx,ny wrote out to file arr = %d,%d\n",nx,ny);
*/
  newsize = nx * ny * sizeof (double);

  /* create temporary f4 array to write out */

  farr = (f4 *) malloc (sizeof (f4) * nx * ny);

  /* now fill the array */

  for (i = 0; i < nx * ny; i++)
    {
      *(farr + i) = (f4) * (arr + i);
    }

  if (ise)
    byteswap ((void *) farr, nx * ny, sizeof (f4));
  /* byteswap if small endian */

  /* write out the array */
  fwrite (farr, sizeof (f4), nx * ny, f1);

  /* free temp. array and close file */
  free (farr);
  fclose (f1);
  ier = 1;
  return ier;
}

i4 write2images (char *fname, double *arr, double *barr, 
         i4 nx, i4 ny, i4 transp)
{

/* Function to write array dimensions, and write out 2 arrays arr and barr
 * while converting them to single precision f4 arrays  */

  FILE *f1;
  i4 newsize;
  i4 ier, i, ise, ibe;
  i4 nxtmp, nytmp, vcid, vcidtmp;
  f4 *farr, *fbarr;
  nxtmp = nx;
  nytmp = ny;
  vcid = 2136967593;
  vcidtmp = vcid;
  ibe = is_big_endian ();
  ise = 0;
  if (ibe == 0) ise = 1;
  if (ise)		/* if small endian, byteswap nxtmp and nytmp */
    {
      byteswap ((void *) &vcidtmp, 1, sizeof (i4));
      byteswap ((void *) &nxtmp, 1, sizeof (i4));
      byteswap ((void *) &nytmp, 1, sizeof (i4));
    }

  /* open the file fname for a binary write */

  f1 = fopen (fname, "wb");
  ier = 0;
  if (f1 == NULL)
    {
      printf ("write2images: cannot open file %s\n", fname);
      exit (1);
    }

  fwrite (&vcidtmp, sizeof (i4), 1, f1);	/* write vel_ccor id flag */

  /* switch role of nx, ny if transp is set (for interaction with IDL) 
   * and write these into the file */

  if (transp)
    {
      fwrite (&nytmp, sizeof (i4), 1, f1);
      fwrite (&nxtmp, sizeof (i4), 1, f1);
    }
  else
    {
      fwrite (&nxtmp, sizeof (i4), 1, f1);
      fwrite (&nytmp, sizeof (i4), 1, f1);
    }

/*
      printf("\n\nnx,ny wrote out to file arr = %d,%d\n",nx,ny);
*/
  newsize = nx * ny * sizeof (double);

  /* create space for the temporary f4 arrays farr and fbarr */

  farr = (f4 *) malloc (sizeof (f4) * nx * ny);
  fbarr = (f4 *) malloc (sizeof (f4) * nx * ny);

  /* fill the temporary arrays */

  for (i = 0; i < nx * ny; i++)
    {
      *(farr + i) = (f4) * (arr + i);
      *(fbarr + i) = (f4) * (barr + i);
    }

  /* now write out the 2 arrays */

  if (ise)	/* byteswap if small endian */
    {
      byteswap ((void *) farr, nx * ny, sizeof (f4));
      byteswap ((void *) fbarr, nx * ny, sizeof (f4));
    }

  fwrite (farr, sizeof (f4), nx * ny, f1);
  fwrite (fbarr, sizeof (f4), nx * ny, f1);

  /* free temp arrays and close file */

  free (farr);
  free (fbarr);
  fclose (f1);
  ier = 1;
  return ier;
}

i4 write3images (char *fname, double *arr, double *barr, double *carr,
	      i4 nx, i4 ny, i4 transp)
{

/* Function to write array dimensions, and write out 3 arrays arr,barr, and carr
 * while converting them to single precision f4 arrays  */

  FILE *f1;
  i4 newsize;
  i4 ier, i, ibe, ise;
  i4 nxtmp, nytmp, vcid, vcidtmp;
  f4 *farr, *fbarr, *fcarr;
  nxtmp = nx;
  nytmp = ny;
  vcid = 2136967593;
  vcidtmp = vcid;
  ibe = is_big_endian ();
  ise = 0;
  if (ibe == 0) ise = 1;	/* test for small endian for doing byteswaps */
  if (ise)			/* byteswap nxtmp, nytmp if small endian */
    {
      byteswap ((void *) &vcidtmp, 1, sizeof (i4));
      byteswap ((void *) &nxtmp, 1, sizeof (i4));
      byteswap ((void *) &nytmp, 1, sizeof (i4));
    }

  /* open the file fname for a binary write */

  f1 = fopen (fname, "wb");
  ier = 0;
  if (f1 == NULL)
    {
      printf ("write3images: cannot open file %s\n", fname);
      exit (1);
    }

  fwrite (&vcidtmp, sizeof (i4), 1, f1);

  /* switch role of nx, ny if transp is set (for interaction with IDL) 
   * and write these into the file */

  if (transp)
    {
      fwrite (&nytmp, sizeof (i4), 1, f1);
      fwrite (&nxtmp, sizeof (i4), 1, f1);
    }
  else
    {
      fwrite (&nxtmp, sizeof (i4), 1, f1);
      fwrite (&nytmp, sizeof (i4), 1, f1);
    }

/*
      printf("\n\nnx,ny wrote out to file arr = %d,%d\n",nx,ny);
*/
  newsize = nx * ny * sizeof (double);

  /* create space for the temporary f4 arrays farr, fbarr, and fcarr */

  farr = (f4 *) malloc (sizeof (f4) * nx * ny);
  fbarr = (f4 *) malloc (sizeof (f4) * nx * ny);
  fcarr = (f4 *) malloc (sizeof (f4) * nx * ny);

  /* fill the temporary arrays */

  for (i = 0; i < nx * ny; i++)
    {
      *(farr + i) = (f4) * (arr + i);
      *(fbarr + i) = (f4) * (barr + i);
      *(fcarr + i) = (f4) * (carr + i);
    }

  if (ise)			/* if small endian, byteswap the arrays */
    {
      byteswap ((void *) farr, nx * ny, sizeof (f4));
      byteswap ((void *) fbarr, nx * ny, sizeof (f4));
      byteswap ((void *) fcarr, nx * ny, sizeof (f4));
    }

  /* now write out the 3 arrays */

  fwrite (farr, sizeof (f4), nx * ny, f1);
  fwrite (fbarr, sizeof (f4), nx * ny, f1);
  fwrite (fcarr, sizeof (f4), nx * ny, f1);

  /* free temp arrays and close file */

  free (farr);
  free (fbarr);
  free (fcarr);
  fclose (f1);
  ier = 1;
  return ier;
}

i4 shift2d (double *arr, i4 nx, i4 ny, i4 ishift, i4 jshift)
{

/* Circular shift of the x,y indices of array *arr by ishift,jshift */
/* This function is similar to the shift function in IDL.  nx, ny
 * are the assumed dimensions of the array */

  double *temp;
  i4 i, j, ii, jj;
  temp = (double *) malloc (sizeof (double) * nx * ny);
  for (i = 0; i < nx; i++)
    {
      ii = (i + ishift) % nx;	/* ii = (i + ishift) modulo nx */

      for (j = 0; j < ny; j++)
	{
	  jj = (j + jshift) % ny;	/* jj = (j+jshift) modulo ny */

	  /* Now members of temp array get shifted: */

	  *(temp + ii * ny + jj) = *(arr + i * ny + j);
	}
    }

  /* Now copy temp array back into arr, then destroy temp and return */

  memcpy ((void *) arr, (void *) temp, nx * ny * sizeof (double));
  free (temp);
  return 0;
}

i4 maxloc (double *arr, i4 xsize)
{

/* finds the location of the maximum of the double array *arr and returns it. */

  i4 i, location;
  double amax;
  /* initialize amax and location to 0th element */
  amax = *(arr + 0);
  location = 0;
  for (i = 1; i < xsize; i++)
    {
      if (*(arr + i) > amax)
	{
	  amax = *(arr + i);
	  location = i;
	}
    }
  return location;
}

i4 imaxloc (i4 * arr, i4 xsize)
{

/* finds the location of the maximum of the i4 array *arr and returns it. */

  i4 i, location;
  i4 amax;
  /* initialize amax and location to 0th element */
  amax = *(arr + 0);
  location = 0;
  for (i = 1; i < xsize; i++)
    {
      if (*(arr + i) > amax)
	{
	  amax = *(arr + i);
	  location = i;
	}
    }
  return location;
}

i4 minloc (double *arr, i4 xsize)
{

/* finds the location of the minimum of the double array *arr and returns it. */

  i4 i, location;
  double amin;
  /* initialize amin and location to 0th element */
  amin = *(arr + 0);
  location = 0;
  for (i = 1; i < xsize; i++)
    {
      if (*(arr + i) < amin)
	{
	  amin = *(arr + i);
	  location = i;
	}
    }
  return location;
}

i4 iminloc (i4 * arr, i4 xsize)
{

/* finds the location of the minimum of the i4 array *arr and returns it. */

  i4 i, location;
  i4 amin;
  /* initialize amin and location to 0th element */
  amin = *(arr + 0);
  location = 0;
  for (i = 1; i < xsize; i++)
    {
      if (*(arr + i) < amin)
	{
	  amin = *(arr + i);
	  location = i;
	}
    }
  return location;
}

i4 interpcc2d (double *fdata, double xmiss, i4 nx, i4 ny, 
    double *xwant, i4 nxinterp, double *ywant, i4 nyinterp, double **finterp)
{
  /*
   * This function does cubic convolution interpolation onto an array 
   * finterp from data defined in array fdata.  nx, ny are the
   * assumed dimensions of fdata, and nxinterp, nyinterp are the
   * assumed dimensions of finterp.  The values of x,y at which
   * the interpolation is desired are passed in through the arrays
   * xwant and ywant, which are dimensioned nxinterp and nyinterp,
   * respectively.  It is assumed that xwant, ywant are in units of
   * the indices of the original data array (fdata), 
   * treated as floating point (double precision, actually) 
   * numbers. Arrays fdata, xwant, and ywant are passed in
   * as pointers; The array finterp is defined in this function
   * as a "double" pointer and the array is created and passed back to
   * the calling function.  In the calling function, finterp is declared
   * as a pointer, but when it is passed into this function as
   * an argument, the address of the pointer is used.
   * 
   * if any of the datapoints within a kernel weighting distance of
   * xwant and ywant are equal to xmiss,
   * the returned value of finterp is also set to xmiss.  xmiss is a user-
   * defineable calling argument.
   */

  double *cdata;
/*  double txt, tyt, xint, yint, ftmp, xmiss = 0.; */
  double txt, tyt, xint, yint, ftmp;

  /* Logic for a user-defined value of xmiss has been added.  Previously
   * was just set to 0 as a local variable */

  double tx, ty, rx, ry;
  i4 i, ii, j, jj, itemp, jtemp, izero, jzero, databad;
/*  i4 transp; */

  /* Now, create the cdata array, bigger by 1 gp than fdata
   * all around the borders: */

  cdata = (double *) malloc (sizeof (double) * (nx + 2) * (ny + 2));

  /* Now fill the interior of cdata with fdata */

  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
	{
	  *(cdata + (i + 1)*(ny + 2) + (j + 1)) = *(fdata + i*ny + j);
	}
    }

  /*
   * The basic concept for filling in edges and corners of cdata is this:
   * The edge point is equal to 3*(value of adjacent point)
   * -3*value(next to adjacent point) + 1*value(3rd point).  This
   * prescription yields an extrapolation which is consistent with
   * a 3rd (or is it 4th?) order Taylor expansion of the function
   * evaluated at the last real gridpoint, and extrapolated to the
   * edge point.  This procedure is followed
   * thoughout here, though I think it isn't really correct for the
   * corner points because there I think an expansion from both
   * both directions should be done.  But no harm seems to be done
   * to the results.
   */

  /* Fill in the edges of cdata: */

  for (j = 0; j < ny; j++)
    {

      /* left and right edges: */

      *(cdata + 0*(ny + 2) + (j+1)) = *(fdata + 2*ny + j)
	- 3. * (*(fdata + 1*ny + j)) + 3. * (*(fdata + 0*ny + j));

      *(cdata + (nx + 1)*(ny + 2) + (j + 1)) = *(fdata + (nx - 3)*ny + j)
	- 3. * (*(fdata + (nx - 2)*ny + j)) + 3. * (*(fdata + (nx - 1)*ny + j));
    }
  for (i = 0; i < nx; i++)
    {

      /* bottom and top edges: */

      *(cdata + (i + 1)*(ny + 2) + 0) = *(fdata + i*ny + 2)
	- 3. * (*(fdata + i*ny + 1)) + 3. * (*(fdata + i*ny + 0));

      *(cdata + (i + 1)*(ny + 2) + ny + 1) = *(fdata + i*ny + ny - 3)
	- 3. * (*(fdata + i*ny + ny - 2)) + 3. * (*(fdata + i*ny + ny - 1));
    }

  /* Now fill in the 4 corners: */

  *(cdata + 0*(nx + 2) + 0) = 
    3. * (*(cdata + 1*(ny + 2) + 0)) -
    3. * (*(cdata + 2*(ny + 2) + 0)) + *(cdata + 3*(ny + 2) + 0);

  *(cdata + (nx + 1)*(ny + 2) + 0) = 
    3. * (*(cdata + nx*(ny + 2) + 0)) -
    3. * (*(cdata + (nx - 1)*(ny + 2) + 0)) + *(cdata +
						    (nx - 2)*(ny + 2) + 0);

  *(cdata + 0*(ny + 2) + ny + 1) = 
    3. * (*(cdata + 0*(ny + 2) + ny)) -
    3. * (*(cdata + 0*(ny + 2) + ny - 1)) + *(cdata + 0*(ny + 2) + ny - 2);

  *(cdata + (nx + 1)*(ny + 2) + ny + 1) =
    3. * (*(cdata + nx*(ny + 2) + ny + 1)) -
    3. * (*(cdata + (nx - 1)*(ny + 2) + ny + 1)) + *(cdata +
						       (nx - 2)*(ny + 2) +
						       ny + 1);

  /* Now create the space for finterp */

  *finterp = (double *) malloc (sizeof (double) * nxinterp * nyinterp);

  /* Now interpolate onto the desired grid */

  for (i = 0; i < nxinterp; i++)
    {
      /* starting the outer loop over x */

      xint = *(xwant + i);

      /* make sure izero is in bounds */

      itemp = ((i4) xint > 0) ? (i4) xint : 0;
      izero = (itemp < (nx - 2)) ? itemp : nx - 2;
      for (j = 0; j < nyinterp; j++)
	{
	  /* starting the outer loop over y */

	  yint = *(ywant + j);
	  if ((yint < 0.) || (yint > (double) (ny - 1))
	      || ((xint < 0) || (xint > (double) (nx - 1))))
	    {
	      /* if data missing, set interp to xmiss */

/* Debug
              printf("interpccd2: i=%d,j=%d gets finterp[i,j] set to xmiss\n",
                    i,j);
*/
	      *(*finterp + i * nyinterp + j) = xmiss;
	    }
	  else
	    {
	      /* make sure jzero is in bounds */

	      jtemp = ((i4) yint > 0) ? (i4) yint : 0;
	      jzero = (jtemp < (ny - 2)) ? jtemp : ny - 2;

	      /* initialize the temporary finterp value */

	      ftmp = (double) 0.;

	      /* start the innermost loops over neighboring
	       * data points*/

              databad=0;
	      for (ii = -1; ii < 3; ii++)
		{
		  txt = xint - (double) (izero + ii);
		  tx = (double) fabs (txt);

		  /* evaluate kernel wt function r(tx): */

		  /* Note no testing for out of range 
		   * values of |tx| or |ty| > 2 --
		   * we assume the tx, ty are properly
		   * computed such that their absolute
		   * value never exceeds 2. */

		  rx = (tx >= (double) 1.0) ?
		    (((((double) (-0.5)) * tx +
		       ((double) 2.5)) * tx) -
		     (double) 4.) * tx + (double) 2. :
		    (((double) (1.5)) * tx -
		     ((double) (2.5))) * tx * tx + (double) 1.;

		  for (jj = -1; jj < 3; jj++)
		    {

		      tyt = yint - (double) (jzero + jj);
		      ty = (double) fabs (tyt);

		      /* evaluate kernel weighting
		       * function r(ty): */

		      ry = (ty >= (double) 1.0) ?
			(((((double) (-0.5)) * ty +
			   ((double) 2.5)) * ty) -
			 (double) 4.) * ty + (double) 2. :
			(((double) (1.5)) * ty -
			 ((double) (2.5))) * ty * ty + (double) 1.;

		      /* do the cubic convolution
		       * over the neighboring data
		       * points, using the x and
		       * y evaluated kernel weighting
		       * functions rx and ry: */

		      ftmp = ftmp +
			*(cdata + (izero + 1 + ii)*(ny + 2)
			  + jzero + 1 + jj) * rx*ry;
                      if( *(cdata+(izero+1+ii)*(ny+2)+jzero+1+jj) == xmiss)
                          databad=1;
		    }
		}
	      /* now assign this value to interpolated
	         array, unless one of the values was xmiss: */
              if(databad)
              {
/* Debug
                 printf("interpcc2d: i=%d,j=%d gives databad\n",i,j);
*/
                 *(*finterp + i*nyinterp + j) = xmiss;
              }
              else
              {
	         *(*finterp + i*nyinterp + j) = ftmp;
              }
	    }
	}
    }


/* DEBUG
  transp=1;
  writeimage("cdata.dat",cdata,nx+2,ny+2,transp);
*/

  /* free the cdata space */
  free (cdata);

  /* we're done */

  return 0;
}

i4 byteswap (unsigned char *arr, i4 arrsize, i4 nbpw)
/* Pretty simple:  arr is input array, which is byte-swapped in place,
   nbpw is the number of bytes per word, and arrsize is the size of the array
   (in units of nbpw bytes).  It is assumed that arr has
   already have been correctly defined and allocated in the calling program. */
{
  i4 i, j;
  unsigned char temp;
  for (i = 0; i < arrsize; i++)	/* the loop over the array elements */
    {
      for (j = 0; j < nbpw/2; j++)/* the loop over bytes in a single element */
	{
	  temp = *(arr + i*nbpw + (nbpw - j - 1));
	  *(arr + i*nbpw + (nbpw - j - 1)) = *(arr + i*nbpw + j);
	  *(arr + i*nbpw + j) = temp;
	}
    }
  return 0;
}

i4 is_big_endian ()
/* This function returns 1 if it is a big endian machine, 0 otherwise */
{
  const unsigned char fakeword[4] = { 0xFF, 0x00, 0xFF, 0x00 };
  i4 realword = *((i4 *) fakeword);
  if (realword == 0xFF00FF00)
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

