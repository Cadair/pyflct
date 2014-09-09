/*

 warp: http://solarmuri.ssl.berkeley.edu/overview/publicdownloads/software.html
 Copyright (C) 2009 Regents of the University of California

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
i4 readimage (char *fname, i4 *nx, i4 * ny, double **arr, i4 transp);
i4 read2images (char *fname, i4 * nx, i4 * ny, double **arr, double **barr,
		i4 transp);
i4 writeimage (char *fname, double *arr, i4 nx, i4 ny, i4 transp);
i4 write2images (char *fname, double *arr, double *barr, i4 nx, i4 ny,
		 i4 transp);
i4 write3images (char *fname, double *arr, double *barr, double *carr,
		 i4 nx, i4 ny, i4 transp);
i4 make_freq(double *k, i4 ndim);
i4 is_big_endian ();
i4 byteswap (unsigned char *arr, i4 arrsize, i4 nbpw);
i4 warp_frac2d(double *arr, double *delx, double *dely, double *outarr,
        i4 nx, i4 ny, i4 transp, i4 quiet);
i4 shift_frac2d(double *arr, double delx, double dely, double *outarr,
        i4 nx, i4 ny, i4 transp, i4 quiet);

/* end function prototypes */

int main (int argc, char *argv[])
{

/* BEGIN MAIN PROGRAM */

  char *version ="1.01   ";
  char ifile[100], sfile[100], ofile[100];
  i4 quiet, verbose;
  i4 nx, nxdel, ny, nydel;
  i4 ier, ibe;
  double *f1, *f2, *delx, *dely;
  i4 transp = 1; /* This flag is nonzero to transpose input/output arrays */

  /* CODE TO READ IN ARGUMENTS AND OPTIONS: */

  /* check to see if number args is in right range - 
   * if not, print syntax & quit */

  if ((argc < 4) || (argc > 6))
    {
      printf("warp: Version %s Copyright: 2009 University of California\n",
          version);
      printf("Authors: G.H. Fisher, B.T. Welsch, UCB Space Sciences Lab\n\n");
      printf
        ("Syntax: %s ifile sfile ofile [-q]\n\n",argv[0]);
      printf("ifile - contains 1 image for shifting\n");
      printf("sfile - contains 2 values (or arrays) of shifts (or warps)\n");
      printf("ofile - contains output image of shifted or warped array\n");
      printf("-q - optional flag to omit all non-error messages\n");
      exit (1);
    }

  /* GET THE 3 REQUIRED ARGUMENTS */

  /* get input image file name */

  strncpy (ifile, argv[1], 99);

  /* get input  file name */
  strncpy (sfile, argv[2], 99);

  strncpy (ofile, argv[3], 99);

  /*  GET OPTIONAL ARGUMENT quiet : */
  quiet=0;
  verbose=1;
  if(argc == 5)
    {
       if(!strncmp("-q",argv[argc-1],2))
       {
          quiet=1;
          verbose=0;
       }

    }

 /* DONE FINDING ARGUMENTS AND OPTIONS */

  /* determine if this is a large endian or small endian platform */
  ibe = is_big_endian ();
  if (verbose)
    {
      printf("warp: Version %s Copyright: 2009 University of California\n",
          version);
      printf("Authors: G.H. Fisher, B.T. Welsch, UCB Space Sciences Lab\n\n");
      if (ibe)
	{
	  printf ("warp: large endian machine; i/o not byteswapped\n");
	}
      else
	{
	  printf ("warp: small endian machine; i/o will be byteswapped\n");
	}
    }

  /* print out arguments and options */

  if (verbose) printf ("warp: ifile = %s\n", ifile);
  if (verbose) printf ("warp: sfile = %s\n", sfile);
  if (verbose) printf ("warp: ofile = %s\n", ofile);


  /*
   * read nx, ny, and return references to nx and ny to main prog. *
   * NOTE -- roles of nx, ny are reversed from IDL!!!! In the C version,
   * must work in transposed space.  Therefore transp is set to 1.
   */

  ier = readimage (ifile, &nx, &ny, &f1, transp);

  if (verbose)
    printf ("warp: from input file, nx = %d, ny = %d\n", nx, ny);


  ier = read2images(sfile, &nxdel, &nydel, &delx, &dely, transp);
  if ((nxdel == nx) && (nydel == ny)) 
  {
     /* in this case we do warping */
     f2=(void *)malloc(sizeof(double) *nx*ny);
     ier=warp_frac2d(f1,delx,dely,f2,nx,ny,transp,quiet);
     free(delx);
     free(dely);
     free(f1);
  }
  else if ((nxdel == 1) && (nydel == 1))
  {
     /* In this case just a uniform shift */
     f2=(void *)malloc(sizeof(double) *nx*ny);
     ier=shift_frac2d(f1,delx[0],dely[0],f2,nx,ny,transp,quiet);
     free(delx);
     free(dely);
     free(f1);
  }
  else
  {
     /* In this case we have an error */
     printf("warp: bad shift arrays: nx=%d, nxdel=%d, ny=%d, nydel=%d\n",
        nx,nxdel,ny,nydel);
     exit(1);
  }

  writeimage (ofile, f2, nx, ny, transp);

  /* free the gaussian mask array, the original images, and the
   * velocity arrays */

  if (verbose)
    printf ("\nwarp: finished\n");

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


i4 make_freq(double *k, i4 ndim)
{
/* k is assumed already allocated in main program, with dimension ndim */
i4 n21,i,inext;
n21=(ndim/2)-1;
/* increasing k for the 1st half */
for (i=0;i<n21+1;i++)
{
	k[i]=(double)i;
}

/* now compute k for the 'middle' points */
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

/* now do negative k's for the upper index range, decreasing in amplitude */
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

i4 readimage (char *fname, i4 * nx, i4 * ny, double **arr,
	     i4 transp)
/* Function to read array dims, create space for 1 double array, read it in.
 * Note the double pointer to the double precision * arrays in the calling 
 * argument.  Note also these are referenced as pointers and returned to the
 * calling function. */
{
  FILE *f1;			/* pointer to input file */
  i4 newsize;			/* size of the new double prec. array to be 
				   read in =nx*ny */
  i4 i, ier, ibe, ise, vcid;
  f4 *farr;
  ibe = is_big_endian ();
  ise = 0;
  if (ibe == 0) ise = 1;	/* set small endian flag if not big  endian */
  f1 = fopen (fname, "rb");	/* open file for binary unformatted read */
  ier = 0;
  if (f1 == NULL)		/* error exit if file can't be opened */
    {
      printf ("readimage: cannot open file %s\n", fname);
      exit (1);
    }

  /* check that files begins with special vel_ccor flag: */
  fread (&vcid, sizeof (i4), 1, f1);
  if (ise) byteswap ((void *) &vcid, 1, sizeof (i4));
  if (vcid != 2136967593)
    {
      printf ("readimage: input file is not a vel_ccor i/o file\n");
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
      printf ("readimage: memory request for arr failed\n");
      exit (1);
    }
/*
	printf("%d bytes of memory now allocated for arr \n",newsize);
*/

  /* now read in the arr array */

  fread (farr, sizeof (f4), (*nx) * (*ny), f1);
  if (ise) byteswap ((void *) farr, (*nx) * (*ny), sizeof (f4));
  /*byteswap if needed */


  /* now transfer data from temp. arrays to arr: */

  for (i = 0; i < (*nx) * (*ny); i++)
    {
      *(*arr + i) = (double) *(farr + i);
    }

  /* now free the temp. arrays and close the files */

  free (farr);
  fclose (f1);
  ier = 1;
  return ier;
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

i4 shift_frac2d(double *arr, double delx, double dely, double *outarr,
        i4 nx, i4 ny, i4 transp, i4 quiet)

{

/* uncomment next line if to be compiled by itself */
/* #include <fftw3.h> */

  i4 i,j, verbose;
  double *kx, *ky, *fftdeltrx, *fftdeltry,*fftdeltix,*fftdeltiy;
  double normfac, dxarg, dyarg,fftdeltr,fftdelti,outar,outai;
  double pi=3.1415926535898;
  fftw_complex *ina, *outa;
  fftw_plan pa, pbacka;
  /* allocate memory for FFT arrays and intermediate arrays */
  outa = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
     nx * ny );
  ina = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * nx * ny);
  kx=(double *) fftw_malloc (sizeof(double)*nx);
  fftdeltrx=(double *) fftw_malloc (sizeof(double)*nx);
  fftdeltix=(double *) fftw_malloc (sizeof(double)*nx);
  ky=(double *) fftw_malloc (sizeof(double)*ny);
  fftdeltry=(double *) fftw_malloc (sizeof(double)*ny);
  fftdeltiy=(double *) fftw_malloc (sizeof(double)*ny);
  /* set up plans for FFTs: */
  pa = fftw_plan_dft_2d (nx, ny, ina, outa, FFTW_FORWARD, FFTW_MEASURE);
  pbacka = fftw_plan_dft_2d (nx, ny,outa, ina, FFTW_BACKWARD, FFTW_MEASURE);
  /* compute spatial frequencies kx and ky */
  make_freq(kx,nx);
  make_freq(ky,ny);
  if(quiet)
  {
    verbose=0; /* no io to stdout */
  }
  else
  {
    verbose=1; /* flag for sending any io to stdout */
  }
 /* calculate normalization factor */
  normfac = (1. / ((double) nx * ny));
 /* assign image fft array (complex) to input data array (double) */
      for (i = 0; i < nx * ny; i++)
	{
	  ina[i][0] = arr[i]; /* real part */
	  ina[i][1] = (double) 0.; /* imag part */
	}

  /* do the forward FFT: */
  
  fftw_execute (pa);
  
  /* get shifts normalized to nx, ny.  Note if transp==1 then roles of delx,
  dely are switched.  This occurs if array is column-major.  */

  if (transp)
  {
     dxarg= -dely/((double) nx);
     dyarg= -delx/((double) ny);
  }
  else
  {
     dxarg= -delx/((double) nx);
     dyarg= -dely/((double) ny);
  }

  /* calculate FFT of delta function at delx, dely, then mpy by FFT of image */

    /* calculate sin, cosine factors first */
    for (i=0;i<nx;i++)
    {
        fftdeltrx[i]=cos(2.*pi*kx[i]*dxarg);
        fftdeltix[i]=sin(2.*pi*kx[i]*dxarg);
    }
    for (j=0;j<ny;j++)
    {
        fftdeltry[j]=cos(2.*pi*ky[j]*dyarg);
        fftdeltiy[j]=sin(2.*pi*ky[j]*dyarg);
    }
    /* now compute fft of shifted image */
    for (i=0;i<nx;i++)
    {
      for (j=0;j<ny;j++)
      {
         /* real part of exp(i kx dxarg + i ky dyarg) */
         fftdeltr=fftdeltrx[i]*fftdeltry[j]-fftdeltix[i]*fftdeltiy[j];
         /* imag part of exp(i kx dxarg + i ky dyarg) */
         fftdelti=fftdeltrx[i]*fftdeltiy[j]+fftdeltix[i]*fftdeltry[j];
         outar=outa[i*ny+j][0];
         outai=outa[i*ny+j][1];
         /* real part of fft of shifted fn */
         outa[i*ny+j][0]=(outar*fftdeltr-outai*fftdelti);
         /* imag part of fft of shifted fn */
         outa[i*ny+j][1]=(outar*fftdelti+outai*fftdeltr);
      }
    }

   /* Do inverse FFT to get shifted image fn back into the ina complex array */

    fftw_execute (pbacka);

   /* Output array is the real part of ina (imag. part should be zero) */

  for (i = 0; i < nx * ny; i++)
    {
      outarr[i]=(double) (ina[i][0])*normfac;
    }

/* Free the plans and locally created arrays */

      fftw_free (outa);
      fftw_free (kx);
      fftw_free (ky);
      fftw_free (fftdeltrx);
      fftw_free (fftdeltix);
      fftw_free (fftdeltry);
      fftw_free (fftdeltiy);
      fftw_free (ina);
      fftw_destroy_plan (pa);
      fftw_destroy_plan (pbacka);

/* done */

return 0;

}

i4 warp_frac2d(double *arr, double *delx, double *dely, double *outarr,
        i4 nx, i4 ny, i4 transp, i4 quiet)

{

/* uncomment if this function will be compiled on its own */
/* #include <fftw3.h> */

  i4 i,j,ii,jj,verbose;
  double *kx, *ky, *snkx, *cskx, *snky, *csky;
  double pi=3.1415926535898;
  double normfac, fftdeltr, fftdelti,outar,outai,totalimij,nxinv,nyinv,
         xarg,yarg,shiftx,shifty;
  fftw_complex *ina, *outa; 
  fftw_plan pa;
  /* allocate FFT input, output arrays, plus temporary arrays */
  outa = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) *
     nx * ny);
  ina = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * nx * ny);
  kx=(double *) fftw_malloc (sizeof(double)*nx);
  ky=(double *) fftw_malloc (sizeof(double)*ny);
  cskx=(double *) fftw_malloc (sizeof(double)*nx);
  snkx=(double *) fftw_malloc (sizeof(double)*nx);
  csky=(double *) fftw_malloc (sizeof(double)*ny);
  snky=(double *) fftw_malloc (sizeof(double)*ny);
  /* set up plan for forward FFT: */
  pa = fftw_plan_dft_2d (nx, ny, ina, outa, FFTW_FORWARD, FFTW_MEASURE);
  /* calculate spatial frequencies */
  make_freq(kx,nx);
  make_freq(ky,ny);
  if(quiet)
  {
    verbose=0; /* no io to stdout */
  }
  else
  {
    verbose=1; /* send io to stdout */
  }
  /* calculate normalization factor */
  normfac = (1. / ((double) nx * ny));
  nxinv = 1./((double) nx);
  nyinv = 1./((double) ny);
  /* copy input array to real part of input fft array */
  for (i = 0; i < nx * ny; i++)
  {
    ina[i][0] = (double) arr[i]; /* real part */
    ina[i][1] = (double) 0.; /* imag part */
  }

  /* actually do the forward FFT: */

  fftw_execute (pa);

  /* outer loop over spatial locations */
    for (i=0;i<nx;i++)
    {
      if(verbose)
      {
         printf ("warp: progress  i = %d out of %d\r", i, nx - 1);
         fflush(stdout);
      }
      for (j=0;j<ny;j++)
      {
         /* Note that if (transp) then array is assumed column major, switch
            roles of delx, dely */
         if(transp)
         {
            shiftx=dely[i*ny+j];
            shifty=delx[i*ny+j];
         }
         else
         {
            shiftx=delx[i*ny+j];
            shifty=dely[i*ny+j];
         }
         xarg=2.*pi*(((double) i) - shiftx)*nxinv;
         yarg=2.*pi*(((double) j) - shifty)*nyinv;
         for (ii=0;ii<nx;ii++)
         {
             /* compute kx-dependent sin, cos terms */
             cskx[ii]=cos(kx[ii]*xarg);
             snkx[ii]=sin(kx[ii]*xarg);
         }
         for (jj=0;jj<ny;jj++)
         {
             /* compute ky-dependent sin, cos terms */
             csky[jj]=cos(ky[jj]*yarg);
             snky[jj]=sin(ky[jj]*yarg);
         }

         /* initialize the integral over wavenumber */
         totalimij=(double) 0.;

         /* inner loop to integrate over kx, ky wavenumbers */
         for (ii=0;ii<nx;ii++)
             {
               for (jj=0;jj<ny;jj++)
                 {
                    /* compute real, im parts of exp (i [kx xarg + ky yarg]) */
                    fftdeltr=cskx[ii]*csky[jj]-snkx[ii]*snky[jj];
                    fftdelti=cskx[ii]*snky[jj]+snkx[ii]*csky[jj];
                    /* extract real, imag parts of fft of original image */
                    outar=outa[ii*ny+jj][0];
                    outai=outa[ii*ny+jj][1];
                    /* add contributions to the total shifted image at i,j 
                       noting that only the real part matters */
                    totalimij+= (outar*fftdeltr-outai*fftdelti);
                 }
             }
         
         /* record shifted image at i,j, after normalizing */
         outarr[i*ny+j] = (double) totalimij*normfac;
      }
    }

/* Free the plans and locally created arrays */

      fftw_free (outa);
      fftw_free (kx);
      fftw_free (ky);
      fftw_free (cskx);
      fftw_free (snkx);
      fftw_free (csky);
      fftw_free (snky);
      fftw_free (ina);
      fftw_destroy_plan (pa);

return 0;
}
