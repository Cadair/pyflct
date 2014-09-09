# Fourier Local Corellation Tracking (FLCT)

April 10, 2009 - George H. Fisher and Brian T. Welsch, SSL, UC Berkeley

## DESCRIPTION OF THE FLCT CODE FUNCTION

The flct code is designed to estimate a 2-d velocity field from two
images, under the assumption that the 2nd image differs from the first one
from changes imposed by the velocity field.  The /doc folder in the distribution
contains files that describe the details of the methods underlying the 
code and how the code is run.  Specifically, the /doc folder includes the 
unix man page flct.1, and postscript and .pdf versions of the man page 
viewable from non-unix (eg MS Windows) systems.  The /doc folder also 
includes a copy of the Fisher and Welsch 2008 PASP paper describing the flct 
code as it existed in the Spring of 2007 (version test_13).

## LICENSE

The flct source code is open-source, under the GPL version 2 license.  You
are free to copy the code and use it as you like (in accordance with GPL v.2), 
but if you find it useful in any work you publish, we would greatly appreciate 
citations to the papers Fisher & Welsch PASP 383, 373, (2008), and 
Welsch et al (ApJ 610, 1148) (2004), and any future papers which describe 
updated versions of flct.  If you incorporate flct source code into any 
other programs, the GPL v. 2 licensing agreement provisions will apply.

## VERSION

This current version of flct is 1.01-1.  The C-source code for flct is
unchanged from version 1.01.  The only significant change from 1.01 is
(1) the addition of a pre-compiled binary for solaris running on x86_64 (i86pc)
hardware; (2) the modification of the vcimage1out.pro, vcimage2out.pro, and 
vcimage3out.pro IDL procedures to allow for the output of single scalar 
values, as well as the output of 2d image arrays; and (3) the addition of the 
source code for the "warp" executable, which is a C version of the IDL 
procedure shift_frac2d.pro.  This program can perform image 
warping significantly faster than shift_frac2d.pro.  Instructions for 
compiling warp are included below in the instructions for compiling flct.  
The documentation for running warp is included in the /docs folder as a man
page, plus .ps and .pdf copies of the man page.  No pre-compiled executables 
of warp are included.

This version (1.01 and 1.01-1) differ from 1.0 (final) in that
the temporary arrays g1tmp and g2tmp were eliminated from the calculation of
the "sliding box" arrays g1 and g1.  This resulted in a significant speedup over
version 1.0.  This version (1.01) supercedes all of the earlier 
versions, including all the test_* versions, and specifically including
test_13 that is described in the Fisher and Welsch paper flct_technique.pdf
in the /doc folder.

The shift_frac2d.pro function included in the IDL-io-procedures folder
uses real_pt instead of real_part so that it can be used without 
licensing issues in GDL.  

## MAIN CHANGES FROM VERSION test_13 (version described in Fisher & Welsch 2008)

When sub-images are extracted from the two images used for correlation tracking,
a local mean subtraction is now done before the sub-images are muliplied by the
gaussian windowing function.  This change greatly improves the behavior with
white-light images and other images with significant non-zero biases, and
in conjunction with the low-pass filtering option, greatly improves the
behavior for very small applied shifts (shifts << 1 pixel).

For example, the Fisher and Welsch paper showed the breakdown of version 
test_13 of flct when a .01 degree rotation was applied to a magnetogram image.  
With version 1.0, an accurate reconstruction of the velocities was found for
a .0001 degree rotation, a factor of 100 smaller than the smallest (.01 degree) 
rotation considered in Fisher and Welsch.

A new option to the code was added to allow for "skipping" points in the
calculation of LCT derived velocities.  This is useful for large images when
it may not be necessary to compute an LCT velocity at every pixel.  This option
also allows for interpolation of the skipped points from the computed points
using cubic convolution interpolation.

We are grateful to Karin Muglach of NRL for pointing out the usefulness of
having a "skipping" option, and for pointing out the desireability of local mean
subtraction.

The IDL i/o procedures and functions designed to be used in conjunction with
flct have been updated to include a few error checks.  In addition, a new
IDL function, shift_frac2d.pro, has been added, to facilitate the accurate
shifting of images by non-integer numbers of pixels.  The syntax of
shift_frac2d is essentially the same as the IDL shift function.
Shift_frac2d also includes the ability to perform non-uniform shifts (warping) 
of images.  While this functionality is handled very accurately, 
it is also slow, and may be impractical for large images.  
On the other hand, shift_frac2d is very fast for applying simple, uniform 
non-integer shifts.  An obvious application of shift_frac2d is the removal 
of secular shifts, such as those due to solar rotation between two images 
on the Sun taken at different times.

## INSTALLATION

The flct code is written in C, although it is designed to be used in
conjunction with IDL for performing the binary I/O that flct uses to read
in the images and write out the velocity fields.

For the convenience of users, we have attempted to pre-compile executable,
stand-alone binaries of flct for a variety of commonly used computing
platforms, which can be found in the /bin directory of the distribution.

However, the most reliable way to install flct is to compile the code 
yourself.  To do that, you first need to make sure that version 3 of fftw is 
installed.

If fftw is not installed, go to http://www.fftw.org and download the 
latest stable version (currently 3.1.2) as a tarball.  Unpack the tarball,
get into the distribution directory, and then type the command ./configure 
followed by make .  Then, as root, type make install.  If you do not have 
root priviledge, you can skip the last step, but will then have to edit 
the flct Makefile as described below.  flct currently assumes that version
3 of fftw will be used, as opposed to the much older version 2.

Once fftw3 has been installed, unpack the flct tarball, and go into the
source directory.  If the fftw3 library is in /usr/local/lib , and the
fftw3 include files are in /usr/local/include , then you should just be able
to type "make" and flct should compile.  If you have root priviledge, and
are happy to install flct into /usr/local/bin and the man page into
/usr/local/man/man1, then become root and type "make install".  Becoming
root and typing "make uninstall" will uninstall flct and the man page.  Note
that the default man page locations are different in Mac OSX (see comments in
the Makefile).  

If the fftw3 library and include files are not located in these places,
you will have to edit Makefile.  If you were unable to install fftw3 as root,
but were able to compile the fftw source code, you can set the LIBFFTW3
variable to point to the .libs subdirectory of the main fftw3 build
directory, and set the INCLUDEFFTW3 variable to point to the api subdirectory
of the build directory.  After this flct should compile.  If you do not
want to put flct in /usr/local/bin and the man page in /usr/local/man/man1,
edit Makefile and change the FLCT_MANDIR and FLCT_BINDIR variables to
point to where you do want them to go.  Then typing make install should
install flct and its man page where you like.  Typing man flct should show 
the flct man page, describing how to run the code.

If you would like to compile warp (not required) after compiling flct, 
type "make warp" after typing make.  To install the warp executable and 
man page, type "make warpinstall".  To uninstall warp and its manpage, type
"make warpuninstall".  The warp executable performs the same function as
the IDL procedure shift_frac2d.pro, but for image warping, is considerably
faster.

If you do not wish to compile flct yourself, we have included pre-compiled
binaries of flct for several different computing platforms, including linux
on x86 and x86_64 processors, MS Windows XP (win32), Mac OSX (x86), 
Mac OSX (G4, G5), and Sun solaris-sparc and Sun solaris-i386pc 
platforms.  Not all of these have been tested, so we cannot vouch for their 
functionality.  To install a pre-compiled binary, go into the bin folder, 
find your platform, and copy the corresponding executable to somewhere in your 
path.  If you then type flct<enter> you should get a page of output 
consisting of a short summary of flct documentation.  You can manually 
copy the flct.1 file in /docs to your default man page location.
No pre-compiled versions of warp are included.

For MS windows users, the flct binary was compiled with a MINGW gcc compiler
on a 32-bit windows xp platform.  We know that this does work with at least
some versions of Vista, but are unsure about older versions of MS Windows at
this time.  It is important to note that flct.exe
in Windows is a console application, not a standard windows application.  It
must be run from the windows command-prompt tool, ie the black DOS-like
window.  The syntax for its use is identical to the unix/linux version.
It will also work if spawned from IDL, but only if flct.exe is installed
somewhere in your default path.  For lazy MS Windows users, this can
be accomplished by copying flct.exe to C:\WINDOWS\SYSTEM32 .  A better
solution is to modify your PATH variable in windows to include the location
where you have copied the file.

## TESTING

The tests/ subdirectory of the flct distribution contains 3 files, 
hashgauss.dat, deltaxygauss.dat, and testgaussvel.dat .  hashgauss.dat
contains two input images, which we'll call f1 and f2, upon which flct 
can be run and tested.  The file deltaxygauss.dat contains two 2-D arrays 
of x and y shifts that were applied to f1 to generate f2.  The file 
testgaussvel.dat contains the output velocity field that flct generated for 
a specific set of parameters, derived from the f1 and f2 images.  
Ideally, the derived velocity field should correspond closely to the 
arrays of x and y shifts.

To ensure that your version of flct is working correctly, you can compare
its performance to the output from the files.

Make sure that all of the IDL functions and procedures in the 
IDL-io-procedures/ folder are placed in your IDL path.

Open an IDL window in the tests folder, and type
vcimage2in,f1,f2,'hashgauss.dat'

If you type help, you should see two floating point arrays, f1 and f2,
dimensioned 201,101.

Now, type 
vcimage2in,deltax,deltay,'deltaxygauss.dat'

If you type help, you should also see floating point arrays deltax, and deltay,
also dimensioned 201,101.  If you type

shade_surf,deltax,chars=2

you should see a gaussian in the center of the image.  If you type
shade_surf,deltay,chars=2

you should see a negative gaussian in the center, in the same location as
the other gaussian.

To perform a test run of flct, enter this command in your IDL window:
$flct hashgauss.dat testvels.dat 1 1 5 -k 0.5

You should see some output printed to the screen, and a dynamically
changing counter while the code is running.  It will tell you when it
is finished.

To read in the output velocity field computed by flct, type this command in
your IDL session:

vcimage3in,vxtest,vytest,vmtest,'testvels.dat'

After typing help, you should now see 3 additional arrays, vxtest,vytest,
vmtest, all dimensioned 201,101.

To compare this with the applied shifts, you can compare with
these commands:

shade_surf,vxtest,chars=2
shade_surf,deltax,chars=2

For a scatterplot of the derived velocities with the applied shifts, you
can type

plot, deltax,vxtest,ps=3
plot, deltay,vytest,ps=3

To better examine the behavior for small shifts, you can look at a log-log
plot:

plot_oo, deltax,vxtest,ps=3

If you want to compare your results with the results we achieved to make
sure that flct is functioning consistently, you can read in the precomputed flct
results as follows:

vcimage3in,vx,vy,vm,'testgaussvel.dat'

and then compare your derived velocities, vxtest and vytest, with vx and vy
as read in from testgaussvel.dat.

How were the test images f1 and f2 generated?  f1 is a field of random numbers
uniformly distributed between 0 and 1.  It was generated in IDL with the
command

f1=randomu(seed,201,101)

We have found this to be an extremely challenging test image, 
and it has posed great difficulties for earlier versions of flct, 
especially for shifts of a fraction of a pixel.  f2 is simply a non-uniformly
shifted version of f1, where the shifts are given by the deltax and deltay
arrays.  The IDL procedure shift_frac2d, located in the IDL-io-procedures 
folder, was used to generate f2 from f1 and deltax, deltay via this command:

f2=shift_frac2d(f1,deltax,deltay,/progress)

