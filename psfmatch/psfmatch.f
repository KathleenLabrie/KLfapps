c PSF-match two IRAF images (real or simulated). Compensates for small 
c shifts, photometric zeropoints, offsets, as well as PSF differences!
c Program makes the first-mentioned image look like the second image,
c or just prints out info about the images if no output image specified.

c Usage:  psfmatch fnm0 fnm1 zmin zmax [fnmout]
c
c where fnm0 = name of first input image (the one to be convolved, normally
c	       the good seeing image)
c 	fnm1 = name of second input image (the bad seeing one, to be matched.
c              this is the reference image that is unchanged)
c	zmin, zmax = peak intensity range for PSF stars, after sky subtraction
c       fnmout = output blurred image name - if not given then the program 
c              computes/prints out information only, without producing
c	       output images. (Convolving the image takes some time.)
c
c Difference image goes to "diff_psfm". (All output images are overwritten if
c they already exist.) Summary output goes to 'psfmatch.log'.

c Hardwired parameters: 
c   nmax = max image size
c   nkxy= kernel size - kernel is nkxy x nkxy
c   nxystar= postage stamps around each star are nxystar x nxystar
c	(has to be bigger than about 2 x nkxy, I think)
c   nstarmax= max # of stars.
c
c	parameter (nmax=4600000, nkxy=11, nxystar=25, nstarmax=1000)
	parameter (nmax=9200000, nkxydef=55, nxystardef=101, 
     &		nstarmax=1000)
c
	character fnm0*80, fnm1*80, fnmout*80, fnmdiff*80, fcoo*80,
     &	   fnmdiffdef*80, str*80, ans*1, imtype0*5, imtype1*5, fkern*80
       character*80 otherarg(8)
	real z0(nmax), z1(nmax), z0s((nxystardef+nkxydef+5)**2),
     &    z1s((nxystardef+nkxydef+5)**2),
     &	   diffs((nxystardef+nkxydef+5)**2), z2(nmax),
     &	   z2s(nxystardef**2), diff(nmax)
	real xstar(nstarmax), ystar(nstarmax), kernel(nkxydef,nkxydef),
     &    kernel1(nkxydef,nkxydef), chisq(nstarmax), sigma(nstarmax),
     &    z0smax(nstarmax), temp(3*nstarmax), fwhm0(nstarmax),
     &    fwhm1(nstarmax), dx(nstarmax), dy(nstarmax)
	real z0s2(nxystardef,nxystardef), z1s2(nxystardef,nxystardef)
	integer iuse(nstarmax)
	integer i,j
	integer nkxy, nxystar
	real foo
	equivalence (z0s2,z0s), (z1s2,z1s)
	logical debug, doit, sym, findstars
CCC KL
	logical inverse
CCC endKL
	common /fit/ debug, nxfit, nyfit, sym
	data fnmdiffdef/'diff_psfm.imh'/, nitermax/2/, debug/.false./
	data fkern/'kernel.dat'/
	data sym/.true./, findstars/.true./
CCC KL
	data inverse/.false./
CCC endKL

	if (debug) call pgbegin (1, '/xt', 2, 2)
cccc	if (debug) call pgask (.false.)

c OK, here is where program data is extracted from command line. Note 
c that the fifth argument is optional - if not given then the program
c computes info about the images, prints out stuff (including summary 
c output to file 'psfmatch.log'), but doesn't output a convolved or 
c difference image.
c ----------------------------------------------------------------------
	narg=iargc()
	if (narg .ge. 6) then
	  i=1
	  j=0
	  do while (i .le. narg)
	     call getarg (i, str)
	     if (str(1:index(str,' ')) .eq. '-f') then
	     	i=i+1
		call getarg (i,fcoo)
		findstars = .false.
	     else if (str(1:index(str,' ')) .eq. '--inverse') then
	        inverse = .true.
	     else
	     	j=j+1
	     	otherarg(j) = str
	     endif
	     i=i+1
	  enddo
	else
	  write(6,'(a)')'Usage:   psfmatch fnm0 fnm1 pkmin pkmax nkxy '//
     &		'nxystar [fnmout] [fnmdiff] [--inverse]'
     	  write(6,'(a)')'	or'
	  write(6,'(a)')'Usage:   psfmatch fnm0 fnm1 -f starcoo nkxy '//
     &		'nxystar [fnmout] [fnmdiff] [--inverse]'
	  stop
	end if

	fnm0 = otherarg(1)
	fnm1 = otherarg(2)
	if (findstars) then
	  zmin0 = str2real (otherarg(3))
	  zmax0 = str2real (otherarg(4))
	  i=5
	else
	  i=3
	end if
	nkxy = nint(str2real (otherarg(i)))
	nxystar = nint(str2real (otherarg(i+1)))
	if (narg .ge. 7) then
	  doit=.true.
	  fnmout = otherarg(i+2)
	  if (narg .ge. 8) then
	    fnmdiff = otherarg(i+3)
	  else
	    fnmdiff=fnmdiffdef
	  endif
	else
	  doit=.false.
	  fnmout=' '
	endif

c -------------------------------------------------------------------

c Initialization ...

	do 71 is=1,nstarmax
	iuse(is)=1
71	chisq(is)=0.
	chisqrejmax=1.e30

c Load initial BIG arrays ...

	write(6,'(a)')'Reading images ...'
c   First image
	ind = rsubstr(fnm0,'.')
	imtype0 = fnm0(ind:)
	if ( imtype0(1:length(imtype0)) .eq. '.imh') then
         call iraf_read_r4 (fnm0, z0, nmax, nx, ny, ixoff, iyoff, ier)
	else if ( imtype0(1:length(imtype0)) .eq. '.fits') then
	  call fits_read_r4 (fnm0, z0, nmax, nx, ny, ixoff, iyoff, ier)
	else
	  write(*,*) 'ERROR: Unknown image format.'
	  goto 99
	endif
c   Second image
	ind = rsubstr(fnm1,'.')
	imtype1 = fnm1(ind:)
	if ( imtype1(1:length(imtype1)) .eq. '.imh') then
         call iraf_read_r4 (fnm1, z1, nmax, nx, ny, ixoff, iyoff, ier)
	else if ( imtype1(1:length(imtype1)) .eq. '.fits') then
	  call fits_read_r4 (fnm1, z1, nmax, nx, ny, ixoff, iyoff, ier)
	else
	  write(*,*) 'ERROR: Unknown image format.'
	  goto 99
	endif

c Get stars for use, once and for all ... 

	call rmnsd_sampled (z0, nx*ny, sky0, skysig0)
	call rmnsd_sampled (z1, nx*ny, sky1, skysig1)
	if (findstars) then
	  write(*,*) 'Searching for isolated stars...'
	  call find_isolated_new (z0, nx, ny, zmin0, zmax0, xstar, ystar,
     &      nstarmax, nstar)
       else
	  write(*,*) 'Reading star coordinates out of ',
     &		fcoo(1:length(fcoo))
	  call read_coo (fcoo, xstar, ystar, nstarmax, nstar)
	endif
	write(6,'(a,f10.2,a,f7.2)')'1st image: sky level =',sky0,
     &    '  sky noise =',skysig0
	write(6,'(a,f10.2,a,f7.2)')'2nd image: sky level =',sky1,
     &    '  sky noise =',skysig1
	write(6,'(a,i5)')'Number of isolated stars =',nstar
	if (nstar .le. 0) stop 'Error: no stars found.'

c Loop over all stars, compute mean kernel ... 

	do 51 iter=1,nitermax
	write(6,'(/a,i2,a)')'Computing kernel - pass',iter,' ...'

c Initialize only ... (only important argument is the last one!)
	write(*,*) 'here 1'

	call getkern (z0s, z1s, nxystar, nxystar, kernel1, nkxy, nkxy, 
     &    dsky, 1)
	nuse=0
	write(*,*) 'here 2'
	
c Add stars to arrays only, and do fit to profile if doit=.false..
c (Used to do this only for doit=.false. to speed things up, but
c this is actually a valuable indicator of whether a star is OK
c or not.) 

	do 11 is=1,nstar
	write(*,*) 'here 3'
	write(str,'(a,i3,a)')'Star',is,' - '
CCC	if ((iter .eq. 1) .and.
	if ((iter .eq. 1 .or. chisq(is) .le. chisqrejmax) .and.
     &    iuse(is) .ge. 1) then
	write(*,*) 'here 4'
     	  debug=.true.
	  if (debug) write(6,'(a,i3,a,i2,a,2i5)') 'Star',is,' pass',
     &      iter,'  (xc,yc)=',nint(xstar(is)), nint(ystar(is))
   	  debug=.false.
	  call subarray (z0, nx, ny, nint(xstar(is)), nint(ystar(is)),
     &      z0s, nxystar, nxystar)
	  call subarray (z1, nx, ny, nint(xstar(is)), nint(ystar(is)),
     &      z1s, nxystar, nxystar)
c Changed - always do fit to profile as it is a useful indicator as to
c whether a star can be used for the fit.
	  call fitmoffat2d (z0s, nxystar, nxystar, sky0, fwhm0(is), 
     &      xc0, yc0, str(1:length(str))//fnm0(1:length(fnm0)), ierr0)
	  call fitmoffat2d (z1s, nxystar, nxystar, sky1, fwhm1(is), 
     &      xc1, yc1, str(1:length(str))//fnm1(1:length(fnm1)), ierr1)
	  dx(is)=xc1-xc0
	  dy(is)=yc1-yc0
CCC	  if (ierr0 .lt. 0 .or. ierr1 .lt. 0) iuse(is)=0
c Only use star if the fit to its profile worked ok ...
CCC	  if (iuse(is) .gt. 0) then
	    call getkern (z0s, z1s, nxystar, nxystar, kernel1, nkxy,  
     &        nkxy, dsky, 2)
	    nuse=nuse+1
CCC	  end if
	end if
11	continue

c Start debug
	write(*,*) 'Number of stars used=',nuse
c End debug

c Compute kernel ... (only important argument is the last one!)
	if (nuse .lt. 1) stop 'Error: no stars to compute kernel!'
	call getkern (z0s, z1s, nxystar, nxystar, kernel1, nkxy, nkxy, 
     &    dsky, 4)

c Stats on "fit" for each star ... Don't need full blown convolution if 
c you're just doing this on a few stars to check quality of fit (but of
c course you'll have to do the big convolution later!). Also note that 
c this is photon-noise limited: noise for a given pixel of a star is
c scaled from sky noise using sqrt(N) stats.

	do 21 is=1,nstar
	if (iuse(is) .le. 0) then
	  chisq(is)=9999.99
	  sigma(is)=9999.99
	  ans='X'
	else
	  nxystar1=nxystar+nkxy
	  call subarray (z0, nx, ny, nint(xstar(is)), nint(ystar(is)), 
     &      z0s, nxystar1, nxystar1)
	  call subarray (z1, nx, ny, nint(xstar(is)), nint(ystar(is)),  
     &      z1s, nxystar1, nxystar1)
	  call conv2d (z0s, nxystar1, nxystar1, kernel1, nkxy, nkxy, z2)
	  do 215 i=1,nxystar1**2
215	  diffs(i)=z1s(i)-z2(i)-dsky
	  call subarray (diffs, nxystar1, nxystar1, nxystar1/2+1, 
     &      nxystar1/2+1, temp, nxystar, nxystar)
	  call mv4 (temp, diffs, nxystar**2)
	  call minmax (z0s, nxystar1**2, z0smin, z0smax(is)) 
	  chisqsum=0.
	  sumsq=0.
	  do 23 i=1,nxystar**2
	  sumsq=sumsq+diffs(i)**2
c Following little kludge is needed to avoid problems with bogus pixels.
c (Otherwise a simple 1 line statement would be OK). This is quite a 
c serious problem otherwise!
	  sig=skysig0
	  if (z0s(i) .gt. sky0 .and. sky0 .gt. 0.) 
     &      sig=sqrt(z0s(i)/sky0)*skysig0
23	  chisqsum=chisqsum+(diffs(i)/sig)**2
	  chisq(is)=chisqsum/nxystar**2
	  sigma(is)=sqrt(sumsq/nxystar**2)
	  ans=' '
	end if
	if (debug) write(6,'(3x,i3,2i6,f10.1,6f9.2,4x,a1)')
     &    is, nint(xstar(is)), nint(ystar(is)), z0smax(is)-sky0,
     &    chisq(is), sigma(is), fwhm0(is), fwhm1(is), dx(is), dy(is), 
     &    ans
21	continue

c Do rejection ... Find median, then reject all stars with chisq > twice
c median chisq, before looping back and redoing kernel determination. 

	if (iter .ne. nitermax) then
	  rmdnchisq=rmedian1(chisq,iuse,nstar,temp)
c	  chisqrejmax=3.*rmdnchisq
	  chisqrejmax=2.*rmdnchisq
	  write(6,'(2(a,f10.3))')'Median chisq =',rmdnchisq,
     &      ',  rejection chisq =',chisqrejmax
	end if

51 	continue

c Print out kernel and sky solution ...

c KL
	call write_array(fkern, kernel1, nkxy, nkxy)
c endKL

	if (debug) write(6,'(a)')'Kernel and sky solution ...'
	sum=0.
	sumsq=0.
	do 12 j=1,nkxy
	if (debug) write(6,'(13f10.3)')(kernel1(i,j),i=1,nkxy)
	do 12 i=1,nkxy
	sum=sum+kernel1(i,j)
12	sumsq=sumsq+kernel1(i,j)**2
	write(6,'(a,f7.4,a,f9.2))')'Kernel sum =',sum,
     &    ',  sky offset  =',dsky
	write(6,'(a,f6.3,a)')'Convolution scales noise by',
     &    sqrt(sumsq),' X'
	write(6,'(a)')'   Star   x     y      peak   chisq/nu'//
     &    ' sigma   fwhm0   fwhm1     dx      dy   used?'
	do 55 is=1,nstar
	ans=' '
	if (chisq(is) .gt. chisqrejmax .or. iuse(is) .le. 0) ans='X'
55	write(6,'(3x,i3,2i6,f10.1,6f8.2,4x,a1)')is,nint(xstar(is)),
     &    nint(ystar(is)),z0smax(is)-sky0,chisq(is),sigma(is),
     &    fwhm0(is), fwhm1(is), dx(is), dy(is), ans

c Write to a file the list of stars that were actually used
	open(9, file='psfstars.log', status='unknown', access=
     &	   'sequential')
       write(9,'(a38,a37)')'#    x     y   id      peak   chisq   ',
     &		'sigma   fwhm0   fwhm1      dx      dy'
	do 56 is=1,nstar
	  if (chisq(is) .gt. chisqrejmax .or. iuse(is) .le. 0) goto 56
	  write(9,'(2i6,1x,i4,f10.1,6f8.2)') nint(xstar(is)),
     &		nint(ystar(is)),is,z0smax(is)-sky0,chisq(is),sigma(is),
     &		fwhm0(is),fwhm1(is),dx(is),dy(is)
56	continue
	close(unit=9)

c Compute median FWHM, centroid for each input image just for info. 

	rmdnfwhm0=rmedian1(fwhm0,iuse,nstar,temp)
	rmdnfwhm1=rmedian1(fwhm1,iuse,nstar,temp)
	dx1=rmedian1(dx,iuse,nstar,temp)
	dy1=rmedian1(dy,iuse,nstar,temp)
		
c Finish up - if "doing it", then convolve entire first image with
c kernel, and write both it and the difference image to disk.

	if (doit) then
	  write(6,'(a)')'Convolving entire first image with kernel ...'
	  call conv2d (z0, nx, ny, kernel1, nkxy, nkxy, z2)
CCC KL: allow user to choose whether to subtract the first, convolved image, z0,
CCC or the second, reference image, z1, from the other.  Default is subtract
CCC convolved from reference.
	  if (inverse) then
	    do 14 i=1,nx*ny
	    z2(i)=z2(i)+dsky
14	    diff(i)=z2(i)-z1(i)
	  else
	    do 15 i=1,nx*ny
	    z2(i)=z2(i)+dsky
15	    diff(i)=z1(i)-z2(i)
	  end if
CCC endKL
	  write(6,'(a)')'Writing out convolved and difference images ...'
	  call imdele (fnmout, ier)
	  call iraf_write_r4 (fnmout, z2, nx, ny, ier)
	  call imdele (fnmdiff, ier)
	  call iraf_write_r4 (fnmdiff, diff, nx, ny, ier)
	end if

c Finally, update logfile 'psfmatch.log'.

	open (7, file='psfmatch.log', status='unknown', access=
     &    'append')
	write(7,'(2f10.5,2f10.2,4f7.2,3(2x,a))')sum,sqrt(sumsq),
     &    sky1-sky0,dsky,dx1,dy1,rmdnfwhm0,rmdnfwhm1,
     &    fnm0(1:length(fnm0)), fnm1(1:length(fnm1)),
     &    fnmout(1:length(fnmout))
c

99	stop
	end

	subroutine find_isolated (z, nx, ny, zmin0, zmax0, xstar, 
     &    ystar, nstarmax, nstar)
c
c Find relatively isolated bright stars ... 1. Find a local max whose
c central intensity is in the range zmin0-zmax0 (above sky). 2. Reject
c if there's another local max within icheck pixels whose brightness is
c > frac X zmin0 (usually 10% of zmin0). 3. Also reject if there are
c any pixels < (1-frac)*sky (i.e. anomalously low) within icheck pixels
c of local max. Seems to work very robustly. 4. Also reject if star is
c within icheck pixels of edge of image (new!).
c
c Better to use "fin_isolated_new" - see below. The new version uses
c sigma(sky) to compute whether a star is accepted or not.
c
	real z(nx,ny), xstar(nstarmax), ystar(nstarmax)
	data icheck/13/, frac/0.1/
c Get sky level first, do some other checks ...
	call rmnsd_sampled (z, nx*ny, sky, skysig)
	zminplsky=zmin0+sky
	zminfracplsky=frac*zmin0+sky
	zmaxplsky=zmax0+sky
	nstar=0
	if (nx .lt. 2*icheck+3 .or. ny .lt. 2*icheck+3) then
	  write(6,'(a)')'Error: array too small to find stars!'
	  return
	end if
c Start looking (ignore icheck pixels around edge) ...
	do 1 j=icheck+1,ny-icheck
	do 1 i=icheck+1,nx-icheck
c Ignore pixels outside desired range ...
	if (z(i,j) .lt. zminplsky .or. z(i,j) .gt. zmaxplsky) go to 1
c Ignore if not a local max ...
	zmaxsurround=amax1(z(i-1,j-1),z(i,j-1),z(i+1,j-1),z(i-1,j),
     &    z(i+1,j),z(i-1,j+1),z(i,j+1),z(i+1,j+1))
	if (z(i,j) .lt. zmaxsurround) go to 1
c Check within icheck pixels in x and y - another local max more than
c frac*zmin0 above sky? Also check for other really bright (brighter than
c local max) or unusually dim (more than a fraction frac fainter than
c mean sky) pixels.
	do 2 jj=max0(2,j-icheck),min0(ny-1,j+icheck)
	do 2 ii=max0(2,i-icheck),min0(nx-1,i+icheck)
	if (z(ii,jj) .gt. z(i,j) .or. z(ii,jj) .lt. (1.-frac)*sky) 
     &    go to 1
	if (z(ii,jj) .ge. zminfracplsky) then
	  zmaxsurround1=amax1(z(ii-1,jj-1),z(ii,jj-1),z(ii+1,jj-1),
     &      z(ii-1,jj),z(ii+1,jj),z(ii-1,jj+1),z(ii,jj+1),z(ii+1,jj+1))
	  if (z(ii,jj) .gt. zmaxsurround1 .and. (ii .ne. i .or. 
     &      jj .ne. j)) go to 1
	end if
2	continue
c OK, it's a local maximum and it's isolated and in a relatively clean
c region, so write it out ...
	nstar=nstar+1
	if (nstar .gt. nstarmax) stop 'nstar > nstarmax.'
	xstar(nstar)=i
	ystar(nstar)=j
	write(19,'(2i6)')i,j
c End of loop (jump to here if criteria not satisfied) ...
1	continue
	close(19)
	return
	end

	subroutine find_isolated_new (z, nx, ny, zmin0, zmax0, xstar, 
     &    ystar, nstarmax, nstar)
c
c Find relatively isolated bright stars ... 1. Find a local max whose
c central intensity is in the range zmin0-zmax0 (above sky). 2. Reject
c if there's another local max within icheck pixels whose brightness is
c > sky+k*skysig. 3. Also reject if there are any pixels < sky-k*skysig
c (i.e. anomalously low) within icheck pixels of local max. Seems to work 
c very robustly. 4. Also reject if star is within icheck pixels of edge
c of image (new!).
c
c This version works better than "find_isolated" when zmin0 starts to
c get small! 
c
	real z(nx,ny), xstar(nstarmax), ystar(nstarmax)
	data icheck/13/, ksiglo/5/, ksighi/10/
c Get sky level first, do some other checks ...
	call rmnsd_sampled (z, nx*ny, sky, skysig)
	zminplsky=zmin0+sky
	zmaxplsky=zmax0+sky
	nstar=0
	if (nx .lt. 2*icheck+3 .or. ny .lt. 2*icheck+3) then
	  write(6,'(a)')'Error: array too small to find stars!'
	  return
	end if
c Start looking (ignore icheck pixels around edge) ...
	do 1 j=icheck+1,ny-icheck
	do 1 i=icheck+1,nx-icheck
c Ignore pixels outside desired range ...
	if (z(i,j) .lt. zminplsky .or. z(i,j) .gt. zmaxplsky) go to 1
c Ignore if not a local max ...
	zmaxsurround=amax1(z(i-1,j-1),z(i,j-1),z(i+1,j-1),z(i-1,j),
     &    z(i+1,j),z(i-1,j+1),z(i,j+1),z(i+1,j+1))
	if (z(i,j) .lt. zmaxsurround) go to 1
c Check within icheck pixels in x and y - another local max more than
c frac*zmin0 above sky? Also check for other really bright (brighter than
c local max) or unusually dim (more than a fraction frac fainter than
c mean sky) pixels.
	do 2 jj=max0(2,j-icheck),min0(ny-1,j+icheck)
	do 2 ii=max0(2,i-icheck),min0(nx-1,i+icheck)
	if (z(ii,jj) .gt. z(i,j) .or. z(ii,jj) .lt. sky-ksiglo*skysig)
     &    go to 1
	if (z(ii,jj) .ge. sky+ksighi*skysig) then
	  zmaxsurround1=amax1(z(ii-1,jj-1),z(ii,jj-1),z(ii+1,jj-1),
     &      z(ii-1,jj),z(ii+1,jj),z(ii-1,jj+1),z(ii,jj+1),z(ii+1,jj+1))
	  if (z(ii,jj) .gt. zmaxsurround1 .and. (ii .ne. i .or. 
     &      jj .ne. j)) go to 1
	end if
2	continue
c OK, it's a local maximum and it's isolated and in a relatively clean
c region, so write it out ...
	nstar=nstar+1
	if (nstar .gt. nstarmax) stop 'nstar > nstarmax.'
	xstar(nstar)=i
	ystar(nstar)=j
	write(19,'(2i6)')i,j
c End of loop (jump to here if criteria not satisfied) ...
1	continue
	close(19)
	return
	end


	subroutine matrixmult8 (a, b, c, n)
	real*8 a(n,n), b(n), c(n)
	do 2 j=1,n
	c(j)=0.d0
	do 2 i=1,n
2	c(j)=c(j)+a(i,j)*b(i)
	return
	end


	subroutine conv2d (z, nx, ny, kernel, nkx, nky, zc)
	real z(nx,ny), kernel(nkx,nky), zc(nx,ny)
	if (mod(nkx,2) .ne. 1 .or. mod(nky,2) .ne. 1) stop
     &    'kernel dimensions must be odd.'
	nkx2=nkx/2
	nky2=nky/2
	do 1 j=1,ny
	do 1 i=1,nx
1	zc(i,j)=0.
	do 2 j0=nky2+1,ny-nky2
	do 2 i0=nkx2+1,nx-nkx2
	r=0.
	do 3 j=1,nky
	idy=j-nky2-1
	do 3 i=1,nkx
	idx=i-nkx2-1
3	r=r+kernel(i,j)*z(i0-idx,j0-idy)
2	zc(i0,j0)=r
	return
	end


	subroutine fitmoffat2d (z, nx, ny, sky, fwhm, xc, yc, toplbl0,
     &    ierr)
c
c Full 2D fit of a circularly symmetric Moffat profile to an image,
c including centroid. Assumes a small subraster containing only 1 object, 
c roughly centred in the subraster. At present returns only the FWHM and 
c centroid, but could be made to return all of the Moffat fcn params.
c ierr=0 - everything OK;  ierr=-1 - couldn't compute FWHM; ierr=-2 - 
c linear equations fail.
c
	parameter (nmax=11000)
	character toplbl0*(*), toplbl*120, chp(10)*4, chNaN*4
	real z(nx,ny), zfit(nmax), wt(nmax), e1(10), e2(10), p(10), 
c Note: xy needed for call to lqf, but its values don't have to be loaded!
c See function aux_moffat2d for explanation!
     &    xy(nmax)
	external aux_moffat2d, aux_amoffat2d
	logical debug, sym
	common /fit/ debug, nxfit, nyfit, sym
	equivalence (chp(1),p(1)), (chNan,rNan)
	data rNaN/'7FC00000'x/
	if (nx*ny .gt. nmax) stop 'fitmoffat: too many data points.'
	toplbl=toplbl0
	nxfit=nx
	nyfit=ny
	call centroid (z, nx, ny, sky, 11., float(nx/2+1), 
     &	  float(ny/2+1), xc, yc)
	if (debug) write(6,'(11i7)')((nint(z(i,j)-sky),i=8,18),j=8,18)
	fwhm=fwhmcrude (z, nx, ny, xc, yc, sky, ierr)
	if (ierr .lt. 0) go to 90
ccc	  if (debug) call pgadvance
ccc	  return
ccc	end if
	if (debug) write(6,'(3x,a,2f6.2,a,f6.2)')
     &     toplbl(1:length(toplbl))//'   init (xc,yc)=',
     &     xc,yc,'   fwhm=',fwhm

c Definition of p (free param of Moffat profile).  Symmetrical.
c   p(1) = sky,  p(2) = Ic,  p(3) = alpha,  p(4) = beta, 
c   p(5) = xc,  p(6) = yc
c Definition of p (free param of Moffat profile).  Asymmetrical.
c   p(1) = sky,  p(2) = Ic,  p(3) = alpha_x,  p(4) = alpha_y, 
c   p(5) = beta,  p(6) = xc,  p(7) = yc
	if (sym) then
	  p(1)=sky
 	  p(2)=z(nint(xc),nint(yc))-sky
	  p(3)=0.88*fwhm
	  p(4)=-2.5
	  p(5)=xc
	  p(6)=yc
	else
	  p(1)=sky
	  p(2)=z(nint(xc),nint(yc))-sky
	  p(3)=0.88*fwhm
	  p(4)=0.88*fwhm
	  p(5)=-2.5
	  p(6)=xc
	  p(7)=yc
	endif

c Do the fit.
	if (sym) then
	  write(*,*) 'pin: ',p(1),p(2),p(3),p(4),p(5),p(6)
	  call lqf (xy, z, zfit, wt, e1, e2, p, 0.0, nx*ny, 6, 20, nd,
     &    	1.e-7, aux_moffat2d)
	  write(*,*) 'pout: ',p(1),p(2),p(3),p(4),p(5),p(6)
       else
	  write(*,*) 'pin: ',p(1),p(2),p(3),p(4),p(5),p(6),p(7)
	  call lqf (xy, z, zfit, wt, e1, e2, p, 0.0, nx*ny, 7, 20, nd,
     &    	1.e-7, aux_amoffat2d)
	  write(*,*) 'pout: ',p(1),p(2),p(3),p(4),p(5),p(6),p(7)
	endif
	ierr=0
	if (nd .eq. 0) ierr=-2
c Some checks here for valid fit - otherwise could get screwed up by
c bad pixels, cosmic rays, etc ... Note that the checks for the NaN
c condition have to be done on character variables.
	if (sym) then
	  if (p(3) .le. 0. .or. p(4) .ge. 0. .or. p(5) .lt. 0. .or.
     &      p(5) .gt. nx .or. p(6) .lt. 0. .or. p(6) .gt. ny .or.
     &      chp(1) .eq. chNaN .or. chp(2) .eq. chNaN .or. chp(3) .eq.
     &      chNaN .or. chp(4) .eq. chNaN .or. chp(5) .eq. chNaN .or. 
     &      chp(6) .eq. chNaN) ierr=-3
	else
	  if (p(3) .le. 0. .or. p(4) .le. 0. .or. p(5) .ge. 0. .or. 
     &	     p(6) .lt. 0. .or. p(6) .gt. nx .or. p(7) .lt. 0. .or. 
     &      p(7) .gt. ny .or. chp(1) .eq. chNaN .or. chp(2) .eq. chNaN
     &      .or. chp(3) .eq. chNaN .or. chp(4) .eq. chNaN .or. 
     &      chp(5) .eq. chNaN .or. chp(6) .eq. chNaN .or. 
     &      chp(7) .eq. chNaN) ierr=-3
	endif

90	if (ierr .ge. 0) then
         if (sym) then
	    fwhm=2.*p(3)*sqrt(amax1(0.,0.5**(1./p(4))-1.))
	    xc=p(5)
	    yc=p(6)
	  else
	    fwhm=2.*((p(3)+p(4))/2.)*sqrt(amax1(0.,0.5**(1./p(5))-1.))
	    xc=p(6)
	    yc=p(7)
	  endif
	else
	  fwhm=0.
	  xc=1.
	  yc=1.
	endif
c Debug stuff here ... Can only plot if no error condition - oterwise
c a serious chance of screwing things up!
	if (debug) write(6,'(3x,a,2f6.2,a,f6.2)')
     &     toplbl(1:length(toplbl))//'  final (xc,yc)=',
     &     xc,yc,'   fwhm=',fwhm
	if (debug) write(6,*)(p(k),k=1,6),ierr
	if (debug) then
	  call pgsch(2.0)
	  call pgenv (0., 13., sky-100., z(nint(xc),nint(yc))+100.,
     &      0, 0)
	  call pglabel ('R [pix]', 'I [du]', toplbl)
	  if (ierr .eq. 0) then
	    do 5 j=1,ny
	    do 5 i=1,nx
	    r=sqrt((i-p(5))**2+(j-p(6))**2)
	    call pgpoint (1, r, z(i,j), 20)
5	    call pgpoint (1, r, zfit((j-1)*nx+i), 2)
	  end if
	end if
	return
	end

	function aux_moffat2d (p, d, x, l)
	real p(6), d(6)
	logical debug
	common /fit/ debug, nxfit, nyfit, sym
	j=(l-1)/nxfit+1
	i=l-(j-1)*nxfit	
	r=sqrt((i-p(5))**2+(j-p(6))**2)
	aux_moffat2d=p(1)+p(2)*(1.+(r/p(3))**2)**p(4)
	d(1)=1.
	d(2)=(aux_moffat2d-p(1))/p(2)
	d(3)=-2.*p(2)*p(4)*((1.+(r/p(3))**2)**(p(4)-1.))*(r**2)/(p(3)**3)
	d(4)=(aux_moffat2d-p(1))*alog(1.+(r/p(3))**2)
	d(5)=-2.*p(2)*p(4)*((1.+(r/p(3))**2)**(p(4)-1.))*(i-p(5))/p(3)**2
	d(6)=-2.*p(2)*p(4)*((1.+(r/p(3))**2)**(p(4)-1.))*(j-p(6))/p(3)**2
	return
	end

	function aux_amoffat2d (p, d, x, l)
	real p(7), d(7), brackett
	logical debug
	common /fit/ debug, nxfit, nyfit, sym
	j=(l-1)/nxfit+1
	i=l-(j-1)*nxfit
	brackett = 1 + ((i-p(6))/p(3))**2 + ((j-p(7))/p(4))**2
	aux_amoffat2d = p(1) + p(2)*(brackett)**p(5)
	d(1) = 1
	d(2) = (aux_amoffat2d-p(1))/p(2)
	d(3) = -2.*p(2)*p(5)*(brackett**(p(5)-1))*((i-p(6))**2)/(p(3)**3)
	d(4) = -2.*p(2)*p(5)*(brackett**(p(5)-1))*((j-p(7))**2)/(p(4)**3)
	d(5) = (aux_amoffat2d-p(1))*alog(brackett)
	d(6) = -2.*p(2)*p(5)*(brackett**(p(5)-1))*(i-p(6))/(p(3)**2)
	d(7) = -2.*p(2)*p(5)*(brackett**(p(5)-1))*(i-p(7))/(p(4)**2)
	return
	end

	subroutine centroid (z, nx, ny, sky, rad, xc0, yc0, xc1, yc1)
c
c Compute centroid within aperture of radius rad with centre (xc0,yc0).
c New centroid is (xc1,yc1). (From photkronnew.f, but modified to handle
c non-zero sky, the level of which is an input parameter.)
c
	real z(nx,ny)
	radsq=rad*rad
	sum=0.
	sumx=0.
	sumy=0.
	do 1 j=1,ny
	do 1 i=1,nx
	rsq=(i-xc0)**2+(j-yc0)**2
	if (rsq .le. radsq) then
	  sum=sum+(z(i,j)-sky)
	  sumx=sumx+(z(i,j)-sky)*i
	  sumy=sumy+(z(i,j)-sky)*j
	end if
1	continue
	xc1=sumx/sum
	yc1=sumy/sum
	return
	end

	function fwhmcrude (z, nx, ny, xc, yc, sky, istat)
c
c Crude estimate of FWHM. Assumes a small subraster roughly centered
c on the object. Use fitmoffat for a better estimate. Needs a sky
c estimate (input) and centroid (input) to work. istat < 0 means
c a problem, in which case fwhmcrude=0. is returned.
c
	real z(nx,ny)
	ixc=nint(xc)
	iyc=nint(yc)
c Check for bogus centroid! Occasionally happens ...
	if (ixc .le. 1 .or. ixc .gt. nx .or. iyc .lt. 1 .or.
     &    iyc .gt. ny) go to 9
	zfwhm=sky+(z(ixc,iyc)-sky)/2.
	i1=0
	do 1 i=ixc-1,1,-1
	if (z(i,iyc) .le. zfwhm) then
	  x1=i+(zfwhm-z(i,iyc))/(z(i+1,iyc)-z(i,iyc))
	  go to 2
	end if
1	continue
	go to 9
2	do 3 i=ixc+1,nx
	if (z(i,iyc) .le. zfwhm) then
	  x2=i-(zfwhm-z(i,iyc))/(z(i-1,iyc)-z(i,iyc))
	  fwhmcrude=abs(x2-x1)
	  istat=0
	  return
	end if
3	continue
c One cause of the following error is that you forgot to use *aligned*
c images! But, it does happen (especially near edges) and is not fatal -
c the star is just ignored.
9	write(6,'(a)')'Error: star not contained in subraster.'//
     &    ' Continuing ...'
	write(6,*)nx,ny,xc,yc,sky,zfwhm
	if (ixc-7 .ge. 1 .and. ixc+7 .le. nx .and. iyc-7 .ge. 1 .and.
     &    iyc+7 .le. ny) then
	  do 10 j=iyc+7,iyc-7,-1
10	  write(6,'(15i6)')(nint(z(i,j)-sky),i=ixc-7,ixc+7)
ccccc	  write(6,'(a,$)')'Hit return to continue ... '
ccccc	  read(5,*,end=99)
	end if
	fwhmcrude=0.
	istat=-1
	return
99	stop
	end

	function rmedian (x, n, temp)
c 
c Compute accurate median using sorting. If temp is the same array
c as x, this is OK, but then of course x is destroyed.
c
	real x(1), temp(1)
	call mv4 (x, temp, n)
	call qsort (temp, n, 10)
	rmedian = 0.5*(temp((n+1)/2)+temp(n/2+1))
	return
	end

	function rmedian1 (x, iuse, n, temp)
c 
c Compute accurate median using sorting. This version uses a weight
c iuse(i) to determine whether a value is used or not. If temp is the 
c same array as x, this is OK, but then of course x is destroyed.
c
	real x(1), temp(1)
	integer iuse(1)
	nn=0
	do 1 i=1,n
	if (iuse(i) .gt. 0) then
	  nn=nn+1
	  temp(nn)=x(i)
	end if
1	continue
	call qsort (temp, nn, 10)
	rmedian1 = 0.5*(temp((nn+1)/2)+temp(nn/2+1))
	return
	end

	subroutine fitmoffat1d (z, nx, ny, sky, fwhm, xc, yc)
c
c Fit a circularly symmetric Moffat profile to a radial profile. Assumes a
c small subraster containing only 1 object, roughly centred in the
c subraster. At present returns only the FWHM and centroid, but could 
c be made to return all of the Moffat fcn params.
c
	parameter (nmax=10000)
	real z(nx,ny), r(nmax), zr(nmax), zrfit(nmax), wt(nmax),  
     &    e1(10), e2(10), p(10)
	external aux_moffat1d
	logical debug
	common /fit/ debug, nxfit, nyfit, sym
	if (nx*ny .gt. nmax) stop 'fitmoffat: too many data points.'
	call centroid (z, nx, ny, sky, 11., float(nx/2+1), 
     &	  float(ny/2+1), xc, yc)
	if (debug) write(6,'(11i7)')((nint(z(i,j)),i=8,18),j=8,18)
	fwhm=fwhmcrude (z, nx, ny, xc, yc, sky, ierr)
	if (debug) write(6,*)xc, yc, fwhm
	nr=0
	do 1 j=1,ny
	do 1 i=1,nx
	nr=nr+1
	r(nr)=sqrt((i-xc)**2+(j-yc)**2)
1	zr(nr)=z(i,j)
	p(1)=sky
 	p(2)=z(nint(xc),nint(yc))-sky
	p(3)=0.88*fwhm
	p(4)=-2.5
	call lqf (r, zr, zrfit, wt, e1, e2, p, 0.0, nr, 4, 20, nd,
     &    1.e-7, aux_moffat1d)
	fwhm=2.*p(3)*sqrt(0.5**(1./p(4))-1.)
	if (debug) call pgenv (0., 13., sky-100., z(nint(xc),
     &    nint(yc))+100., 0, 0)
	if (debug) call pgpoint (nr, r, zr, 20)
	if (debug) call pgpoint (nr, r, zrfit, 2)
	return
	end

	function aux_moffat1d (p, d, r, l)
	real p(4), d(4)
	aux_moffat1d=p(1)+p(2)*(1.+(r/p(3))**2)**p(4)
	d(1)=1.
	d(2)=(aux_moffat1d-p(1))/p(2)
	d(3)=-2.*p(2)*p(4)*((1.+(r/p(3))**2)**(p(4)-1.))*(r**2)/(p(3)**3)
	d(4)=(aux_moffat1d-p(1))*alog(1.+(r/p(3))**2)
	return
	end
