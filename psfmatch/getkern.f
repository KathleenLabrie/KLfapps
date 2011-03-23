	subroutine getkern (z0, zc, nx, ny, kernel, nkx, nky, offset,
     &    iop)
c
c Find a kernel that, when convolved with image z0, optimally (in a least 
c squares sense) transforms it to image zc. This will account for PSF
c changes, small shifts, and photometric matching, and, now, even sky level
c offsets!
c
c iop=1-7 - i.e. 3 bits can be set.
c   1 bit - initialize (not set = don't initialize)
c   2 bit - include/process star (not set = don't)
c   4 bit - compute (not set = don't)
c e.g. 7 = initialize, process star, and compute (i.e. just 1 star used)
c e.g. 2 = just add a new star

c Note: you can now compute the solution each time a new star is added,
c without reinitializing.
c 
c The following nmax (=100000) will allow for a 16x16 kernel - i.e. 17**4
c array elements. The extra element is for the sky determination.
c
	parameter (nmax=9200000)
c
	real z0(nx,ny), zc(nx,ny), kernel(1)
	real*8 a, b, aa(nmax), bb(nmax)
	real*8 junk
	common /getkernarrays/ a(nmax), b(nmax)

	junk=(nkx*nky)**2
	write(*,*) nkx,nky,junk,nmax
	if ((nkx*nky)**2 .gt. nmax) stop 'nmax too small.'
	write(*,*) 'getkern: here 1'
	if (mod(nkx,2) .ne. 1 .or. mod(nky,2) .ne. 1) stop
     &    'kernel dimensions must be odd.'

	iinit=iand(iop,1)
	iproc=iand(iop,2)
	icomp=iand(iop,4)
	
	nkx2=nkx/2
	nky2=nky/2
	ndim=nkx*nky+1

c ---------------------- 1. initialize if requested ---------------------

	if (iinit .eq. 1) then
	  do 120 i=1,nmax
	  a(i)=0.d0
120	  b(i)=0.d0
	end if

c --------- 2. process (i.e. include another star) if requested ---------

	if (iproc .ne. 0) then
c Loop over first nkx*nky lines of matrix ...
	  irow=0
	  do 1 l=-nky2,nky2
	  do 1 k=-nkx2,nkx2
	  irow=irow+1
c Make RHS for this line ...
	  do 2 j=nky2+1,ny-nky2
	  do 2 i=nkx2+1,nx-nkx2
2	  b(irow)=b(irow)+zc(i,j)*z0(i-k,j-l)
	  icol=0
c Make first nkx*nky elements of matrix
	  do 3 l1=-nky2,nky2
	  do 3 k1=-nkx2,nkx2
	  icol=icol+1
	  ia=(irow-1)*ndim+icol
	  do 4 j=nky2+1,ny-nky2
	  do 4 i=nkx2+1,nx-nkx2
4	  a(ia)=a(ia)+z0(i-k1,j-l1)*z0(i-k,j-l)
3	  continue
c Sky term at end of each line ...
	  ia=ia+1
	  do 6 j=nky2+1,ny-nky2
	  do 6 i=nkx2+1,nx-nkx2
6	  a(ia)=a(ia)+z0(i-k,j-l)
c End of loop
1	  continue
c Finally make last line ...
	  irow=irow+1
	  icol=0
	  do 8 l1=-nky2,nky2
	  do 8 k1=-nkx2,nkx2
	  icol=icol+1
	  ia=(irow-1)*ndim+icol
	  do 9 j=nky2+1,ny-nky2
	  do 9 i=nkx2+1,nx-nkx2
9	  a(ia)=a(ia)+z0(i-k1,j-l1)
8	  continue
	  ia=ia+1
	  a(ia)=a(ia)+(ny-2*nky2)*(nx-2*nkx2)
	  do 10 j=nky2+1,ny-nky2
	  do 10 i=nkx2+1,nx-nkx2
10	  b(irow)=b(irow)+zc(i,j)
c Check sizes ...
	  if (ia .ne. ndim**2)stop 'something wrong with dimensions of a!'
	  if (irow .ne. ndim) stop 'something wrong with dimensions of b!'
	end if

c ---------------- 3. Do computation and return results ---------------
c Note that arrays are moved into temporary arrays before solution of 
c linear equations, because dlineq messes up arrays!

	if (icomp .ne. 0) then
	  call mv8 (a, aa, ndim**2)
	  call mv8 (b, bb, ndim)
	  call dlineq (ndim, aa, bb)
	  do 11 i=1,ndim-1
11	  kernel(i) = bb(i)
	  offset = bb(ndim)
	end if

	return
	end
