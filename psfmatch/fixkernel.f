	program fixkernel
	implicit none
c
c  Using the difference image obtain with 'psfmatch', estimate 
c  the kernel representing the residual.  Then fix the kernel found
c  by 'psfmatch'.
c
c  Usage: fixkernel fkern1 fdiff fcoo nkxy naxy fkern2
c
	integer nkxydef, naxydef, nareamax
	parameter (nkxydef=70, naxydef=150, nareamax=20)

	character fkern1*80, fdiff*80, fkern2*80, fcoo*80, str*20
	character imtype*5
	integer nkxy, naxy, nmax, ixoff, iyoff, ier, nx, ny, narea
	integer i, j, narg, ind, ida
	real kern1(nkxydef,nkxydef), kdiff(nkxydef,nkxydef),
     &		kern2(nkxydef, nkxydef),diff(nkxydef*nkxydef)
	real xarea(nareamax), yarea(nareamax), subdiff(naxydef*naxydef),
     &		subflat(naxydef*naxydef), bkg, bkgsig, dbkg
c Functions
	real str2real
	integer length, rsubstr

c  Read command line

	narg = iargc()
	if (narg .eq. 6) then
	  call getarg (1, fkern1)
	  call getarg (2, fdiff)
	  call getarg (3, fcoo)
	  call getarg (4, str)
	  nkxy = nint(str2real (str))
	  call getarg (5, str)
	  naxy = nint(str2real (str))
	  call getarg (6, fkern2)
	else
	  write(6,'(a)')'Usage:  fixkernel fkern1 fdiff fcoo nkxy naxy'//
     &		' fkern2'
	endif

c  Initialize
	nmax = nkxydef*nkxydef


c  Load fdiff
	ind = rsubstr(fdiff,'.')
	imtype = fdiff(ind:)
	if ( imtype(1:length(imtype)) .eq. '.imh') then
	  call iraf_read_r4 (fdiff, diff, nmax, nx, ny, ixoff, iyoff, ier)
	else if ( imtype(1:length(imtype)) .eq. '.fits') then
	  call fits_read_r4 (fdiff, diff, nmax, nx, ny, ixoff, iyoff, ier)
	else
	  write(*,*)'ERROR: Unknown image format.'
	  goto 99
	endif

c  Read areas' coordinates
	call read_coo (fcoo, xarea, yarea, nareamax, narea)

c  Estimate average bkg flux in 'diff'
	call rmnsd_sampled (diff, nx*ny, bkg, bkgsig)

c  Initialize the flat to the average bkg flux
	do 20 j=1,naxy,1 
	do 20 i=1,naxy,1
20	  subflat(naxy*(j-1) + i) = bkg

c  Initialize kernel  (only important argument is the last one)
	call getkern (subdiff, subflat, naxy, naxy, kdiff, nkxy, nkxy, 
     &		dbkg, 1)

c  Loop over the areas
	do 30 ida=1,narea,1
	  call subarray (diff, nx, ny, nint(xarea(ida)), nint(yarea(ida)),
     &		subdiff, naxy, naxy)
     	  call getkern (subdiff, subflat, naxy, naxy, kdiff, nkxy, nkxy,
     &		dbkg, 2)
30	continue

c  Compute kernel
	call getkern (subdiff, subflat, naxy, naxy, kdiff, nkxy, nkxy,
     &		dbkg, 4)

c  Read kernel to be fixed
	call read_array (fkern1, kern1, nkxy, nkxy, 0, 0)

c  Fix kernel
	do 40 j=1,nkxy,1
	do 40 i=1,nkxy,1
40	  kern2(i,j) = kern1(i,j)*kdiff(i,j)

c  Write new kernel
	call write_array (fkern2, kern2, nkxy, nkxy)

c  we're done (or something went wrong)
99	stop
	end
