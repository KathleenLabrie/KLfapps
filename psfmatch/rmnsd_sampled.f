	subroutine rmnsd_sampled (z, n, rmn, rsd)
	parameter (nsample=1000, frac=0.01)
	real z(1), z1(nsample)
	logical print
	data print/.false./
	call get_seed0()
	if (n .lt. nint(nsample/frac)) then
	  call rmnsd_rmv1 (z, n, 2.3, 20, print, rmn, rsd)
	else
	  do 1 i=1,nsample
	  iz=nint(2+ran0(0)*(n-4))
1	  z1(i)=z(iz)
	  call rmnsd_rmv1 (z1, nsample, 2.3, 20, print, rmn, rsd)
	end if
	return
	end
