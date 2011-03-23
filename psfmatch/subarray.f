	subroutine subarray (z, nx, ny, ixc, iyc, s, nxs, nys)
	real z(nx,ny), s(nxs,nys)
c Initialize (in case of out-of-bounds problem) ...
	do 1 j=1,nys
	do 1 i=1,nxs
1	s(i,j)=0.
c OK do it ...
	j1=iyc-nys/2
	j2=j1+nys-1
	i1=ixc-nxs/2
	i2=i1+nxs-1
	do 5 j=max0(j1,1),min0(j2,ny)
	jj=j-j1+1
	do 5 i=max0(i1,1),min0(i2,nx)
	ii=i-i1+1
5	s(ii,jj)=z(i,j)
	return
	end
