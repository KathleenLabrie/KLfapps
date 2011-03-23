	program array2list
	implicit none
c
c  Transform a table into a 3 columns list (col,row,value)
c
c  Usage: array2list fin fout nx ny [xoff yoff] 
c
	real str2real
	integer length

	integer nxmax, nymax
	parameter (nxmax=1000, nymax=1000)

	character fin*80, fout*80, str*20
	integer nx,ny,xoff,yoff
	integer narg,i,j,n,l,k
	real array(nxmax, nymax)

c  Initialize
	nx=1
	ny=1
	xoff=0
	yoff=0

c  Read command line
	narg=iargc()
	if (narg .ge. 4) then
	  call getarg (1,fin)
	  call getarg (2,fout)
	  call getarg (3,str)
	  nx = nint(str2real (str))
	  call getarg (4,str)
	  ny = nint(str2real (str))
	  if (narg .eq. 6) then
	    call getarg (5,str)
	    xoff = nint(str2real (str))
	    call getarg (6,str)
	    yoff = nint(str2real (str))
	  else if (narg .ne. 4) then
	    write(6,'(a)')'Usage: array2list fin fout nx ny [xoff yoff]'
	    stop
	  else
	    continue
	  end if
	else
	  write(6,'(a)')'Usage: array2list fin fout nx ny [xoff yoff]'
	  stop
	end if

c  Check inputs
	if ((nx .gt. nxmax) .or. (ny .gt. nymax)) then
	  write(*,*) 'nx or ny are too large!'
	  stop
	end if

c  Open table, read it, close it
	call read_array (fin, array, nx, ny, xoff, yoff)
	n=nx*ny
	do 1 j=ny,1,-1
	do 1 i=nx,1,-1
	  l=(n-1)/nxmax + 1
	  k=n-(l-1)*nxmax
	  array(i,j)=array(k,l)
1	n=n-1

c  Create list.
	open(10, file=fout, status='new', err=91, access='sequential',
     &		form='formatted')
	do 10 j=1,ny,1
	do 10 i=1,nx,1
10	write(10,*,err=92) i,j,array(i,j)
	close(unit=10, err=93)
	goto 99

c  Error messages
91     write(*,*) 'ERROR: Unable to open file ',fout(1:index(fout,' '))
       write(*,*) 'Aborting...'
       goto 99
92     write(*,*) 'ERROR: Unable to write to file ',
     &							fout(1:index(fout,' '))
       write(*,*) 'Aborting...'
       goto 99
93     write(*,*) 'ERROR: Unable to close file ',
     &							fout(1:index(fout,' '))
       write(*,*) 'Continuing'

99     stop
       end
