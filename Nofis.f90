program rod 

  use mcnp_random 
  implicit none
! This code assumes the geometry starts at 0, with no fission, a
  integer, parameter :: I8 = selected_int_kind(18)
  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: n0 = 100000
  integer, parameter :: nrun = 1
  integer, parameter :: loc0 = 0
  integer, parameter :: length = 10
  real(R8), parameter :: binwidth =0.1
  integer  :: i,absorb,leak
  real(R8) :: total, locn, rn, direction,rx,den,Navg,wght
  integer :: bins, bin0, n,run,pos
  real(R8) :: sigc, sigs,sigf, siga, sigt,nubar,move,omega
  real(R8), dimension (0:99) ::fbank, fbanknew,xpos, catfun
  real(R8), dimension (0:99,0:nrun) :: xposs
  den=19.1 !g/cc
  Navg=6.022e23 !avagodros number
  wght=235.0439    !molecular weight
  sigc= 98.86e-24 !cm^-2
  sigf= 584.994e-24
  sigs= 14.942e-24
  sigc=den*Navg*sigc/(wght)
  sigf=0!den*Navg*sigf/(wght)
  sigs=den*Navg*sigs/(wght)*1000
  siga= sigc+sigf
  sigt= siga+sigs 
  bins= length/binwidth
  bin0= loc0/binwidth
  nubar= 2.4367
  leak=0
  absorb=0
  xpos(0:(bins))=0
  n=n0
  run=0
  locn=loc0
print *, sigc

do while (run<nrun)
  do while(n>0)
    move=-log(rang())/sigt
    omega=2*rang()-1;direction=omega/abs(omega) 
    move=move*omega !movement is either positive or negative
    if (locn==loc0) then
      move=abs(move)
    end if
    if (length<locn+move ) then     !escapes out right
      xpos(nint(locn/binwidth):bins)=xpos(nint(locn/binwidth):bins)+1
      n=n-1
      locn=loc0
    else if (locn+move<0) then     !Escapes out the left
      xpos(0:nint(locn))=xpos(0:nint(locn))+1
      n=n-1
      locn=loc0
    
    else 
      rn=rang()
        if (move>0) then
          if (rn<sigc/sigt) then !capture
            xpos(nint(locn/binwidth):nint((locn+move)/binwidth))=xpos(nint(locn/binwidth):nint((locn+move)/binwidth))+1
            n=n-1
            locn=loc0
            absorb=absorb+1
          else 
            xpos(nint(locn/binwidth):nint((locn+move)/binwidth))=xpos(nint(locn/binwidth):nint((locn+move)/binwidth))+1
            locn=locn+move
          end if
        else
          if (rn<sigc/sigt) then
            xpos(nint((locn+move)/binwidth):nint(locn/binwidth))=xpos(nint((locn+move)/binwidth):nint(locn/binwidth))+1
            n=n-1
            locn=loc0
          else 
            xpos(nint((locn+move)/binwidth):nint(locn/binwidth))=xpos(nint((locn+move)/binwidth):nint(locn/binwidth))+1
            locn=locn+move
          end if
    
        end if  
    end if
  end do
!  xposs(nrun,0:99)=xpos
  run=run+1
end do  
print *, sum(xpos)

!i=0
!do while (i<=bins) 
!catfun(i)=dexp(-sigc*i*binwidth)
!i=i+1
!end do


 OPEN(UNIT=12, FILE="Nofis.txt", ACTION="write", STATUS="replace")
 write(12,*), xpos/maxval(xpos)
 close(12)
end program rod  
