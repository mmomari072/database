program rod 

  use mcnp_random 
  implicit none

  integer, parameter :: I8 = selected_int_kind(18)
  integer, parameter :: R8 = selected_real_kind(15,307)
  real(R8), parameter :: n0 = 100000d+0
  integer, parameter :: nrun = 500
  integer, parameter :: bins= 101
  real(R8), parameter :: length =3.87426d+0
  real(R8), parameter :: binwidth=length/(bins-1) 
  real(R8), parameter :: loc0 = 2d+0
  real(R8), parameter :: nubar = 2.4367d+0
  real(R8), parameter :: error = .001
  integer, parameter :: source=100
  real(R8), parameter :: source_pos=2d+0
  integer  :: i, absorb, leak
  real(R8) :: total, locn, rn, direction, rx, den, Navg, wght, kinf, L2,D, kpath,kpathold, pos, j,xpos_ratio, kact
  integer :: bin0, n,run, particle, nleft, coll, inscat, mult,k
  real(R8) :: sigc, sigs,sigf, siga, sigt, move, omega, PI, Clength, keff, kpathrat, errork,nmax,nmin,nsig,navgn
  real(R8), dimension (0:bins-1) ::fbank, fbanknew, xpos, xposnew, flux, xposmid 
  real(R8), dimension (0:7) :: fis_dist
  !real(R8), dimension (0:nint((bins-1)/25)) :: shannnon

!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
 
    den=19.1d+0           !g/cc
    Navg=6.022d+23       !avagodros number
    wght=238.0289d+0       !molecular weight
    sigc= 3.4d-24       !cm^2
    sigf= 4.19d-24
    sigs= 10.0d-24
    sigc=den*Navg*sigc/(wght)
    sigf=den*Navg*sigf/(wght)
    sigs=den*Navg*sigs/(wght)
    siga= sigc+sigf
    sigt= siga+sigs
    fis_dist=(/ 0.0317223, 0.2034294, 0.5396285, 0.843598, 0.9705439, 0.9972232, 0.9998554, 1.00/)
!-----------------------------------------------------Actual keff-------------------------------------------------------
!Assumes a bare slab
    kinf=nubar*sigf/siga
    D=1.0d+0/(sigs*(1.0d+0-2.0/(3.0*nint(wght))))
    L2=D/siga
    PI=4.D0*DATAN(1.D0)
    Clength=PI*(L2/(kinf-1d+0))**0.5d+0-2.0d+0*2.1312d+0*D
    kact=(nubar*sigf)/(siga+D*(PI/(length+2.1312*2.0*D))**2.0d+0)
    
!////////////////////////////////////////////   Initializes components   /////////////////////////////////////////////////////  

  bin0= loc0/binwidth
  n=n0
  locn= loc0
  fbank(0:bins-1)= 0
  fbank(0:bins-1)= floor(real(n0/bins))
  fbanknew(0:bins-1)=fbank(0:bins-1)
  xposmid(0:bins-1)=0
  run=0

print *, "kact= ", kact, Clength

!////////////////////////////////////////  The beginning of the successive runs   ////////////////////////////////////////////

call RN_init_problem( 1234567_I8, 1 )

do while (run<nrun)                                 !ensures the number of iterations occur
    xpos(0:bins-1)=0
    xpos_ratio=20.0d+0
    kpathrat=10.5d+0
    kpathold = 200d+0

    do while (abs(kpathrat)>error .and. xpos_ratio>error)                !ensures k caluclated is converged before the next run
      leak= 0; absorb=0
      mult=0; coll=0
      i=maxval(maxloc(fbank))
      call RN_init_particle( int(i,I8) )
      do while(sum(fbank)>0)
        if (i>bins) then                          !index remornalization
            i=0
        end if

        if (fbank(i)<1) then                        !No particle present
            i=i+1
        else                                        !particle exist and moves
            locn=i*binwidth
            move=-log(rang())/sigt;                 !how far the particle moves
            omega=sign(1.0d+0,2.0*rang()-1.0) 
            move=move*(2.0d+0*rang()-1d+0)                        !direction of movement

            if (length<locn+move) then              !leaks out right
                leak=leak+1
                fbank(i)=fbank(i)-1
                fbanknew(i)=fbanknew(i)-1
                xpos(nint(locn/binwidth):bins)=xpos(nint(locn/binwidth):bins)+1
            else if(locn+move<0) then               !leaks out left
                leak=leak+1
                fbank(i)=fbank(i)-1
                fbanknew(i)=fbanknew(i)-1
                xpos(0:nint(locn/binwidth))=xpos(0:nint(locn/binwidth))+1 
            else

! track the particles movement for flux calculation
                if (move>0) then
                    xpos(nint(locn/binwidth):nint((locn+move)/binwidth))=xpos(nint(locn/binwidth):nint((locn+move)/binwidth))+1
                else if (move==0) then
                    xpos(nint(locn/binwidth))=xpos(nint(locn/binwidth))
                else
                    xpos(nint((locn+move)/binwidth):nint(locn/binwidth))=xpos(nint((locn+move)/binwidth):nint(locn/binwidth))+1
                end if
!---------------------------------------Determine the type of reaction-----------------------------------------------

                rn=rang();
                if (rn<sigc/sigt) then               !capture
                    absorb =absorb+1
                    fbank(i)=fbank(i)-1
                    fbanknew(i)=fbanknew(i)-1

                else if (rn<(sigc+sigf)/sigt) then   !fission
                    rx=rang();
                    absorb=absorb+1;
                    k=0
                    do while (rx>fis_dist(k))
                        k=k+1
                    end do
                        fbank(i)=fbank(i)-1
                        fbanknew(i)=fbanknew(i)-1   
                        fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+k    
                        mult=mult+k             
                        
                else                                  !scatter
                    fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+1
                    fbank(nint(((move+locn)/binwidth)))=fbank(nint(((move+locn)/binwidth)))+1 
                    fbank(i)=fbank(i)-1
                    fbanknew(i)=fbanknew(i)-1
                    coll=coll+1 
                end if
!--------------------------------------------------------------------------------------------------------------------------
            end if             !movement
        end if                 !particle present
      end do                   !sum(fbank)=0

!-------------------------------------------------Calculated keff---------------------------------------------------------      

       kpath=sum(fbanknew)/n
       kpathrat=(kpath/kpathold-1)*(kpath/kpathold)**run
       kpathold=kpath 
       
!nmax=maxval(fbanknew); nmin=minval(fbanknew); navgn=(nmax-nmin)/2d+0
!nsig=sqrt(sum((fbanknew-navgn)**2d+0)/n0)

!----------------------------------------This Renormalizes the fission bank-----------------------------------------------

  fbanknew=(fbanknew/kpath)  

  do i=0,100,1
    if (fbanknew(i)-floor(fbanknew(i))>rang()) then
      fbanknew(i)=ceiling(fbanknew(i))
    else
      fbanknew(i)=floor(fbanknew(i))
    end if
  end do
  n=sum(fbanknew)

    fbank(0:bins-1)=fbanknew(0:bins-1)
!-------------------------------------------------------------------------------------------------------------------------

!i=0
!do while (i<bins-1)
xpos_ratio=xpos_ratio+(abs(xposmid(i)-xpos(i)/maxval(xpos)))
xposmid=xpos/maxval(xpos)
!end do

end do                         !Keff converged

!---------------------------------------xpos averaging at 1/2 the total runs--------------------------------------------- 
    if (run>nrun/2) then
      xposnew(0:bins-1)=(xpos(0:bins-1)/maxval(xpos)+xposnew(0:bins-1))/2d+0
      keff=(keff+kpath)/2d+0
    else
      xposnew(0:bins-1)=xpos(0:bins-1)/maxval(xpos)
      keff=kpath
    end if

  run=run+1
    
end do                                  !ends the number of runs

!---------------------------------------------------flux calc--------------------------------------------------------
i=1
flux(0:bins-1)=0
pos=0
flux(0)=dsin(PI*(pos+2.1312*D)/(length+2*2.1312*D))
do while(i<bins)
    pos=pos+binwidth
    flux(i)=dsin(PI*(pos+2.1312*D)/(length+2*2.1312*D))
    i=i+1
end do

!----------------------------------------------------output-------------------------------------------------------------
 print *, "keff=", keff, "dkeff[$]=",abs(kpath-keff)/.0065
 print *, maxloc(xposnew), maxloc(flux) 

 OPEN(UNIT=12, FILE="rod.txt", ACTION="write", STATUS="replace")
 write(12,*), xposnew/maxval(xposnew)
 write(12,*), flux 
 close(12)

end program rod
