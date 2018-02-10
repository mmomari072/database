program rod 

  use mcnp_random 
  implicit none

  integer, parameter :: I8 = selected_int_kind(18)
  integer, parameter :: R8 = selected_real_kind(15,307)
  integer, parameter :: n0 = 1000
  integer, parameter :: nrun = 250
  integer, parameter :: bins= 101
  real(R8), parameter :: length =4.3933
  real(R8), parameter :: binwidth=length/(bins-1) 
  real(R8), parameter :: loc0 = 2
  real(R8), parameter :: nubar = 2.4367
  real(R8), parameter :: error = .0001
  integer  :: i, absorb, leak
  real(R8) :: total, locn, rn, direction, rx, den, Navg, wght, kinf, L2,D, kpath,kpathold, pos, j
  integer :: bin0, n,run, particle, nleft, coll, inscat, mult
  real(R8) :: sigc, sigs,sigf, siga, sigt, move, omega, PI, Clength, keff, kpathrat, errork
  real(R8), dimension (0:bins-1) ::fbank, fbanknew, xpos, xposnew, flux 

!///////////////////////////////////////////////////////////////////////////////////////////////   
 
    den=19.1            !g/cc
    Navg=6.022e23       !avagodros number
    wght=238.0289       !molecular weight
    sigc= 3.4e-24       !cm^-2
    sigf= 4.19e-24
    sigs= 10.0e-24
    sigc=den*Navg*sigc/(wght)
    sigf=den*Navg*sigf/(wght)
    sigs=den*Navg*sigs/(wght)
    siga= sigc+sigf
    sigt= siga+sigs

!-----------------------------------------------------Actual keff-------------------------------------------------------
    
    kinf=nubar*sigf/siga
    D=1/(3*sigs*(1-2/(3*nint(wght))))
    L2=D/siga
    PI=4.D0*DATAN(1.D0)
    Clength=PI*(L2/(kinf-1))**0.5-2*2.1312*D
    keff=(nubar*sigf)/(siga+D*(PI/(length+2.1312*2*D))**2)
    
!////////////////////////////////////////////   Initializes components   /////////////////////////////////////////////////////  

  bin0= loc0/binwidth
  n=n0
  locn= loc0
  fbank(0:bins-1)= 0
  fbank(0:bins-1)= floor(real(n0/bins))
  absorb= 0

  coll=0
  fbanknew(0:bins-1)=fbank(0:bins-1)
  run=0

print *, "kact=", keff, (L2)**.5

!////////////////////////////////////////  The beginning of the successive runs   ////////////////////////////////////////////

call RN_init_problem( 1234567_I8, 1 )
do while (run<nrun)                                 !ensures the number of iterations occur
    xpos(0:bins-1)=0
    kpathrat=1.5
    kpathold = 20
    do while (abs(kpathrat-1)>error)                !ensures k caluclated is converged before the next run
      leak= 0; absorb=0
      mult=0; coll=0
      i=maxval(maxloc(fbank))
!      call RN_init_particle( int(i,I8) )
      do while(sum(fbank)>0)
        if (i>bins) then                          !index remornalization
            i=0
        end if

        if (fbank(i)<1) then                        !No particle present
            i=i+1
        else                                        !particle exist and moves
            locn=i*binwidth
            move=-log(rang())/sigt;                 !how far the particle moves
            omega=2*rang()-1;direction=omega/abs(omega) 
            move=move*omega                         !movement is either positive or negative

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
                    rx=rang();absorb=absorb+1
                    if (rx<0.0317223) then
                        fbank(i)=fbank(i)-1
                        fbanknew(i)=fbanknew(i)-1
                    else if (rx<0.2034294) then
                        fbank(i)=fbank(i)-1
                        fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+1
                        fbanknew(i)=fbanknew(i)-1
                        mult=mult+1
                    else if (rx<0.5396285) then
                        fbank(i)=fbank(i)-1
                        fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+2
                        fbanknew(i)=fbanknew(i)-1
                        mult=mult+2
                    else if (rx<0.843598) then
                        fbank(i)=fbank(i)-1
                        fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+3
                        fbanknew(i)=fbanknew(i)-1
                        mult=mult+3
                    else if (rx<0.9705439) then
                        fbank(i)=fbank(i)-1
                        fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+4
                        fbanknew(i)=fbanknew(i)-1
                        mult=mult+4
                    else if (rx<0.9972232) then
                        fbank(i)=fbank(i)-1
                        fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+5
                        fbanknew(i)=fbanknew(i)-1
                        mult=mult+5
                    else if (rx<0.9998554) then
                        fbank(i)=fbank(i)-1
                        fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+6
                        fbanknew(i)=fbanknew(i)-1
                        mult=mult+6
                    else
                        fbank(i)=fbank(i)-1
                        fbanknew(nint(((move+locn)/binwidth)))=fbanknew(nint(((move+locn)/binwidth)))+7
                        fbanknew(i)=fbanknew(i)-1
                        mult=mult+7
                    end if
                        
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

       kpath=sum(fbanknew)/n0
       kpathrat=kpath/kpathold
       kpathold=kpath 
!print *, sum(fbanknew)

!----------------------------------------This Renormalizes the fission bank-----------------------------------------------

fbanknew=nint(fbanknew/kpath)

!Supercritical            
    do while (sum(fbanknew)>n0)
        rx=rang() 
        if (fbanknew(nint(bins*rx))>0) then
            fbanknew(nint(bins*rx))=fbanknew(nint(bins*rx))-1
        end if
    end do
    
!Subcritical    
    do while (sum(fbanknew)<n0)
        rx=rang()        
        fbanknew(nint(bins*rx))=fbanknew(nint(bins*rx))+1
    end do  
    
    fbank(0:bins-1)=fbanknew(0:bins-1)


end do                         !Keff converged
!print *, kpath  
!---------------------------------------xpos averaging at 1/2 the total runs--------------------------------------------- 
    if (run>nrun/2) then
      xposnew(0:bins-1)=(xpos(0:bins-1)/maxval(xpos)+xposnew(0:bins-1))/2

    else
      xposnew(0:bins-1)=xpos(0:bins-1)/maxval(xpos)
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
 print *, "keff=", kpath, "dkeff[$]=",abs(kpath-keff)/.0065
 print *, maxloc(xposnew), maxloc(flux)

 OPEN(UNIT=12, FILE="rod.txt", ACTION="write", STATUS="replace")
 write(12,*), xposnew/maxval(xposnew)
 write(12,*), flux 
 close(12)

end program rod
