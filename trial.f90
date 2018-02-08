program trial

use mcnp_random
implicit none

integer :: i
real(8)  :: rx,j
  integer, parameter :: I8 = selected_int_kind(18)
  integer, parameter :: R8 = selected_real_kind(15,307)


                    j=0;i=0
call RN_init_problem( 1234567_I8, 1 )

                do while (i<10000000) 
                    call RN_init_particle( int(i,I8) )
                    rx=rang()
                    if (rx<0.0317223) then
                        j=j
                        
                    else if (rx<0.2034294) then
                        j=j+1
                        
                    else if (rx<0.5396285) then
                        j=j+2
                        
                    else if (rx<0.843598) then
                        j=j+3

                    else if (rx<0.9705439) then
                        j=j+4

                    else if (rx<0.9972232) then
                        j=j+5

                    else if (rx<0.9998554) then
                        j=j+6

                    else
                        j=j+7
                    end if
                i=i+1
                end do
print *, "nubar=", real(j/10000000)

end program trial
