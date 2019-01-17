program main
use m_input      
use m_global_data 
use m_alloc       
use m_init_test   
use m_init        
use m_frobenius_norm
use m_iter_jacobi
use m_iter_gauss_siedel
use m_output
implicit none

!global_data
print*,'Pi = ',Pi

!input
call input(N,X_start, X_end, Y_start, Y_end,k_max,d)
print*,'N = ',N
print*,'X_start = ',X_start
print*,'X_end = ',X_end
print*,'Y_start = ',Y_start
print*,'Y_end = ',Y_end
print*,'k_max = ',k_max
print*,'d = ',d

!alloc
call alloc(N,T_old,info)
!print*,'T_old is allocated',info !debug
call alloc(N,T_new,info)
!print*,'T_new is allocated',info !debug
call alloc(N,f,info)
!print*,'f is allocated',info !debug
call alloc(N,f_theo,info)
!print*,'f is allocated',info !debug

!init
!call init(N,T_old,T_new,f,X_start,X_end,Y_start,Y_end,A_X,B_X,A_Y,B_Y)
call init_test(N,T_old,T_new,f,X_start,X_end,Y_start,Y_end,A_X,B_X,A_Y,B_Y)

!iter
call iter_jacobi(N,T_old,T_new,f,k_max,d,X_start, X_end, Y_start, Y_end)
!call iter_gauss_siedel(N,T_old,T_new,f,k_max,d,X_start, X_end, Y_start, Y_end)

!result
call output(N,T_old,T_new,f,f_theo,A_X,B_X,A_Y,B_Y)

!dealloc !add info and print
deallocate(T_old,STAT=info)
!print*,'T_old is deallocated',info !debug
deallocate(T_new,STAT=info)
!print*,'T_new is deallocated',info !debug
deallocate(f,STAT=info)
!print*,'f is deallocated',info !debug
deallocate(f_theo,STAT=info)
!print*,'f_theo is deallocated',info !debug

end program main
