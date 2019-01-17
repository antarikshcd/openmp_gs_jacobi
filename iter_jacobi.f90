module m_iter_jacobi
use m_global_data, only: mk, mkd
use m_frobenius_norm
contains
    subroutine iter_jacobi(N,T_old,T_new,f,k_max,d,X_start, X_end, Y_start, Y_end)
    implicit none
    integer(mk) :: i, j, k, k_max, count0, count1, c_rate
    integer(mk) :: N
    real(mk) :: elap_time, norm, sum, d, delta_X, delta_Y, X_start, X_end, Y_start, Y_end,norm1,norm2, aii
    real(mk),dimension(:,:) :: T_new, T_old, f


    delta_X = (X_end - X_start) / real(N+1, mk)
    delta_Y = (Y_end - Y_start) / real(N+1, mk)
    print*, delta_X
    print*, delta_Y
    norm = 50
    k = 0
    aii = 1.0/4.0

    ! call system clock
    call system_clock(count=count0, count_rate=c_rate)

    do while (norm > d .and. k < k_max)
        sum = 0.0
        do j=2,N+1
            do i=2,N+1
                T_new(i,j) = aii * (T_old(i-1,j) +&
                                        T_old(i+1,j) +&
                                        T_old(i,j-1) +&
                                        T_old(i,j+1) +&
                                        delta_X*delta_Y * f(i,j))
                sum = sum + (T_new(i,j) - T_old(i,j))*(T_new(i,j) - T_old(i,j))
            enddo
        enddo
    k = k + 1
    
    ! calculate the norm
    norm = sum ** 0.5
    print*,'k=',k
    print*,'norm=',norm
    !print*,'d',d !debug
    T_old = T_new
    enddo
    
    ! time
    call system_clock(count=count1)
    elap_time = real(count1 - count0)/real(c_rate)
    print*,'Wall time (fortran):', elap_time,'[s]' 
    print*,'Iteration rate=', k/elap_time,'[iteration/s]'  
    !print*,'delta_X*(N+1)=', delta_X*(N+1) !debug
    !print*,'delta_Y*(N+1)=', delta_Y*(N+1) !debug
    end subroutine iter_jacobi
end module
