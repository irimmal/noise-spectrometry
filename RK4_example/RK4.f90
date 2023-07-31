!  FILE: RK4.f90
module RK4_func
    public :: K1, K2, K3, K4, RK4, calculate_RK4
    ! public :: RK4, calculate_RK4
    ! private :: K1, K2, 
    ! integer :: dim=2

    contains

    function K1(psi_n, H, dt)
        complex (kind=4) :: psi_n(2), K1(2), H(2,2)
        real :: dt
        K1 = 1/(0,1) * matmul(H, psi_n)
    end function K1

    function K2(psi_n, H, dt)
        complex (kind=4) :: psi_n(2), K2(2), H(2,2)
        real (kind=4) :: dt
        K2 = 1/(0,1) * matmul(H, psi_n + 0.5 * dt * K1(psi_n, H, dt))
    end function K2

    function K3(psi_n, H, dt)
        complex (kind=4) :: psi_n(2), K3(2), H(2,2)
        real (kind=4) :: dt
        K3 = 1/(0,1) * matmul(H, psi_n + 0.5 * dt * K2(psi_n, H, dt))
    end function K3

    function K4(psi_n, H, dt)
        complex (kind=4) :: psi_n(2), K4(2), H(2,2)
        real (kind=4) :: dt
        K4 = 1/(0,1) * matmul(H, psi_n + dt * K3(psi_n, H, dt))
    end function K4

    function RK4(psi_n, H, dt) result(psi_n_1)
        complex (kind=4) :: psi_n(2), psi_n_1(2), H(2,2)
        real (kind=4) :: dt
        psi_n_1 = psi_n + dt/6*(K1(psi_n, H, dt) + 2*K2(psi_n, H, dt) + 2*K3(psi_n, H, dt) + K4(psi_n, H, dt))
    end function RK4

    ! final calcualation
    ! subroutine calculate_RK4(dt, N, psi_0, H, psi)
    subroutine calculate_RK4(psi, psi_0, H, dt, N)
        implicit none
        real (kind=4), intent(in) :: dt
        integer, intent(in) :: N
        complex (kind=4), intent(in) :: psi_0(2), H(2,2)
    
        integer :: i ! declaration must precede all codes
        complex (kind=4), intent(inout) :: psi(2, N)

        psi = 0

        ! print *, psi_0
        psi(:,1) = psi_0
    
        do i = 1, N-1
            psi(:, i+1) = RK4(psi(:, i), H, dt)
        end do
    
        ! print *, psi ! WHY does this fail???
    
    end subroutine calculate_RK4
    
end module RK4_func

! ! to test the program in fortran
! program name
!     use RK4_func
!     implicit none
!     complex (kind=4) :: psi(2, 10), H(2,2)
!     ! integer :: N

!     ! external calculate_RK4
!     psi = 0
!     ! N = 10
!     H = reshape((/(1,0), (0,1), (0,-1) ,(1,0)/), (/2,2/))
!     call calculate_RK4(1., 10, [(1,0), (0,0)], H, psi)
!     ! print *, psi

! end program name

! end FILE RK4.f90
