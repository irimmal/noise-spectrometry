module rk4mod

    implicit none
    real :: lendg, magin, m_t0(3), T1, var
    real, parameter :: pi = 3.1415926535897932384626433832795, gamma = 8.8 ! 8.8e10 C/kg
    real, parameter :: h_bar = 1.05457266  ! e-34J*s, 6.58211928e-16 eV*s
    public :: cross_3d, f_t_m, constant_get, rk4_calc, gamma, lendg, magin, m_t0, T1, var

contains

    subroutine rk4_calc(delta_t, y0, ndim, step_num)
        ! 4阶龙格库塔法
        implicit none
        integer :: i
        integer, intent(in) :: ndim, step_num ! 维数，步数
        real, intent(in) :: delta_t ! 步长
        real, intent(in) :: y0(ndim) ! 磁矩初始值
        real, dimension(ndim) :: y, k1, k2, k3, k4
        real, dimension(step_num) :: t_list
        real, dimension(step_num, ndim) :: y_list
        real :: t
        ! 初始值
        t = 0.0
        y = y0
        !B = ?
        !delta_t = ?
        do i = 1, step_num
            open(unit=1, file='rk4.txt', status='unknown')
            write(1, *) t, y
            t_list(i) = t
            y_list(i, :) = y
            k1 = delta_t*f_t_m(t, y)
            k2 = delta_t*f_t_m(t + delta_t/2.0, y + k1/2.0)
            k3 = delta_t*f_t_m(t + delta_t/2.0, y + k2/2.0)
            k4 = delta_t*f_t_m(t + delta_t, y + k3)
            y = y + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
            t = t + delta_t
        end do
        close(unit=1)
    end subroutine rk4_calc

    subroutine constant_get()
        ! 获取常数
        implicit none
        real :: constant(4)
        constant = [lendg, magin, T1, var]
        print *, constant
    end subroutine constant_get

    function cross_3d(a, b)
        ! 3维矢量叉乘
        implicit none
        real, dimension(3) :: a, b, cross_3d
        cross_3d = [a(2)*b(3) - a(3)*b(2), &
                    a(3)*b(1) - a(1)*b(3), &
                    a(1)*b(2) - a(2)*b(1)]
    end function cross_3d

    function f_t_m(t, m) result(f)
        ! 磁矩随时间变化率方程
        implicit none
        real :: t, m(3), f(3), kai(3), m0, B_ag(3), m0_ag(3)
        B_ag = [0.0, 0.0, magin]
        m0 = sqrt(3.0)*lendg*h_bar*gamma/2.0 ! h_bar*gamma = bohr磁子
        m0_ag = [0.0, 0.0, m0]
        kai = random_kai(t, var)
        f = gamma * cross_3d(m, B_ag) - (m - m0_ag)/T1 + kai
    end function f_t_m
    
    function random()
        ! 生成[0,1)的随机数
        implicit none
        integer, save :: mark = 0
        real :: random
        if (mark == 0) then
            call random_seed()
            mark = 1
        end if
        call random_number(random)
    end function random

    function gauss_random(var_in)
        ! box-muller法产生高斯随机数
        implicit none
        real :: var_in, gauss_random
        real :: random1, random2
        random1 = random()
        random2 = random()
        gauss_random = sqrt(-2*var_in*log(random1))*cos(2*pi*random2)
    end function gauss_random

    function random_kai(t,var_in)
        ! 生成白噪声随机矢量
        implicit none
        real :: t, var_in, random_kai(3)
        random_kai(1) = gauss_random(var_in)
        random_kai(2) = gauss_random(var_in)
        random_kai(3) = gauss_random(var_in)
    end function random_kai

end module rk4mod

! test:
!program main
!    use rk4mod
!    implicit none
!    real :: delta_t
!    integer :: ndim, step_num
!    ! 状态量
!    lendg = 2
!    magin = 0.1
!    T1 = 10.0
!    var = 0.1
!    ! 初始值
!    delta_t = 0.001
!    ndim = 3
!    step_num = 1000
!    m_t0 = [0.0, 0.0, 1.0]
!    call rk4_calc(delta_t, m_t0, ndim, step_num)
!end program main
