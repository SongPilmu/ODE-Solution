program test
    
implicit none

    ! Varibale setting
    

    double precision, parameter :: g_Na = 30.0d0
    double precision, parameter :: g_K = 12.0d0
    double precision, parameter :: E_Na = 54.2d0
    double precision, parameter :: E_K = -74.7d0
    double precision, parameter :: E_rest = -68.0d0
    double precision, parameter :: C = 1.0d0
    
    double precision :: an, am, ah, bn, bm, bh, Vm, dt, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16
    double precision :: n, m, h
    double precision :: I_Na, I_K, I_ion, time
    integer :: i, k
    double precision :: i_cur(1000)
    i_cur(:) = 0.0d0
    i_cur(1:80) = 75.0d0
    
    ! File name
    open(unit = 1, file = "Vm_RK4.txt", form="formatted")
    open(unit = 2, file = "Potassium_RK4.txt", form = "formatted")
    open(unit = 3, file = "Sodium_RK4.txt", form = "formatted")

    ! Inital condition
    
    Vm = 0.0d0    
    
    ! Time step
    time = 0.0d0
    dt = 0.01d0

    ! RK_4 algorithm
    
    ! Initial n,m,h
    an = 0.01d0 * (10.0d0 - Vm) / (exp(10.0d0 - Vm/10.0d0) - 1.0d0)
    am = 0.1d0 * (25.0d0 - Vm) / (exp(25.0d0 - Vm/10.0d0) - 1.d0)
    ah = 0.07d0 * exp(-Vm/20.0d0)
    bn = 0.125d0 * exp(-Vm/80.0d0)
    bm = 4.0d0 * exp(-Vm/18.0d0)
    bh = 1.0d0 / (exp(30.0d0 - Vm/10.0d0) + 1.0d0)
    
    n = 0.3d0
    m = 0.065d0
    h = 0.6d0
    write (0, '(a)'), "iter   I_cur    Vm      K     Na"

    ! RK4 method
    
    do i = 1,999
        
        
        an = 0.01d0 * (10.0d0 - Vm) / (exp((10.0d0 - Vm)/10.0d0) - 1.0d0)
        am = 0.1d0 * (25.0d0 - Vm) / (exp((25.0d0 - Vm)/10.0d0) - 1.0d0)
        ah = 0.07d0 * exp(- Vm / 20.0d0)
        bn = 0.125d0 * exp(- Vm / 80.0d0)
        bm = 4.0d0 * exp(- Vm / 18.0d0)
        bh = 1.0d0 / (exp((30.0d0 - Vm) / 10.0d0) + 1.0d0)
        
        ! I_Na = m**3 * g_Na * h * (Vm - E_Na)
        ! I_K = n**4 * g_K * (Vm - E_K)
        ! I_ion = i_cur(i) - I_K - I_Na
        ! 
        ! Vm = Vm + dt * I_ion/C
        
        k1 = dt * f1(i_cur(i),Vm)
        k2 = dt * f1(i_cur(i) + dt/2.0d0, Vm + k1/2.0d0)
        k3 = dt * f1(i_cur(i) + dt/2.0d0, Vm + k2/2.0d0)
        k4 = dt * f1(i_cur(i) + dt, Vm + k3)
        Vm = Vm + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
        
        k5 = dt * f2(Vm,n)
        k6 = dt * f2(Vm + dt/2.0d0, n + k5/2.0d0)
        k7 = dt * f2(Vm + dt/2.0d0, n + k6/2.0d0)
        k8 = dt * f2(Vm + dt, n + k7)
        n = n + (k5 + 2.0d0*k6 + 2.0d0*k7 + k8)/6.0d0
        
        k9 = dt * f3(Vm,m)
        k10 = dt * f3(Vm + dt/2.d0, m + k9/2.0d0)
        k11 = dt * f3(Vm + dt/2.0d0, m + k10/2.d0)
        k12 = dt * f3(Vm + dt, m + k11)
        m = m + (k9 + 2.0d0*k10 + 2.0d0*k11 + k12)/6.0d0
        
        k13 = dt * f4(Vm,h)
        k14 = dt * f4(Vm + dt/2.0d0, h + k13/2.0d0)
        k15 = dt * f4(Vm + dt/2.0d0, h + k14/2.0d0)
        k16 = dt * f4(Vm + dt, h + k15)
        h = h + (k13 + 2.0d0*k14 + 2.0d0*k15 + k16)/6.0d0
        
        time = time + dt
           
        if ( mod(i,10) == 0) then
           write (0, '(f, f, f, f, f)'), time, I_cur(i), Vm-70.0d0, g_K * n**4, g_Na * h * m**3
           write (1,*), time, Vm-70.0d0
           write (2,*), time, g_K * n**4
           write (3,*), time, g_Na * h * m**3
        endif
        
    end do
    
contains    
    
    ! Vm function
    function f1(i_cur, Vm)
        double precision :: f1
        double precision, intent(in) :: i_cur, Vm
        f1 = i_cur - (n**4 * g_K * (Vm-E_K)) - (m**3 * g_Na * h * (Vm - E_Na))
    end function f1
    
    ! n function
    function f2(Vm, n) result(var)
        double precision, intent(in) :: Vm, n
        
        double precision :: var
        
        var = (0.01d0 * (10.0d0 - Vm) / (exp((10.0d0 - Vm)/10.0d0) - 1.0d0)) * (1.0d0 - n) - (0.125d0 * exp(-Vm/80.0d0)) * n
    end function f2
    
    ! m function
    function f3(Vm, m)
        double precision :: m, f3
        double precision, intent(in) :: Vm
        f3 = (0.1d0 * (25.0d0 - Vm) / (exp((25.0d0 - Vm)/10.0d0) - 1.0d0)) * (1.0d0-m) - (4.0d0 * exp(-Vm/18.0d0)) * m
    end function f3
    
    ! h function
    function f4(Vm, h)
        double precision :: h, f4
        double precision, intent(in) :: Vm
        f4 = (0.07d0 * exp(-Vm/20.0d0)) * (1.0d0-h) - (1.0d0 / (exp((30.0d0 - Vm)/10.0d0) + 1.0d0)) * h
    end function f4
    
end program test    