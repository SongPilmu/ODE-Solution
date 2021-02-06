program test
    
Implicit none
    
    REAL :: an, am, ah, bn, bm, bh, Vm, Delta_t
    REAL :: g_Na, g_K, E_Na, E_K, E_rest, C, n, m, h
    REAL :: I_Na, I_K, I_ion, time
    integer :: k
    double precision :: I(10000)

I(:) = 0.0d0
I(1:80) = 75.0d0

! file name
open(unit = 1, file = "Vm.txt", form="formatted")
open(unit = 2, file = "Potassium.txt", form = "formatted")
open(unit = 3, file = "Sodium.txt", form = "formatted")

! inital condition
g_Na = 30.0d0
g_K = 12.0d0
E_Na = 54.2d0
E_K = -74.7d0
E_rest = -68.0d0
C = 1.0d0
Vm = 0.0d0

! Time step
time = 0.0d0
Delta_t = 0.01d0

! n,m,h initial
   

    an = 0.01d0 * (10.0d0 - Vm) / (exp(10.0d0 - Vm/10.0d0) - 1.0d0)
    am = 0.1d0 * (25.0d0 - Vm) / (exp(25.0d0 - Vm/10.0d0) - 1.d0)
    ah = 0.07d0 * exp(-Vm/20.0d0)
    bn = 0.125d0 * exp(-Vm/80.0d0)
    bm = 4.0d0 * exp(-Vm/18.0d0)
    bh = 1.0d0 / (exp(30.0d0 - Vm/10.0d0) + 1.0d0)
    
    n = 0.3d0
    m = 0.065d0
    h = 0.6d0
    write (0, '(a)'), "iter   Vm      K     Na"
    
! Euler forward method
do k = 1,999 ! Total time = 10s, Time step = 0.001
    
    an = 0.01d0 * (10.0d0 - Vm) / (exp((10.0d0 - Vm)/10.0d0) - 1.0d0)
    am = 0.1d0 * (25.0d0 - Vm) / (exp((25.0d0 - Vm)/10.0d0) - 1.0d0)
    ah = 0.07d0 * exp(-Vm/20.0d0)
    bn = 0.125d0 * exp(-Vm/80.0d0)
    bm = 4.0d0 * exp(-Vm/18.0d0)
    bh = 1.0d0 / (exp((30.0d0 - Vm)/10.0d0) + 1.0d0)
    
    I_Na = m**3 * g_Na * h * (Vm - E_Na)
    I_K = n**4 * g_K * (Vm-E_K)
    I_ion = I(k) - I_K - I_Na

    Vm = Vm + Delta_t * I_ion/C
    n = n + Delta_t * (an * (1.0d0 - n) - bn * n)
    m = m + Delta_t * (am * (1.0d0-m) - bm * m)
    h = h + Delta_t * (ah * (1.0d0-h) - bh * h)
    time = time + Delta_t
        
    if( mod(k,10)==0 ) then
        write (0, '(f, f, f, f)'), time, Vm-70.0d0, g_K * n**4, g_Na * h * m**3
        write (1,*), time, Vm-70.0d0
        write (2,*), time, g_K * n**4
        write (3,*), time, g_Na * h * m**3
    endif
        
end do



end program test