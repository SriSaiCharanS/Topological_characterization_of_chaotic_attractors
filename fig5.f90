program elnino
implicit none

integer :: i, i1
real*8 :: h
real*8 :: k11,k12,k13,k21,k22,k23,k31,k32,k33
real*8 :: k41,k42,k43,k51,k52,k53,k61,k62,k63
real*8 :: mu,Tr0,z0,alpha,p1,r,Tr,hs,hp,gt,eps
real*8 :: h10,T10,T20,h1,T1,T2
real*8 :: t0
real*8 :: y1,y2,y3

open(1,file='enso_bifur_last12M.dat')

do i1 = 0,499

    eps = 0.09780d0 + (1.0d-3/499.0d0)*i1
    write(*,*) i1

    ! Parameters
    mu    = 0.0026d0
    Tr0   = 16.0d0
    z0    = 75.0d0
    alpha = 1.0d0/180.0d0
    p1    = 22.0d0
    r     = 1.0d0/400.0d0
    Tr    = 29.5d0
    hs    = 62.0d0
    hp    = 100.0d0
    gt    = 1.3d0

    ! Initial conditions
    h10 = 70.01d0
    T10 = 20.01d0
    T20 = 10.198d0

    h  = 0.01d0
    t0 = 0.0d0

    y1 = T20
    y2 = T20
    y3 = T20

    do i = 1,120000000

        ! ---- RK45 stages ----

        k11 = h*(-r*(h10+(p1*(T20-T10))/2.d0))
        k12 = h*(-alpha*(T10-Tr)-eps*mu*(T20-T10)**2)
        k13 = h*(-alpha*(T20-Tr)+gt*mu*(T20-T10)*(T20-Tr + &
              (Tr-Tr0)*(1.d0-tanh((h10+p1*(T20-T10)+hp-z0)/hs))/2.d0))

        h1 = h10 + k11/4.d0
        T1 = T10 + k12/4.d0
        T2 = T20 + k13/4.d0

        k21 = h*(-r*(h1+(p1*(T2-T1))/2.d0))
        k22 = h*(-alpha*(T1-Tr)-eps*mu*(T2-T1)**2)
        k23 = h*(-alpha*(T2-Tr)+gt*mu*(T2-T1)*(T2-Tr + &
              (Tr-Tr0)*(1.d0-tanh((h1+p1*(T2-T1)+hp-z0)/hs))/2.d0))

        h1 = h10 + (3.d0*k11/32.d0) + (9.d0*k21/32.d0)
        T1 = T10 + (3.d0*k12/32.d0) + (9.d0*k22/32.d0)
        T2 = T20 + (3.d0*k13/32.d0) + (9.d0*k23/32.d0)

        k31 = h*(-r*(h1+(p1*(T2-T1))/2.d0))
        k32 = h*(-alpha*(T1-Tr)-eps*mu*(T2-T1)**2)
        k33 = h*(-alpha*(T2-Tr)+gt*mu*(T2-T1)*(T2-Tr + &
              (Tr-Tr0)*(1.d0-tanh((h1+p1*(T2-T1)+hp-z0)/hs))/2.d0))

        h1 = h10 + (1932.d0*k11/2197.d0) - (7200.d0*k21/2197.d0) + &
              (7296.d0*k31/2197.d0)
        T1 = T10 + (1932.d0*k12/2197.d0) - (7200.d0*k22/2197.d0) + &
              (7296.d0*k32/2197.d0)
        T2 = T20 + (1932.d0*k13/2197.d0) - (7200.d0*k23/2197.d0) + &
              (7296.d0*k33/2197.d0)

        k41 = h*(-r*(h1+(p1*(T2-T1))/2.d0))
        k42 = h*(-alpha*(T1-Tr)-eps*mu*(T2-T1)**2)
        k43 = h*(-alpha*(T2-Tr)+gt*mu*(T2-T1)*(T2-Tr + &
              (Tr-Tr0)*(1.d0-tanh((h1+p1*(T2-T1)+hp-z0)/hs))/2.d0))

        h1 = h10 + (439.d0*k11/216.d0) - (8.d0*k21) + &
              (3680.d0*k31/513.d0) - (845.d0*k41/4104.d0)
        T1 = T10 + (439.d0*k12/216.d0) - (8.d0*k22) + &
              (3680.d0*k32/513.d0) - (845.d0*k42/4104.d0)
        T2 = T20 + (439.d0*k13/216.d0) - (8.d0*k23) + &
              (3680.d0*k33/513.d0) - (845.d0*k43/4104.d0)

        k51 = h*(-r*(h1+(p1*(T2-T1))/2.d0))
        k52 = h*(-alpha*(T1-Tr)-eps*mu*(T2-T1)**2)
        k53 = h*(-alpha*(T2-Tr)+gt*mu*(T2-T1)*(T2-Tr + &
              (Tr-Tr0)*(1.d0-tanh((h1+p1*(T2-T1)+hp-z0)/hs))/2.d0))

        h1 = h10 - (8.d0*k11/27.d0) + (2.d0*k21) - &
              (3544.d0*k31/2565.d0) + (1859.d0*k41/4104.d0) - &
              (11.d0*k51/40.d0)
        T1 = T10 - (8.d0*k12/27.d0) + (2.d0*k22) - &
              (3544.d0*k32/2565.d0) + (1859.d0*k42/4104.d0) - &
              (11.d0*k52/40.d0)
        T2 = T20 - (8.d0*k13/27.d0) + (2.d0*k23) - &
              (3544.d0*k33/2565.d0) + (1859.d0*k43/4104.d0) - &
              (11.d0*k53/40.d0)

        k61 = h*(-r*(h1+(p1*(T2-T1))/2.d0))
        k62 = h*(-alpha*(T1-Tr)-eps*mu*(T2-T1)**2)
        k63 = h*(-alpha*(T2-Tr)+gt*mu*(T2-T1)*(T2-Tr + &
              (Tr-Tr0)*(1.d0-tanh((h1+p1*(T2-T1)+hp-z0)/hs))/2.d0))

        h10 = h10 + (16.d0*k11/135.d0) + (6656.d0*k31/12825.d0) + &
               (28561.d0*k41/56430.d0) - (9.d0*k51/50.d0) + &
               (2.d0*k61/55.d0)

        T10 = T10 + (16.d0*k12/135.d0) + (6656.d0*k32/12825.d0) + &
               (28561.d0*k42/56430.d0) - (9.d0*k52/50.d0) + &
               (2.d0*k62/55.d0)

        T20 = T20 + (16.d0*k13/135.d0) + (6656.d0*k33/12825.d0) + &
               (28561.d0*k43/56430.d0) - (9.d0*k53/50.d0) + &
               (2.d0*k63/55.d0)

        t0 = t0 + h

        ! --- maxima detection ---
        y3 = y2
        y2 = y1
        y1 = T20

        if (i > 100000000) then
            if ((y2 > y1) .and. (y2 >= y3)) then
                write(1,*) eps, y2
            end if
        end if

    end do

end do

end program elnino
