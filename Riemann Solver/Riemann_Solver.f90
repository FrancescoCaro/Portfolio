program Riemann_Solver

! CONDIZIONI INIZIALI

    REAL a1,a4,u1,u4,p1,p4,gamm,delt,z,u0,err ! Dichiarazione delle variabili   
    
    ! Condizioni iniziali a sinistra
    a1 = 343  ! Velocità del suono nella zona di sinistra in [m/s]
    u1 = 1  ! Velocità del fluido nella zona di sinistra [m/s]
    p1 = 1  ! Pressione del fluido nella zona di sinistra 
    
    ! Condizioni iniziali a destra
    a4 = 343  ! Velocità del suono nella zona di destra in [m/s]
    u4 = 1  ! Velocità del fluido nella zona di destra [m/s]
    p4 = 1  ! Pressione del fluido nella zona di destra
    
    ! Parametri del problema
    gamm = 1.4
    delt = (1.4- 1)*0.5
    
! ALGORITMO RISOLUTIVO DEL PROBLEMA
    
    ! 1 - Calcolo della soluzione di first guess
    z = (a4/a1)*((p1/p4)** ((gamm-1)/2*gamm))
    u0 = (z*(a1+delt*u1)-(a4-delt*u4))/(delt*(1+z))
    
    ! 2 - Onda di sinistra (u-a)
    
    
    ! 3 - Onda di destra (u+a)
    
    
    ! 4 - Iterazione successiva
    
    
    
    
end program Riemann_Solver
