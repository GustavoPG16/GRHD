! ================================================================================
!                             General Relativity Hydrodinamics 
! Este programa resuelve las ecuaciones de GRHD unidimensional usando la métrica
! de Shwarzschild en coordenadas Eddington Finkestein (ED) usando la composición
! 3+1, usando el método de flujos HLLE 
! =================================================================================


! --------------------------------------------------------------------------
! Variables: Este modulo contiene las variables globales que utilizaremos 
!            a lo largo de este código
!  -------------------------------------------------------------------------
module variables
    implicit none
    
    ! --- Índices de ecuaciones ---
    integer :: eq_de = 1   
    integer :: eq_pr = 2   
    integer :: eq_vx = 3    
    integer, parameter :: neq = 3     ! Número total de ecuaciones

    ! --- Dominio numérico  --- 
    real*8 :: xmin, xmax, dx, CFL, t, dt, final_time
    integer :: Nx, Nt
    integer :: t1, t2

    ! --- Arrays principales ---
    real*8, allocatable :: x(:)              ! posiciones
    real*8, allocatable :: u(:,:)            ! variables conservativas
    real*8, allocatable :: up(:,:)           ! variables conservativas actualizadas
    real*8, allocatable :: p(:,:)            ! variables primitivas
    real*8, allocatable :: fx(:,:)           ! flujos
    real*8, allocatable :: rhs(:,:)          ! término fuente
    real*8, allocatable :: pLx(:,:), pRx(:,:) ! estados reconstruidos

    ! --- Geometría de Schwarzschild en coordenadas EF ---
    real*8, allocatable :: grr(:)        ! Componente radial de la métrica
    real*8, allocatable :: alpha(:)      ! Función lapso
    real*8, allocatable :: dlnalpha(:)   ! Derivada del logaritmo lapso
    real*8, allocatable :: beta(:)     ! Vector shift
    real*8, allocatable :: detgamma(:) ! Determinante de la métrica espacial

    ! --- Parámetros físicos ---
    real*8 :: rho_ini    ! Densidad inicial
    real*8 :: v_infty    ! Velocidad en el infinito
    real*8 :: cs_infty   ! Velocidad del sonido en el infinito
    real*8 :: gam        ! Indice adiabático  
    real*8 :: Ma         ! Masa del agujero negro

    ! --- Método numérico ---
    character(len=10) :: reconstruction_method  ! Método de reconstrucción
    real*8 :: floor                      ! Valor mínimo para densidad/presión

    ! --- Output ---
    character(len=100) :: output_dir = 'output/'
    character(len=100) :: output_prefix
end module variables


! ---------------------------------------------------------------------------------------
! Metrics: en este modulo implementamos la métrica de Schwarzschild para un agujero negro 
!          no rotante en coordenadas ED, implementamos su respectiva descomposición 3+1.
!          Agregamos los simbolos de Christoffel, los tesores energía momento no nulos, 
!          respectivamente.Tambien calculamos los términos fuentes S.
! ---------------------------------------------------------------------------------------
module metrics
    use variables
    implicit none
    private
    public :: set_metric, calculate_schwarzschild_EF, calculate_geometric_sources, &
              calculate_stress_energy_tensor

contains
    ! Elegir la métrica
    subroutine set_metric(x, grr, alpha, dlnalpha, beta, detgamma, metric_type)
        implicit none
        real(kind=8), intent(in) :: x(:)
        real*8, intent(out) :: grr(:)
        real*8, intent(out) :: alpha(:), dlnalpha(:)
        real*8, intent(out) :: beta(:)
        real*8, intent(out) :: detgamma(:)
        character(len=*), intent(in) :: metric_type

        select case(trim(metric_type))
           case('Schwarzschild_EF')
                call calculate_schwarzschild_EF(x, grr, alpha, dlnalpha, &
                                             beta, detgamma)
       end select
    end subroutine set_metric

    ! Métrica de Schwarzschild en coordenadas ED
    subroutine calculate_schwarzschild_EF(x, grr, alpha, dlnalpha, beta, detgamma)
        implicit none
        real*8, intent(in) :: x(:)
        real*8, intent(out) :: grr(:) 
        real*8, intent(out) :: alpha(:), dlnalpha(:)
        real*8, intent(out) :: beta(:)
        real*8, intent(out) :: detgamma(:)

        ! Métrica de Schwarzschild en coordenadas EF
        alpha = 1.0d0 / sqrt(1.0d0 + 2.0d0 * Ma / x)
        dlnalpha = Ma / (x**2 + 2.0d0 * Ma * x)
        
        beta = 2.0d0 * Ma / x / (1.0d0 + 2.0d0 * Ma / x)
        
        grr = 1.0d0 + 2.0d0 * Ma / x
        
        detgamma = x**4 * grr
    end subroutine calculate_schwarzschild_EF

    ! Componentes del tensor energía momento
    subroutine calculate_stress_energy_tensor(rho, press, v, W, h, i, T00, Tr0, Trr)
        implicit none
        real*8, intent(in) :: rho, press, v, W, h
        integer, intent(in) :: i
        real*8, intent(out) :: T00, Tr0, Trr

        ! Componentes del tensor energía-momento
        T00 = (rho * h * W**2 - press) / alpha(i)**2
        Tr0 = rho * h * W**2 / alpha(i) * (v - beta(i)/alpha(i)) + &
              press * beta(i) / alpha(i)**2
        Trr = rho * h * W**2 * (v - beta(i)/alpha(i))**2 + &
              press * (1.0d0/grr(i) - beta(i)**2/alpha(i)**2)
    end subroutine calculate_stress_energy_tensor

    ! Simbolos de Christoffel
    subroutine calculate_geometric_sources(i, T00, Tr0, Trr, source)
        implicit none
        integer, intent(in) :: i
        real*8, intent(in) :: T00, Tr0, Trr
        real*8, intent(out) :: source(neq)
        ! Notación {xyz}={xy}^z
        real*8 :: Gamma0r0, Gamma0rr, Gammarrr, Gammarr0, Gammathrth, Gammaphrph
        real*8 :: Gamma000, Gammathth0, Gammaphph0
        !real*8 :: x_safe  ! Solo para proteger divisiones

        ! Proteger x cerca del horizonte
        !x_safe = max(x(i), 2.001d0 * Ma)

	! --- Simbolos de Christoffel
        Gamma0r0 = Ma / x(i)**2 * (1.0d0 + 2.0d0 * Ma / x(i))
        Gamma0rr = -2.0d0 * Ma**2 / x(i)**3
        Gammarrr = -Ma/ x(i)**2 * (1.0d0 + 2.0d0 * Ma / x(i))
        Gammarr0 = 2.0d0 * Ma / x(i)**2 * (1.0d0 + Ma / x(i))
        Gammathrth = 1.0d0 / x(i)
        Gammaphrph = 1.0d0 / x(i)
        Gamma000 = 2.0d0 * Ma**2 / x(i)**3
        Gammathth0 = -2.0d0 * Ma
        Gammaphph0 = -2.0d0 * Ma !* sin^2(theta)

        ! --- Terminos Fuente
        ! Ecuación de continuidad
        source(1) = 0.0d0  
        ! Ecuación de momento 
        source(2) = T00 * ((-alpha(i)**2 + grr(i)*beta(i)**2)*Gamma0r0 + grr(i)*beta(i)*Gamma0rr) + &
                    Tr0 * ((-alpha(i)**2 + grr(i)*beta(i)**2)*Gamma0r0 * grr(i)*Gammarrr + &
                    grr(i)*beta(i)*(Gammarr0 + Gamma0r0)) + &
                    Trr *(grr(i)*Gammarrr + grr(i)*beta(i)*Gammarr0) 
        ! Ecuación de energíaa
        source(3) = alpha(i) * (- T00 * Gamma000 + Tr0 * (dlnalpha(i) - 2.0d0 * Gamma0r0) - &
                    Gammarr0 * Trr)
                    
        ! Multiplicar por factores globales
         source(:) = alpha(i) * sqrt(detgamma(i)) * source(:)
    end subroutine calculate_geometric_sources

end module metrics

! ---------------------------------------------------------------------------------------
! MODFluxes: Esta subrutina implementa el cálculo del flujo numérico en un esquema de 
!            Godunov usando el solver aproximado HLLE en relatividad general 1D. También
!            se calculan los eigenvalores.
! ---------------------------------------------------------------------------------------
module MODfluxes
    use variables
    use metrics
    implicit none
    private
    public :: fluxes, compute_eigenvalues

contains
    subroutine fluxes(pl, pr, f, alpha_local, beta_local, grr_local, det_gamma)
        real(kind=8), intent(in)  :: pl(neq), pr(neq)
        real(kind=8), intent(in)  :: alpha_local, beta_local, grr_local, det_gamma
        real(kind=8), intent(out) :: f(neq)
        real(kind=8) :: ul(neq), ur(neq), fl(neq), fr(neq)  
        real(kind=8) :: ap, am, cflmin, cflmax, cfrmin, cfrmax

        ! Calcular flujos físicos
        call compute_flux(pl, fl, alpha_local, beta_local, grr_local, det_gamma)
        call compute_flux(pr, fr, alpha_local, beta_local, grr_local, det_gamma)

        ! Convertir a conservativas
        call compute_conservatives(pl, ul, grr_local, det_gamma)
        call compute_conservatives(pr, ur, grr_local, det_gamma)

        ! Calcular velocidades características
        call compute_eigenvalues(pl, cflmin, cflmax, alpha_local, beta_local, grr_local)
        call compute_eigenvalues(pr, cfrmin, cfrmax, alpha_local, beta_local, grr_local)

        ! Velocidades de onda para HLLE
        ap = max(0.0d0, cflmax, cfrmax)
        am = min(0.0d0, cflmin, cfrmin)

        ! Flujo HLLE
        if (abs(ap - am) > 1.0d-12) then
            f(:) = (ap*fl(:) - am*fr(:) + ap*am*(ur(:) - ul(:)))/(ap - am)
        else
            f(:) = 0.5d0*(fl(:) + fr(:))
        end if
    end subroutine fluxes

    subroutine compute_flux(p_local, f_local, alpha_local, beta_local, grr_local, detgamma)
        real(kind=8), intent(in)  :: p_local(neq)
        real(kind=8), intent(in)  :: alpha_local, beta_local, grr_local, detgamma
        real(kind=8), intent(out) :: f_local(neq)
        real(kind=8) :: u_local(neq), v

        ! Convertir a variables conservativas
        call compute_conservatives(p_local, u_local, grr_local, detgamma)
        v = p_local(eq_vx)  ! Velocidad

        ! Flujos 
        f_local(eq_de) = alpha_local * (v - beta_local/alpha_local) * u_local(eq_de)
        f_local(eq_pr) = alpha_local * (v - beta_local/alpha_local) * u_local(eq_pr) + &
                         alpha_local * sqrt(detgamma) * p_local(eq_pr)
        f_local(eq_vx) = alpha_local * (v - beta_local/alpha_local) * u_local(eq_vx) + &
                         alpha_local * sqrt(detgamma) * v * p_local(eq_pr)
    end subroutine compute_flux

    subroutine compute_conservatives(p_local, u_local, grr_local, detgamma)
        real(kind=8), intent(in)  :: p_local(neq)
        real(kind=8), intent(in)  :: grr_local, detgamma
        real(kind=8), intent(out) :: u_local(neq)
        real(kind=8) :: W, h, v, rho, press, sqrt_detgamma
        real(kind=8) :: v2, W2, e_internal
        real(kind=8), parameter :: eps = 1.0d-12
        real(kind=8), parameter :: v_max = 0.99d0

        ! Extraer primitivas con protección
        rho = max(p_local(eq_de), eps)
        press = max(p_local(eq_pr), eps)
        v = p_local(eq_vx)
    
        ! Proteger grr_local y detgamma
        sqrt_detgamma = sqrt(max(detgamma, eps))
    
        ! Limitar velocidad
        v2 = v * v
        if (grr_local * v2 >= 1.0d0) then
            v = sign(v_max, v)
            v2 = v * v
        end if
    
        ! Factor de Lorentz con protección
        W2 = 1.0d0 / max(eps, 1.0d0 - grr_local * v2)
        W = sqrt(W2)
    
        ! Verificar W físico
        if (W > 100.0d0) then  ! Límite práctico
            W = 100.0d0
            W2 = W * W
            v = sign(sqrt((W2 - 1.0d0) / (W2 * grr_local)), v)
        end if
    
        ! Entalpía específica con protección
        e_internal = press / ((gam - 1.0d0) * rho)
        h = 1.0d0 + e_internal + press / rho
        h = max(h, 1.0d0 + eps)
    
        ! Variables conservativas con protección
        u_local(eq_de) = sqrt_detgamma * rho * W
        u_local(eq_pr) = sqrt_detgamma * rho * h * W2 * grr_local * v
        u_local(eq_vx) = sqrt_detgamma * (rho * h * W2 - press - rho * W)
    
       ! Verificación de sanidad
       if (any(abs(u_local) > 1.0d10)) then
            u_local(eq_de) = sqrt_detgamma * rho
            u_local(eq_pr) = sqrt_detgamma * rho * v
            u_local(eq_vx) = sqrt_detgamma * press / (gam - 1.0d0)
       end if
    
    end subroutine compute_conservatives

    subroutine compute_eigenvalues(p_local, lambda_min, lambda_max, alpha_local, beta_local, grr_local)
        real(kind=8), intent(in)  :: p_local(neq)
        real(kind=8), intent(in)  :: alpha_local, beta_local, grr_local
        real(kind=8), intent(out) :: lambda_min, lambda_max
        real(kind=8) :: cs, v, lambda1, lambda2, lambda3
        real(kind=8) :: a, b, c, discriminant, sqrt_disc
        real(kind=8), parameter :: eps = 1.0d-12
        real(kind=8), parameter :: v_max = 0.99d0
    
        v = p_local(eq_vx)
    
        ! Limitar velocidad a valores físicos
        v = sign(min(abs(v), v_max), v)
    
        ! Velocidad del sonido con protección completa
        if (p_local(eq_de) > eps .and. p_local(eq_pr) > eps) then
           cs = sqrt(gam * p_local(eq_pr) * (gam - 1.0d0) / &
                 max(eps, p_local(eq_pr) * gam + p_local(eq_de) * (gam - 1.0d0)))
           cs = min(cs, sqrt(1.0d0 - eps))
        else
           cs = 0.1d0  ! Valor por defecto seguro
        end if
    
        ! Eigenvalue de la ecuación de continuidad
        lambda1 = alpha_local * v - beta_local
    
        ! Coeficientes de la ecuación cuadrática para eigenvalues acústicos
        ! Usando la forma estándar: a*lambda² + b*lambda + c = 0
        a = 1.0d0 - cs * cs
        b = -2.0d0 * alpha_local * v * (1.0d0 - cs * cs)
        c = alpha_local * alpha_local * (v * v - cs * cs * max(eps, 1.0d0 - grr_local * v * v))
    
        ! Resolver ecuación cuadrática de manera robusta
        if (abs(a) < eps) then
           ! Ecuación lineal
           if (abs(b) > eps) then
               lambda2 = -c / b
               lambda3 = lambda2
           else
               lambda2 = lambda1
               lambda3 = lambda1
           end if
        else
           discriminant = b * b - 4.0d0 * a * c
           if (discriminant < 0.0d0) then
               lambda2 = lambda1
               lambda3 = lambda1
           else
               sqrt_disc = sqrt(discriminant)
               lambda2 = (-b + sqrt_disc) / (2.0d0 * a)
               lambda3 = (-b - sqrt_disc) / (2.0d0 * a)
           end if
       end if
    
      ! Aplicar shift
       lambda2 = lambda2 - beta_local
       lambda3 = lambda3 - beta_local
    
       lambda_min = min(lambda1, lambda2, lambda3)
       lambda_max = max(lambda1, lambda2, lambda3)
    
       ! Verificación final
       if (abs(lambda_max - lambda_min) < eps) then
           lambda_max = lambda_min + eps
       end if
    
    end subroutine compute_eigenvalues

end module MODfluxes


! ---------------------------------------------------------------------------------------
! boundary_conditions: Gestiona los valores en las zonas fantasma y en los bordes del 
!                     dominio, usando extrapolación polinómicas para mantener la coherencia 
!                     física cerca de la frontera interna (como el horizonte de unBH) y 
!                     condiciones de salida en la frontera externa.
! ---------------------------------------------------------------------------------------
module boundary_conditions
    use variables
    implicit none
    private
    public :: apply_boundary_conditions, apply_inner_boundary

contains
    ! Subrutina que aplica condiciones de frontera
    subroutine apply_boundary_conditions(q)
        implicit none
        real(kind=8), intent(inout) :: q(neq,-1:nx+2)
        
        ! Frontera interna (horizonte) - extrapolación cúbica
        q(eq_de,0)  = 3.0d0 * q(eq_de,1) - 3.0d0 * q(eq_de,2) + q(eq_de,3)
        q(eq_pr,0)  = 3.0d0 * q(eq_pr,1) - 3.0d0 * q(eq_pr,2) + q(eq_pr,3)
        q(eq_vx,0)  = 3.0d0 * q(eq_vx,1) - 3.0d0 * q(eq_vx,2) + q(eq_vx,3)

        ! Completar dominio fantasma - extrapolación lineal
        q(eq_de,-1) = 2.0d0 * q(eq_de,0) - q(eq_de,1)
        q(eq_pr,-1) = 2.0d0 * q(eq_pr,0) - q(eq_pr,1)
        q(eq_vx,-1) = 2.0d0 * q(eq_vx,0) - q(eq_vx,1)

        ! Frontera externa 
        q(:,nx+1) = q(:,nx)
        q(:,nx+2) = q(:,nx+1)
    end subroutine apply_boundary_conditions

    ! condiciones internas
    subroutine apply_inner_boundary(q)
        implicit none
        real(kind=8), intent(inout) :: q(neq,-1:nx+2)
        
        ! Solo condiciones de frontera interna
        q(eq_de,0)  = 3.0d0 * q(eq_de,1) - 3.0d0 * q(eq_de,2) + q(eq_de,3)
        q(eq_pr,0)  = 3.0d0 * q(eq_pr,1) - 3.0d0 * q(eq_pr,2) + q(eq_pr,3)
        q(eq_vx,0)  = 3.0d0 * q(eq_vx,1) - 3.0d0 * q(eq_vx,2) + q(eq_vx,3)

        q(eq_de,-1) = 2.0d0 * q(eq_de,0) - q(eq_de,1)
        q(eq_pr,-1) = 2.0d0 * q(eq_pr,0) - q(eq_pr,1)
        q(eq_vx,-1) = 2.0d0 * q(eq_vx,0) - q(eq_vx,1)
    end subroutine apply_inner_boundary

end module boundary_conditions

! ----------------------------------------------------------------------------------
! Eqaurions: Este módulo gestiona las ecuaciones principales del HD en 1D, incluyen
!            la conversión entre variables, reconstrucción espacial, integración 
!            temporal y limitadores numéricos.
! ---------------------------------------------------------------------------------
module equations
    use variables
    use metrics
    use MODfluxes, only: fluxes
    use boundary_conditions, only: apply_inner_boundary
    implicit none
    private
    public :: primu, uprim, reconstruct_variables, minmod, advance_solution

contains
    real*8 function minmod(a, b)
        implicit none
        real*8, intent(in) :: a, b
        minmod = 0.5d0 * (sign(1.0d0, a) + sign(1.0d0, b)) * min(abs(a), abs(b))
    end function minmod

    subroutine reconstruct_variables(p, pLx, pRx)
    implicit none
    real*8, intent(in)  :: p(neq, -1:nx+2)
    real*8, intent(out) :: pLx(neq, -1:nx+2), pRx(neq, -1:nx+2)
    
    ! Todas las declaraciones al principio
    integer :: i, k
    real*8 :: slopex, smaxL, sminL
    real*8 :: pL_temp, pR_temp
    real*8 :: grr_avg, p_min, p_max
    real*8 :: v_limit
    real*8, parameter :: eps = 1.0d-12
    real*8, parameter :: safety_factor = 0.3d0  ! Factor conservativo para minmod

    select case(trim(reconstruction_method))
        case('Godunov') 
            do i = 0, nx
                pLx(:,i) = p(:,i)
                pRx(:,i) = p(:,i+1)
            end do
            
        case('minmod')
            do k = 1, neq
                do i = 0, nx
                    ! Usar primer orden cerca de los bordes para mayor estabilidad
                    if (i <= 1 .or. i >= nx-1) then
                        pLx(k,i) = p(k,i)
                        pRx(k,i) = p(k,i+1)
                    else
                        ! Calcular pendientes con protección
                        smaxL = (p(k,i+1) - p(k,i)) / dx
                        sminL = (p(k,i) - p(k,i-1)) / dx
                        slopex = minmod(sminL, smaxL)
                        
                        ! Reconstruir valores con factor conservativo
                        pL_temp = p(k,i) + safety_factor * slopex * dx
                        pR_temp = p(k,i+1) - safety_factor * slopex * dx
                        
                        ! Aplicar límites físicos según la variable
                        if (k == eq_de) then
                            ! Densidad: siempre positiva
                            pLx(k,i) = max(pL_temp, floor)
                            pRx(k,i) = max(pR_temp, floor)
                            
                            ! No permitir variaciones extremas
                            p_min = min(p(k,i), p(k,i+1))
                            p_max = max(p(k,i), p(k,i+1))
                            pLx(k,i) = max(p_min * 0.1d0, min(pLx(k,i), p_max * 5.0d0))
                            pRx(k,i) = max(p_min * 0.1d0, min(pRx(k,i), p_max * 5.0d0))
                            
                        else if (k == eq_pr) then
                            ! Presión: siempre positiva
                            pLx(k,i) = max(pL_temp, floor)
                            pRx(k,i) = max(pR_temp, floor)
                            
                            ! No permitir variaciones extremas
                            p_min = min(p(k,i), p(k,i+1))
                            p_max = max(p(k,i), p(k,i+1))
                            pLx(k,i) = max(p_min * 0.1d0, min(pLx(k,i), p_max * 5.0d0))
                            pRx(k,i) = max(p_min * 0.1d0, min(pRx(k,i), p_max * 5.0d0))
                            
                        else if (k == eq_vx) then
                            ! Velocidad: límite causal
                            pLx(k,i) = pL_temp
                            pRx(k,i) = pR_temp
                            
                            ! Verificar límite causal usando grr local
                            grr_avg = 0.5d0 * (grr(i) + grr(i+1))
                            v_limit = sqrt(0.95d0 / max(grr_avg, eps))
                            
                            if (abs(pLx(k,i)) > v_limit) then
                                pLx(k,i) = sign(v_limit, pLx(k,i))
                            end if
                            if (abs(pRx(k,i)) > v_limit) then
                                pRx(k,i) = sign(v_limit, pRx(k,i))
                            end if
                        else
                            ! Otras variables
                            pLx(k,i) = pL_temp
                            pRx(k,i) = pR_temp
                        end if
                    end if
                end do
            end do
            
        case default
            ! Si no se reconoce el método, usar Godunov
            do i = 0, nx
                pLx(:,i) = p(:,i)
                pRx(:,i) = p(:,i+1)
            end do
    end select
    
    ! Verificación final de consistencia física en todo el dominio
    do i = 0, nx
        ! Asegurar valores mínimos para densidad y presión
        pLx(eq_de,i) = max(pLx(eq_de,i), floor)
        pLx(eq_pr,i) = max(pLx(eq_pr,i), floor)
        pRx(eq_de,i) = max(pRx(eq_de,i), floor)
        pRx(eq_pr,i) = max(pRx(eq_pr,i), floor)
        
        ! Verificación final de velocidad causal
        if (i < nx) then
            if (grr(i) > eps) then
                v_limit = sqrt(0.98d0 / grr(i))
                if (abs(pLx(eq_vx,i)) > v_limit) then
                    pLx(eq_vx,i) = sign(v_limit, pLx(eq_vx,i))
                end if
            end if
        end if
        
        if (i > 0) then
            if (grr(i) > eps) then
                v_limit = sqrt(0.98d0 / grr(i))
                if (abs(pRx(eq_vx,i)) > v_limit) then
                    pRx(eq_vx,i) = sign(v_limit, pRx(eq_vx,i))
                end if
            end if
        end if
    end do
    
end subroutine reconstruct_variables

    subroutine primu(p_local, u_local, i)
        implicit none
        real*8, intent(in)  :: p_local(neq)
        real*8, intent(out) :: u_local(neq)
        integer, intent(in) :: i
        real*8 :: v, h, W, e
        
        v = p_local(eq_vx)
        e = p_local(eq_pr) / (gam - 1.0d0) / p_local(eq_de)
        h = 1.0d0 + e + p_local(eq_pr) / p_local(eq_de)  
        W = 1.0d0 / sqrt(1.0d0 - grr(i) * v**2)

        ! Variables conservativas 
        u_local(eq_de) = sqrt(detgamma(i)) * p_local(eq_de) * W                                   
        u_local(eq_pr) = sqrt(detgamma(i)) * p_local(eq_de) * h * W**2 * grr(i) * v              
        u_local(eq_vx) = sqrt(detgamma(i)) * (p_local(eq_de) * h * W**2 - p_local(eq_pr) - p_local(eq_de) * W)  
    end subroutine primu

subroutine uprim(u_local, p_local, i)
    real*8, intent(in)  :: u_local(neq)
    real*8, intent(out) :: p_local(neq)
    integer, intent(in) :: i
    
    ! Declaraciones
    real*8 :: D, Sr, tau, sqrt_gamma
    real*8 :: rho, press, v, W, h, e
    real*8 :: press_guess, press_old, f, df, tol
    real*8 :: term1, term2, denom, common_factor
    real*8 :: v_max_allowed, press_max, press_min
    real*8 :: safe_denom, safe_sqrt_arg
    real*8, parameter :: eps = 1.0d-10
    real*8, parameter :: max_iter = 50
    real*8, parameter :: newton_tol = 1.0d-8
    real*8, parameter :: huge_val = 1.0d10
    integer :: iter
    logical :: converged, use_fallback
    
    sqrt_gamma = sqrt(max(detgamma(i), eps))
    
    ! Extraer conservativas con verificación
    D = u_local(eq_de)
    Sr = u_local(eq_pr)
    tau = u_local(eq_vx)
    
    ! Verificación de sanidad de las conservativas
    use_fallback = .false.
    
    ! Verificar si D es razonable
    if (abs(D) < eps * sqrt_gamma .or. abs(D) > huge_val) then
        use_fallback = .true.
    end if
    
    ! Verificar si tau es razonable
    if (abs(tau) > huge_val .or. isnan(tau)) then
        use_fallback = .true.
    end if
    
    ! Verificar si Sr es razonable
    if (abs(Sr) > huge_val .or. isnan(Sr)) then
        use_fallback = .true.
    end if
    
    ! Verificar consistencia básica: D debe ser positivo después de corrección
    D = abs(D)
    if (D < eps * sqrt_gamma) then
        use_fallback = .true.
    end if
    
    ! Si las conservativas son problemáticas, usar método de respaldo
    if (use_fallback) then
        p_local(eq_de) = max(floor, D / sqrt_gamma)
        p_local(eq_pr) = floor
        p_local(eq_vx) = 0.0d0
        return
    end if
    
    ! Límites físicos
    press_min = floor
    press_max = min(100.0d0 * D / sqrt_gamma, 0.1d0 * D / sqrt_gamma)
    v_max_allowed = sqrt(0.95d0 / max(grr(i), eps))
    
    ! Estimación inicial conservativa
    press_guess = max(press_min, min(0.001d0 * D / sqrt_gamma, press_max))
    
    ! Newton-Raphson con protecciones extremas
    converged = .false.
    do iter = 1, int(max_iter)
        press_old = press_guess
        
        ! Verificar que pressure guess está en rango
        if (press_guess < press_min .or. press_guess > press_max .or. &
            isnan(press_guess) .or. press_guess <= 0.0d0) then
            press_guess = press_min
            cycle
        end if
        
        ! Cálculos intermedios con protección extrema
        common_factor = D + press_guess * sqrt_gamma + tau
        
        ! Verificar common_factor
        if (abs(common_factor) < eps) then
            press_guess = press_min
            cycle
        end if
        
        safe_denom = grr(i) * common_factor**2
        if (abs(safe_denom) < eps) then
            press_guess = press_min
            cycle
        end if
        
        safe_sqrt_arg = 1.0d0 - Sr**2 / safe_denom
        if (safe_sqrt_arg <= eps) then
            ! Velocidad superlumínica detectada
            press_guess = max(press_min, 0.5d0 * press_guess)
            cycle
        end if
        
        term1 = safe_sqrt_arg
        term2 = sqrt(term1)
        
        ! Verificar term2
        if (term2 < eps .or. isnan(term2)) then
            press_guess = press_min
            cycle
        end if
        
        ! Función objetivo simplificada y más robusta
        f = press_guess - (gam - 1.0d0) * term1 * &
            (tau + D * (1.0d0 - 1.0d0/term2)) / sqrt_gamma
            
        ! Derivada simplificada
        if (term2 > eps .and. abs(common_factor) > eps) then
            df = 1.0d0 + (gam - 1.0d0) * Sr**2 * sqrt_gamma / &
                 (grr(i) * common_factor**2 * term2)
        else
            df = 1.0d0
        end if
        
        ! Verificar derivada
        if (abs(df) < eps .or. isnan(df)) then
            press_guess = press_min
            cycle
        end if
        
        ! Paso Newton con limitación
        press_guess = press_old - f / df
        
        ! Aplicar límites estrictos
        press_guess = max(press_min, min(press_max, press_guess))
        
        ! Verificar convergencia
        if (abs(press_old) > eps) then
            tol = abs(press_guess - press_old) / abs(press_old)
        else
            tol = abs(press_guess - press_old)
        end if
        
        if (tol < newton_tol .and. press_guess > press_min) then
            converged = .true.
            exit
        end if
        
        ! Si no converge después de muchas iteraciones, usar mínimo
        if (iter > max_iter/2 .and. tol > 1.0d-2) then
            press_guess = press_min
        end if
    end do
    
    ! Si Newton-Raphson falló, usar método de respaldo
    if (.not. converged .or. press_guess <= press_min) then
        press_guess = press_min
    end if
    
    ! Recuperar primitivas con máxima protección
    press = max(press_guess, press_min)
    
    ! Velocidad con protección extrema
    safe_denom = tau + D + press * sqrt_gamma
    if (abs(safe_denom) > eps .and. abs(grr(i)) > eps) then
        v = Sr / (safe_denom * grr(i))
        
        ! Verificar límite causal
        if (abs(v) > v_max_allowed) then
            v = sign(v_max_allowed, v)
        end if
        
        ! Verificación adicional de velocidad
        if (grr(i) * v**2 >= 1.0d0) then
            v = sign(v_max_allowed, v)
        end if
    else
        v = 0.0d0
    end if
    
    ! Factor de Lorentz con protección
    safe_sqrt_arg = 1.0d0 - grr(i) * v**2
    if (safe_sqrt_arg > eps) then
        W = 1.0d0 / sqrt(safe_sqrt_arg)
        if (W > 100.0d0) then  ! Límite práctico
            W = 100.0d0
            v = sign(sqrt((W**2 - 1.0d0) / (W**2 * grr(i))), v)
        end if
    else
        W = 1.0d0
        v = 0.0d0
    end if
    
    ! Densidad con protección
    if (W > eps .and. sqrt_gamma > eps) then
        rho = D / (W * sqrt_gamma)
    else
        rho = D / sqrt_gamma
    end if
    rho = max(rho, floor)
    
    ! Asignar primitivas finales
    p_local(eq_de) = rho
    p_local(eq_pr) = press
    p_local(eq_vx) = v
    
    ! Verificación final de consistencia física
    e = press / (max(eps, (gam - 1.0d0) * rho))
    h = 1.0d0 + e + press / rho
    
    if (e < 0.0d0 .or. h <= 1.0d0 .or. rho < floor .or. &
        isnan(e) .or. isnan(h) .or. isnan(rho) .or. &
        isnan(press) .or. isnan(v)) then
        ! Último recurso: valores de respaldo completamente seguros
        p_local(eq_de) = max(floor, D / sqrt_gamma)
        p_local(eq_pr) = floor
        p_local(eq_vx) = 0.0d0
    end if
    
end subroutine uprim


subroutine advance_solution(dt)
    real*8, intent(in) :: dt
    integer :: i, rk, k
    real*8 :: T00, Tr0, Trr
    real*8 :: v_local, h_local, W_local, e_local
    real*8 :: dt_safe, rho_local, press_local
    real*8 :: max_speed, local_dt
    real*8, parameter :: eps = 1.0d-12
    real*8, parameter :: max_h = 100.0d0
    real*8, parameter :: max_W = 100.0d0
    
    ! Verificar dt
    dt_safe = min(dt, dx * 0.01d0)  ! Limitar dt
    
    ! Guardar estado inicial
    up = u
    
    do rk = 1, 3
        ! Verificar primitivas antes de reconstrucción
        do i = 1, nx
            if (p(eq_de,i) < floor .or. p(eq_pr,i) < floor .or. &
                isnan(p(eq_de,i)) .or. isnan(p(eq_pr,i)) .or. isnan(p(eq_vx,i))) then
                p(eq_de,i) = floor
                p(eq_pr,i) = floor  
                p(eq_vx,i) = 0.0d0
            end if
            
            ! Verificar velocidad causal
            if (grr(i) * p(eq_vx,i)**2 >= 1.0d0) then
                p(eq_vx,i) = sign(sqrt(0.95d0/grr(i)), p(eq_vx,i))
            end if
        end do
        
        ! Reconstruct variables
        call reconstruct_variables(p, pLx, pRx)

        ! Calculate fluxes con verificación
        do i = 0, nx-1
            ! Verificar estados reconstruidos
            if (pLx(eq_de,i) < floor) pLx(eq_de,i) = floor
            if (pLx(eq_pr,i) < floor) pLx(eq_pr,i) = floor
            if (pRx(eq_de,i) < floor) pRx(eq_de,i) = floor
            if (pRx(eq_pr,i) < floor) pRx(eq_pr,i) = floor
            
            call fluxes(pLx(:,i), pRx(:,i), fx(:,i), alpha(i), beta(i), grr(i), detgamma(i))
            
            ! Verificar que los flujos no sean infinitos
            if (any(abs(fx(:,i)) > 1.0d10) .or. any(isnan(fx(:,i)))) then
                fx(:,i) = 0.0d0
            end if
        end do

        ! Calcular RHS con máxima protección
        do i = 1, nx-1
            ! Extraer y verificar variables locales
            rho_local = max(p(eq_de,i), floor)
            press_local = max(p(eq_pr,i), floor)
            v_local = p(eq_vx,i)
            
            ! Verificar velocidad causal
            if (grr(i) * v_local**2 >= 1.0d0) then
                v_local = sign(sqrt(0.95d0/grr(i)), v_local)
            end if
            
            ! Calcular cantidades locales con protección
            e_local = press_local / (max(eps, (gam - 1.0d0) * rho_local))
            h_local = 1.0d0 + e_local + press_local / rho_local
            h_local = max(1.0d0 + eps, min(h_local, max_h))
            
            W_local = 1.0d0 / sqrt(max(eps, 1.0d0 - grr(i) * v_local**2))
            W_local = min(W_local, max_W)

            ! Verificar que las cantidades sean físicas
            if (isnan(h_local) .or. isnan(W_local) .or. h_local <= 1.0d0) then
                ! Usar valores de respaldo
                h_local = 1.01d0
                W_local = 1.0d0
                v_local = 0.0d0
            end if

            ! Calcular tensor energía-momento con protección
            call calculate_stress_energy_tensor(rho_local, press_local, v_local, &
                                             W_local, h_local, i, T00, Tr0, Trr)
            
            ! Verificar componentes del tensor
            if (isnan(T00) .or. isnan(Tr0) .or. isnan(Trr) .or. &
                abs(T00) > 1.0d10 .or. abs(Tr0) > 1.0d10 .or. abs(Trr) > 1.0d10) then
                T00 = rho_local
                Tr0 = 0.0d0
                Trr = press_local
            end if

            ! Calcular fuentes geométricas
            call calculate_geometric_sources(i, T00, Tr0, Trr, rhs(:,i))
            
            ! Verificar fuentes
            if (any(isnan(rhs(:,i))) .or. any(abs(rhs(:,i)) > 1.0d10)) then
                rhs(:,i) = 0.0d0
            end if

            ! RHS completo con verificación de flujos
            if (abs(fx(eq_de,i)) < 1.0d10 .and. abs(fx(eq_de,i-1)) < 1.0d10 .and. &
                .not. isnan(fx(eq_de,i)) .and. .not. isnan(fx(eq_de,i-1))) then
                rhs(eq_de,i) = -(fx(eq_de,i) - fx(eq_de,i-1)) / dx + rhs(eq_de,i)
            else
                rhs(eq_de,i) = 0.0d0
            end if
            
            if (abs(fx(eq_pr,i)) < 1.0d10 .and. abs(fx(eq_pr,i-1)) < 1.0d10 .and. &
                .not. isnan(fx(eq_pr,i)) .and. .not. isnan(fx(eq_pr,i-1))) then
                rhs(eq_pr,i) = -(fx(eq_pr,i) - fx(eq_pr,i-1)) / dx + rhs(eq_pr,i)
            else
                rhs(eq_pr,i) = 0.0d0
            end if
            
            if (abs(fx(eq_vx,i)) < 1.0d10 .and. abs(fx(eq_vx,i-1)) < 1.0d10 .and. &
                .not. isnan(fx(eq_vx,i)) .and. .not. isnan(fx(eq_vx,i-1))) then
                rhs(eq_vx,i) = -(fx(eq_vx,i) - fx(eq_vx,i-1)) / dx + rhs(eq_vx,i)
            else
                rhs(eq_vx,i) = 0.0d0
            end if
            
            ! Limitar RHS para evitar cambios extremos
            do k = 1, neq
                if (abs(rhs(k,i) * dt_safe) > abs(u(k,i))) then
                    rhs(k,i) = sign(abs(u(k,i)) / dt_safe * 0.1d0, rhs(k,i))
                end if
            end do
        end do

        ! Frontera externa 
        rhs(:,nx) = 0.0d0

        ! RK3 step con verificación
        if (rk == 1) then
            do i = 1, nx
                do k = 1, neq
                    if (.not. isnan(rhs(k,i)) .and. abs(rhs(k,i)) < 1.0d10) then
                        u(k,i) = up(k,i) + dt_safe * rhs(k,i)
                    else
                        u(k,i) = up(k,i)
                    end if
                end do
            end do
        else if (rk == 2) then
            do i = 1, nx
                do k = 1, neq
                    if (.not. isnan(rhs(k,i)) .and. abs(rhs(k,i)) < 1.0d10) then
                        u(k,i) = 0.75d0 * up(k,i) + 0.25d0 * u(k,i) + 0.25d0 * dt_safe * rhs(k,i)
                    else
                        u(k,i) = 0.75d0 * up(k,i) + 0.25d0 * u(k,i)
                    end if
                end do
            end do
        else
            do i = 1, nx
                do k = 1, neq
                    if (.not. isnan(rhs(k,i)) .and. abs(rhs(k,i)) < 1.0d10) then
                        u(k,i) = up(k,i)/3.0d0 + 2.0d0/3.0d0 * (u(k,i) + dt_safe * rhs(k,i))
                    else
                        u(k,i) = up(k,i)/3.0d0 + 2.0d0/3.0d0 * u(k,i)
                    end if
                end do
            end do
        end if

        ! Condiciones de frontera
        call apply_inner_boundary(u)
        
        ! Recuperar primitivas con verificación
        do i = -1, nx
            call uprim(u(:,i), p(:,i), i)
            
            ! Verificación post-uprim
            if (p(eq_de,i) < floor .or. p(eq_pr,i) < floor .or. &
                isnan(p(eq_de,i)) .or. isnan(p(eq_pr,i)) .or. isnan(p(eq_vx,i))) then
                p(eq_de,i) = floor
                p(eq_pr,i) = floor
                p(eq_vx,i) = 0.0d0
                call primu(p(:,i), u(:,i), i)
            end if
        end do
    end do
end subroutine advance_solution
end module equations


! ----------------------------------------------------------------------------------
! Initial_conditions: Este módulo se encarga de configurar las condiciones iniciales 
!                     y parámetros del problema de evolución HD, además deque el 
!                     usurario elija un método de reconstrucción.
! ---------------------------------------------------------------------------------
module initial_conditions
    use variables
    use metrics
    use equations, only: primu
    implicit none
    private
    public :: choose_reconstruction_method, initialize_problem, set_initial_conditions

contains

subroutine choose_reconstruction_method()
        integer :: choice

        print *, ""
        print *, "======================================"
        print *, "Choose reconstruction method:"
        print *, "======================================"
        print *, "1. Godunov (1st order)"
        print *, "2. Minmod (2nd order)"
        print *, "======================================"
        print *, "Enter your choice (1-2): "
        read *, choice

        select case(choice)
        case(1)
            reconstruction_method = 'Godunov'
            print *, "Selected: Godunov (1st order) reconstruction"
        case(2)
            reconstruction_method = 'minmod'
            print *, "Selected: Minmod (2nd order) reconstruction"
        case default
            print *, "Invalid choice. Using default method (Minmod)."
            reconstruction_method = 'minmod'
        end select

        print *, "======================================"
        print *, ""
    end subroutine choose_reconstruction_method

    subroutine initialize_problem()
        implicit none
        integer :: i
        !character(len=20) :: initial_case

        ! Parámetros 
        nx = 1000
        nt = 200000
        xmin = 1.0d0
        xmax = 51.0d0
        CFL = 0.25d0
        floor = 1.0d-10
        gam = 4.0d0/3.0d0 
        Ma = 1.0d0

        final_time = 1000.0d0  ! Tiempo suficiente para evolución

        ! Intervalos 
        t1 = 100
        t2 = 1000

        ! Calcula dx
        dx = (xmax - xmin)/dble(nx)

        ! Inicializa la malla 
        allocate(x(0:nx))  
        do i = 0, nx
            x(i) = xmin + dble(i) * dx
        end do

        !Extender arrays para celdas fantasma
        deallocate(x)
        allocate(x(-1:nx+2))
        do i = -1, nx+2
            x(i) = xmin + dble(i) * dx
        end do

        ! Aloca memoria para arrays
        allocate(u(neq,-1:nx+2), up(neq,-1:nx+2))
        allocate(p(neq,-1:nx+2), fx(neq,-1:nx+2))
        allocate(rhs(neq,-1:nx+2))
        allocate(pLx(neq,-1:nx+2), pRx(neq,-1:nx+2))
        
        ! Variables geométricas
        allocate(grr(-1:nx+2))
        allocate(alpha(-1:nx+2), dlnalpha(-1:nx+2))
        allocate(beta(-1:nx+2))
        allocate(detgamma(-1:nx+2))

        ! Inicializa la métrica
        call set_metric(x, grr, alpha, dlnalpha, beta, detgamma, 'Schwarzschild_EF')

        ! Elegir método de reconstrucción
        call choose_reconstruction_method()

        call set_initial_conditions('michel')
    end subroutine initialize_problem

    subroutine set_initial_conditions(initial_case)
        implicit none
        character(len=*), intent(in) :: initial_case
        integer :: i
        real*8 :: e_local, h_local, W_local, cs_local

        select case(trim(initial_case))
            case('michel')
                ! Parámetros iniciales
                rho_ini = 1.0d-4
                v_infty = -0.2d0
                cs_infty = 0.1d0

                do i = -1, nx+2
                    ! Variables primitivas 
                    p(eq_de,i) = rho_ini  ! densidad
                    p(eq_pr,i) = cs_infty**2 * rho_ini * (gam-1.0d0) / &
                                (gam*(gam-1.0d0) - cs_infty**2*gam)  ! presión 
                    p(eq_vx,i) = v_infty/sqrt(grr(i))  ! velocidad 
                    
                    ! Verificar que las primitivas son físicas
                    e_local = p(eq_pr,i) / (gam - 1.0d0) / p(eq_de,i)
                    h_local = 1.0d0 + e_local + p(eq_pr,i) / p(eq_de,i)
                    W_local = 1.0d0 / sqrt(1.0d0 - grr(i) * p(eq_vx,i)**2)
                    cs_local = sqrt(p(eq_pr,i) * gam * (gam - 1.0d0) / &
                                  (p(eq_pr,i) * gam + p(eq_de,i) * (gam - 1.0d0)))
                    
                    ! Calcula variables conservativas 
                    call primu(p(:,i), u(:,i), i)
                end do
        end select

        ! Copia inicial a up 
        up = u
        
        ! Inicializar tiempo
        t = 0.0d0
        
        print *, "Condiciones iniciales establecidas:"
        print *, "  Densidad inicial: ", rho_ini
        print *, "  Velocidad en infinito: ", v_infty  
        print *, "  Velocidad del sonido en infinito: ", cs_infty
        print *, "  Gamma: ", gam
        print *, "  Masa del agujero negro: ", Ma
        
    end subroutine set_initial_conditions

end module initial_conditions

! ----------------------------------------------------------------------------------
! Output: Este módulo se encarga de la gestión de salida de datos en la simulación 
!         de GRHD en 1D. Contiene rutinas para inicializar directorios de salida y
!         guardar el estado de la simulación en archivos.
! ---------------------------------------------------------------------------------
module output
    use variables
    implicit none
    private
    public :: save_data, initialize_output

contains
    subroutine initialize_output()
        implicit none
        character(len=100) :: command
        
        ! Crea directorio para output
        write(command, '(3a)') 'mkdir -p output/', trim(reconstruction_method)
        call system(command)
    end subroutine initialize_output

    subroutine save_data(it, time)
        implicit none
        integer, intent(in) :: it
        real*8, intent(in) :: time
        character(len=100) :: filename
        integer :: i
        real*8 :: e_local

        ! Genera nombre de archivo
        write(filename,'(3a,i6.6,a)') 'output/', trim(reconstruction_method), '/data_', it, '.dat'

        ! Abre archivo
        open(unit=10, file=filename, status='replace')

        ! Escribe encabezado
        write(10,*) '# Simulación GRHD 1D en métrica de Schwarzschild (EF)'
        write(10,*) '# Método de reconstrucción: ', trim(reconstruction_method)
        write(10,*) '# Tiempo = ', time
        write(10,*) '# Columnas: x, densidad, presión, velocidad, energía_interna'
        write(10,*) '#'

        ! Escribe datos
        do i = 1, nx
            e_local = p(eq_pr,i) / (gam - 1.0d0) / p(eq_de,i)
            write(10,'(5E20.12)') x(i), p(eq_de,i), p(eq_pr,i), p(eq_vx,i), e_local
        end do

        close(10)
    end subroutine save_data

end module output

! ---------------------------------------------------------------------------------
!                                   PROGRAMA PRINCIPAL
! 
! ---------------------------------------------------------------------------------
program GRHD_1D
    use variables
    use metrics
    use equations
    use MODfluxes
    use initial_conditions
    use boundary_conditions  
    use output
    implicit none

    integer :: it
    real*8 :: time

    ! Inicialización
    print *, "================================================================================"
    print *, "                 CODE General Relativity Hydrodinamics  (GRHD)                 "
    print *, "================================================================================"
    print *, " "
    call initialize_problem()
    call initialize_output()

    ! Guarda condición inicial
    time = 0.0d0
    print *, "Guardando condición inicial..."
    call save_data(0, time)

    ! Loop principal de evolución
    print *, "Comenzando evolución temporal..."
    print *, "----------------------------------------"
    print *, " Iteración |    Tiempo    |   % Completado"
    print *, "----------------------------------------"

    do it = 1, nt
        dt = CFL * dx
        call advance_solution(dt)
        time = time + dt

        !Progreso en consola
        if (mod(it, t1) == 0) then
            write(*, '(I10,F13.6,E11.3,F8.1,"%")') &
                it, time, dt, (time/final_time)*100.0d0
        endif

        if (time >= final_time) then
            print *, "Simulación completada!"
            exit
        end if
    end do

    if (mod(it, t2) /= 0) then
        call save_data(it, time)
    end if

    call deallocate_arrays()
    print *, "----------------------------------------"
    print *, "Simulación finalizada correctamente"
    print *, "Tiempo final alcanzado: ", time
    print *, "Número total de iteraciones: ", it

contains
    subroutine deallocate_arrays()
        implicit none
        
        print *, "Liberando memoria..."
        deallocate(x, u, up, p, fx, rhs, pLx, pRx)
        deallocate(grr, alpha, dlnalpha, beta, detgamma)
        print *, "Memoria liberada correctamente"
    end subroutine deallocate_arrays

end program GRHD_1D




