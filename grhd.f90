! ================================================================================
!                             General Relativity Hydrodinamics 
! Este programa resuelve las ecuaciones de GRHD unidimensional usando la métrica
! de Shwarzschild en coordenadas de Eddington Finkestein (EF) usando la descomposición
! 3+1, usando el método de flujos HLLE 
! =================================================================================


! --------------------------------------------------------------------------
! Variables: Este modulo contiene las variables globales que utilizaremos 
!            a lo largo de este código
!  -------------------------------------------------------------------------


module variables
  implicit none
  public

  ! --- Índices de ecuaciones ---
  integer, parameter :: neq = 3   ! equations: 1=de,2=pr,3=vx,4=vy (vy no usada en 1D)
  integer :: eq_de = 1
  integer :: eq_pr = 2
  integer :: eq_vx = 3

  ! --- Dominio numérico---
  integer :: nx
  real*8 :: x_min, x_max, dx, t, dt, final_time
  real*8 :: CFL      ! Número de Courant-Friedrichs-Lewy

  ! --- Arrays principales ---
  real*8, allocatable :: x(:)             ! Posiciones de los centros de las celdas
  real*8, allocatable :: u(:,:)           ! Variables conservativas
  real*8, allocatable :: up(:,:)          ! Variables conservatidas actualizadas 
  real*8, allocatable :: p(:,:)           ! Variables primitivas  
  real*8, allocatable :: fx(:,:)          ! Flujos en dirección x
  real*8, allocatable :: s(:,:)           ! Fuentes 
  real*8, allocatable :: rhs(:,:)         ! Términos de la derecha
  real*8, allocatable :: pL(:,:), pR(:,:) ! Variables primitivas reconstruidas

  ! --- Geometría de Schwarzschild en coordenadas EF ---
  ! real*8, allocatable :: grr(:)        ! Componente radial de la métrica
  ! real*8, allocatable :: alpha(:)      ! Función lapso
  ! real*8, allocatable :: dlnalpha(:)   ! Derivada del logaritmo del lapso
  ! real*8, allocatable :: beta(:)       ! Vector shift
  ! real*8, allocatable :: detgamma(:)   ! Determinante de la métrica espacial
  ! real*8, allocatable :: gamma(:,:,:)  ! Simbolos de Christoffel - gamma^a_{bc}

  !-- Parámetros físicos  ---
  real*8 :: adb_idx  ! Índice adiabático
  real*8 :: g1       ! Factor de compresión adiabática 
  real*8 :: Ma       ! Masa del agujero

  ! --- Método de reconstrucción ---
  character(len=10) :: reconstruction_method ! Método de reconstrucción seleccionado (Godunov, TVD)
  character(len=10) :: tvd_limiter           ! Limiter para el método TVD (Minmod, Superbee, MC)

  ! --- Output ---
  character(len=50) :: output_prefix, output_folder

  ! --- Variables de guardado ---
  real*8 :: save_interval 
  real*8 :: last_save_time = 0.0d0
  real*8 :: next_save_time

  ! --- Variables de condición inicial ---
  character(len=50) :: cond_ini  ! Condición inicial seleccionada

end module variables

! ---------------------------------------------------------------------------------------
! Metrics: en este modulo implementamos la métrica de Schwarzschild para un agujero negro 
!          no rotante en coordenadas EF, implementamos su respectiva descomposición 3+1.
!          Agregamos los simbolos de Christoffel, los tesores energía momento no nulos, 
!          respectivamente. Tambien calculamos los términos fuentes S.
! ---------------------------------------------------------------------------------------
module metrics
  use variables
  implicit none
  private
  public :: calculate_metric, calculate_christoffel_symbols, calculate_stress_energy_tensor
  ! public :: set_metric, calculate_schwarzschild_EF

contains
  ! Elegir la métrica
  ! subroutine set_metric(xi, grr, alpha, dlnalpha, beta, detgamma, metric_type)
  !   implicit none
  !   real*8, intent(in) :: xi(:)
  !   real*8, intent(out) :: grr(:), alpha(:), dlnalpha(:)
  !   real*8, intent(out) :: beta(:), detgamma(:)
  !   character(len=*), intent(in) :: metric_type

  !   select case(trim(metric_type))
  !     case('Schwarzschild_EF')
  !       call calculate_schwarzschild_EF(xi, grr, alpha, dlnalpha, beta, detgamma)
  !     end select
  !   end subroutine set_metric

    ! Métrica de Schwarzschild en coordenadas ED
    subroutine calculate_metric(ri, grr, alpha, dlnalpha, beta, detgamma, gthth, gphph)
      implicit none
      real*8, intent(in) :: ri
      real*8, intent(out), optional :: grr, alpha, dlnalpha, beta, detgamma, gthth, gphph

      ! Métrica de Schwarzschild en coordenadas EF
      if (present(alpha))     alpha = 1.0d0 / sqrt(1.0d0 + 2.0d0 * Ma / ri)      ! alpha
      if (present(dlnalpha))  dlnalpha = Ma / (ri**2 + 2.0d0 * Ma * ri)          ! d(ln(alpha))/dr
      if (present(beta))      beta = 2.0d0 * Ma / ri / (1.0d0 + 2.0d0 * Ma / ri) ! beta^r
      if (present(grr))       grr = 1.0d0 + 2.0d0 * Ma / ri                      ! g_rr = gamma_rr
      if (present(detgamma))  detgamma = ri**4 * grr                             ! Determinante de la métrica espacial sqrt(gamma)
      if (present(gthth))     gthth = ri**2                                      ! g_thth
      if (present(gphph))     gphph = gthth                                      ! g_phph
    end subroutine calculate_metric

    ! ! Métrica de Schwarzschild en coordenadas ED
    ! subroutine calculate_schwarzschild_EF(xi, grr, alpha, dlnalpha, beta, detgamma)
    !   implicit none
    !   real*8, intent(in) :: xi
    !   real*8, intent(out), dimension(:) :: grr, alpha, dlnalpha, beta, detgamma

    !   ! Métrica de Schwarzschild en coordenadas EF
    !   alpha = 1.0d0 / sqrt(1.0d0 + 2.0d0 * Ma / x)     ! alpha
    !   dlnalpha = Ma / (x**2 + 2.0d0 * Ma * x)          ! d(ln(alpha))/dr
    !   beta = 2.0d0 * Ma / x / (1.0d0 + 2.0d0 * Ma / x) ! beta^r
    !   grr = 1.0d0 + 2.0d0 * Ma / x                     ! g_rr = gamma_rr
    !   detgamma = x**4 * grr                            ! Determinante de la métrica espacial sqrt(gamma)
    ! end subroutine calculate_schwarzschild_EF

    ! Componentes del tensor energía momento
    subroutine calculate_stress_energy_tensor(q, ri, T00, Tr0, Trr, Tthth, Tphph)
      implicit none
      real*8, intent(in) :: q(:), ri
      real*8, intent(out) :: T00, Tr0, Trr, Tthth, Tphph
      real*8 :: grr, alpha, beta
      real*8 :: vv, h, W, v

      ! q = Variables primitivas
      ! ri = posición de la celda actual
      ! T** = Componentes del tensor energía-momento

      call calculate_metric(ri, grr=grr, alpha=alpha, beta=beta)

      vv = q(eq_vx)**2 * grr
      h = 1.0d0 + g1 * q(eq_pr) / q(eq_de)
      W = 1.0d0 / sqrt(1.0d0 - vv)
      v = q(eq_vx) - beta / alpha! Velocidad relativa

      ! Componentes del tensor energía-momento
      T00 = (q(eq_de) * h * W**2 - q(eq_pr)) / alpha**2 !u(eq_pr)/alpha^2
      Tr0 = q(eq_de) * h * W**2 / alpha * v + q(eq_pr) * beta / alpha**2
      Trr = q(eq_de) * h * W**2 * v**2 + q(eq_pr) * (1/grr - beta**2/alpha**2)
      Tthth = q(eq_pr)/ri**2
      Tphph = Tthth ! En 1D, Tthth = Tphph
    end subroutine calculate_stress_energy_tensor

    ! Simbolos de Christoffel
    subroutine calculate_christoffel_symbols(ri, chr)
      implicit none
      real*8, intent(in) :: ri
      real*8, intent(out) :: chr(0:3,0:3,0:3) !Notacion gamma^a_{bc}

      ! ri = posición de la celda actual
      ! chr = Simbolos de Christoffel

      chr = 0.0d0 ! Inicializa los Christoffel a cero

      ! Simbolos de Christoffel para la métrica de Schwarzschild en coordenadas EF
      chr(0,0,0) = 2.0d0 * Ma**2 / ri**3
      chr(0,0,1) = Ma * (1.0d0 + 2.0d0 * Ma / ri) / ri**2
      chr(0,1,1) = 2.0d0 * Ma * (1.0d0 + Ma / ri) / ri**2
      chr(1,0,0) = Ma / ri**2 * (1.0d0 - 2.0d0 * Ma / ri)
      chr(1,0,1) = -chr(0,0,0)
      chr(1,1,1) = -chr(0,0,1)
      chr(2,1,2) = 1.0d0 / ri
      chr(3,1,3) = chr(2,1,2)
      chr(0,2,2) = -2.0d0 * Ma
      chr(0,3,3) = chr(0,2,2)
      chr(1,2,2) = 2.0d0 * Ma - ri
    end subroutine calculate_christoffel_symbols

end module metrics

module conditions
  use variables
  use metrics
  implicit none
  private
  public :: initial_conditions, boundary_conditions

  contains

  subroutine initial_conditions()
    integer :: i
    real*8  :: rho_ini  = 1.0d-4 ! Densidad inicial
    real*8  :: v_infty  = -0.2d0 ! Velocidad del flujo asintótico
    real*8  :: cs_infty = 0.1d0  ! Velocidad del sonido asintótica
    real*8  :: grr

    select case(trim(cond_ini))
    case("1")

      do i = 0, nx
        call calculate_metric(x(i), grr=grr)
        p(eq_de,i) = rho_ini
        p(eq_vx,i) = v_infty/sqrt(grr)
        p(eq_pr,i) = rho_ini * cs_infty**2 / (adb_idx - cs_infty**2*g1)
      end do

    ! case("2") ! Bondi Hoyle accretion
    !   p(eq_de) = 1.0d0
    !   do i = 1, nx
    !     p(eq_vx) = v_infty/sqrt(grr)
    !     p(eq_pr) = 1.0d0 * cs_infty**2 / (adb_idx - cs_infty**2*g1)
    !   end do
    end select

  end subroutine initial_conditions

  subroutine boundary_conditions(q)
    implicit none
    real*8, intent(inout) :: q(neq, -1:nx+1)
    integer :: i

    ! Extrapolación lineal en la frontera interna (i=0,-1)
    q(:,0) = 3.0d0*q(:,1) - 3.0d0*q(:,2) + q(:,3)
    ! q(:,-1) = 3.0d0*q(:,0) - 3.0d0*(q(:,1) + q(:,2))

    ! Extrapolación simple en la frontera externa (i=nx+1, nx+2)
    ! q(:,nx+1) = q(:,nx)
    ! q(:,nx+2) = q(:,nx-1)

    ! Condición de frontera en la frontera externa (i=nx)
    ! rhs(Nx) = 0.0d0 -- Inyección permanente

  end subroutine boundary_conditions

end module conditions

module initialization
  use variables
  use metrics
  use conditions
  implicit none

contains

  subroutine choose_reconstruction_method()
    integer :: choice, choice_limiter

    ! print *, "Choose reconstruction method:"
    ! print *, "1. Godunov (1st order)"
    ! print *, "2. TVD (2nd order)"
    ! read *, choice
    choice = 1 ! Descomentar para entrada manual

    select case(choice)
    case(1)
      reconstruction_method = 'Godunov'
    case(2)
      reconstruction_method = 'TVD'
      ! print *, "Choose TVD limiter:"
      ! print *, "1. Minmod"
      ! print *, "2. Superbee"
      ! print *, "3. Monotonized Centered"
      ! read *, choice_limiter
      choice_limiter = 1

      select case(choice_limiter)
      case(1)
        tvd_limiter = 'minmod'
      case(2)
        tvd_limiter = 'superbee'
      case(3)
        tvd_limiter = 'mc'
      case default
        print *, "Invalid choice. Using minmod."
        tvd_limiter = 'minmod'
      end select
    case default
      print *, "Invalid choice. Using default method (TVD - Minmod)."
      reconstruction_method = 'TVD'
      tvd_limiter = 'minmod'
    end select

    ! print *, "Selected reconstruction method: ", trim(reconstruction_method)
  end subroutine choose_reconstruction_method


  subroutine initialize_problem()
    integer :: i
    character(len=50) :: case_choice

    ! print *, "Select test case:"
    ! print *, "1. Michel Accretion"
    ! print *, "2. Bondi Hoyle"
    ! read *, case_choice
    case_choice = "1" ! Descomentar para entrada manual
    cond_ini = case_choice

    select case(cond_ini)
    case("1")
      call init_MichelAccretion()
    ! case("2")
    !   call init_BondiHoyle()
    case default
      print *, "Invalid selection. Exiting."
      stop
    end select

    g1 = adb_idx / (adb_idx - 1d0)

    dx = (x_max - x_min) / nx

    ! Inicializa la malla
    allocate(x(0:nx)) ! Puntos centrales de las celdas

    !$OMP PARALLEL DO PRIVATE(i)
    do i = 0, nx
      x(i) = x_min + dx*i 
    end do
    !$OMP END PARALLEL DO

    ! Asigna memoria para arrays
    allocate(u(neq,-1:nx+1))
    allocate(up(neq,-1:nx+1))
    allocate(p(neq,-1:nx+1))
    allocate(fx(neq,-1:nx+1))
    allocate(s(neq,-1:nx+1))
    allocate(rhs(neq,-1:nx+1))
    allocate(pL(neq,-1:nx+1))
    allocate(pR(neq,-1:nx+1))

    ! ! Inicializa la métrica
    ! call set_metric(x, grr, alpha, dlnalpha, beta, detgamma, 'Schwarzschild_EF')

    ! Create output folder if it doesn't exist
    call system('mkdir -p ' // trim(output_folder))
    ! Create subfolders for each reconstruction method
    call system('mkdir -p ' // trim(output_folder) // '/Godunov')
    call system('mkdir -p ' // trim(output_folder) // '/TVD')

    ! Set initial conditions
    call initial_conditions()

    ! Set reconstruction method
    call choose_reconstruction_method()

  end subroutine initialize_problem

  subroutine init_MichelAccretion()
    ! Initialize parameters for Michel Accretion
    nx = 1000
    x_min = 1.0d0
    x_max = 51.0d0
    adb_idx = 4.0d0/3.0d0
    Ma = 1.0d0
    final_time = 500.0d0
    output_prefix = 'MA'
    output_folder = 'MichelAccretion/'
    CFL = 0.25d0
    save_interval = 50.0d0
  end subroutine init_MichelAccretion

  ! subroutine init_BondiHoyle()
  !   ! Initialize parameters for Bondi Hoyle accretion
  !   nx = 100
  !   x_min = 1.0d0
  !   x_max = 10.0d0
  !   adb_idx = 5.0d0/3.0d0
  !   final_time = 1.0d0
  !   output_prefix = 'BondiHoyle'
  !   output_folder = 'output/BondiHoyle'
  !   CFL = 0.5d0
  !   save_interval = 0.1d0
  ! end subroutine init_BondiHoyle
end module initialization

module reconstruction
  use variables
  implicit none
  private
  public :: reconstruct_variables

  contains

  real*8 function minmod(a,b)
    implicit none
    real*8, intent(in) :: a, b
    minmod = 0.5d0 * (sign(1.0d0,a) + sign(1.0d0,b)) * min(abs(a),abs(b))
  end function minmod

  real*8 function minmod3(a,b,c)
    implicit none
    real*8, intent(in) :: a,b,c
    minmod3 = sign(1.0d0,a) * max(0.0d0, min(abs(a), sign(1.0d0,a)*b, sign(1.0d0,a)*c))
  end function minmod3

  ! Add to equations_1d module
  function tvd_slope(a, b, limiter) result(slope)
    real*8, intent(in) :: a, b
    character(len=*), intent(in) :: limiter
    real*8 :: slope

    select case(trim(limiter))
    case('minmod')
      slope = minmod(a, b)
    case('superbee')
      if (a*b <= 0.0d0) then
        slope = 0.0d0
      else
        slope = sign(1.0d0, a) * max(min(2.0d0*abs(a), abs(b)), min(abs(a), 2.0d0*abs(b)))
      end if
    case('mc')
      slope = minmod3(2.0d0*b,0.5d0*(a+b), 2.0d0*a)
    case default
      slope = minmod(a, b)
    end select
  end function tvd_slope

  subroutine reconstruct_variables(q, qL, qR)
    implicit none
    real*8, intent(in)                           :: q(neq,-1:nx+1)
    real*8, intent(out), dimension(neq, -1:nx+1) :: qL, qR
    integer                                      :: i, k
    real*8                                       :: smaxL, sminL, smaxR, sminR

    select case(trim(reconstruction_method))
    case('Godunov')
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 0, nx-1
        qL(:,i) = q(:, i )
        qR(:,i) = q(:,i+1)
        ! print*, qL(:,i), i
      end do
      !$OMP END PARALLEL DO

    case('TVD')
      !$OMP PARALLEL DO PRIVATE(k,i,sminL,smaxL,sminR,smaxR)
      ! q(:,-1) = 0.0d0
      do k = 1, neq
        do i = 0, nx-1
          ! Left state reconstruction
          sminL = q(k,i) - q(k,i-1)
          smaxL = q(k,i+1) - q(k,i)
          qL(k,i) = q(k,i) + 0.5d0 * tvd_slope(sminL, smaxL, tvd_limiter)

          ! Right state reconstruction
          sminR = q(k,i+1) - q(k,i)
          smaxR = q(k,i+2) - q(k,i+1)
          qR(k,i) = q(k,i+1) - 0.5d0 * tvd_slope(sminR, smaxR, tvd_limiter)
        end do
      end do
      !$OMP END PARALLEL DO

    end select

  end subroutine reconstruct_variables

end module reconstruction

module equations
  use variables
  use metrics
  implicit none
  private
  public :: primu, uprim, primfx, source

contains

  subroutine primu(q,y,xi) 
    implicit none
    real*8, intent(in) :: q(:), xi
    real*8, intent(out) :: y(:)
    real*8 :: vv, h, W, grr, detgamma

    ! q = Variables primitivas
    ! y = Variables conservativas
    ! xi = posición de la celda actual

    call calculate_metric(xi, grr=grr, detgamma=detgamma)

    vv = q(eq_vx)**2*grr
    h = 1d0 + g1 * q(eq_pr)/q(eq_de)
    W = 1d0 / sqrt(1d0 - vv)

    y(eq_de) = q(eq_de) * W 
    y(eq_pr) = q(eq_de) * h * W**2 - q(eq_pr)
    y(eq_vx) = q(eq_de) * h * W**2 * q(eq_vx) * grr

    y = y*sqrt(detgamma)
  end subroutine primu

  subroutine uprim(y,q,xi)
    implicit none
    real*8, intent(in) :: y(:), xi
    real*8, intent(out) :: q(:)
    real*8 :: lor, V, fn, dV
    real*8 :: grr, detgamma
    real*8 :: floor = 1d-16

    ! y = Variables conservativas
    ! q = Variables primitivas
    ! xi = posición de la celda actual

    lor = rtsafe(y,xi)
    call VdV(lor, y, V, dV)
    call calculate_metric(xi, grr=grr, detgamma=detgamma)

    q(eq_de) = y(eq_de) / lor
    q(eq_vx) = y(eq_vx) / (V * grr)
    q(eq_pr) = V - y(eq_pr)

    q(eq_de:eq_pr) = q(eq_de:eq_pr) / sqrt(detgamma)

    q(eq_de) = max(q(eq_de), floor)
    q(eq_pr) = max(q(eq_pr), floor)

    if (q(eq_de) <= 0d0 .or. q(eq_pr) <= 0d0 .or. abs(q(eq_vx)) >= 1d0) then
      print *, "Error: Valores no físicos en uprim", q, lor
      stop
    end if 
  end subroutine uprim

  subroutine VdV(W,y,V,dV)
    real*8, intent(in) :: W, y(:)
    real*8, intent(out) :: V, dV
    real*8 :: denom

    ! W = Factor de Lorentz 
    ! y = Variables conservativas
    ! V = V(W)
    ! dV = dV/dW

    denom = g1 * W**2 - 1d0
    V = (g1 * W**2 * y(eq_pr) - y(eq_de) * W) / denom
    dV = (g1 * W * (y(eq_de) * W - 2d0 * y(eq_pr)) + y(eq_de)) / denom**2
  end subroutine VdV

  subroutine funcd(W, y, x_i, fn, df)
    real*8, intent(in) :: W, y(:), x_i
    real*8, intent(out) :: fn, df
    real*8 :: V, dV, grr

    ! W = Factor de Lorentz
    ! y = Variables conservativas
    ! x_i = posicion de la celda actual
    ! fn = f(W)
    ! df = f'(W)

    call calculate_metric(x_i, grr=grr)
    call VdV(W,y(:),V,dV)

    fn = V**2*(W**2-1d0)/W**2-y(eq_vx)**2/grr
    df = 2*V*(dV*(W**3 - W) + V)/(W**3)
  end subroutine funcd

  function rtsafe(y,x_i) result(root)
    implicit none
    real*8, intent(in) :: y(:), x_i
    real*8 :: x1 = 1.0d0, x2 = 1.0d3, xacc = 1.0d-8
    real*8 :: root
    integer, parameter :: MAXIT = 100000
    integer :: j
    real*8 :: df, diffx, diffxold, f, fh, fl, temp, xh, xl
    
    call funcd(x1, y, x_i, fl, df)
    call funcd(x2, y, x_i, fh, df)

    if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) then
      print *, fl, fh, 'root must be bracketed in rtsafe'
      stop
    end if
    
    if (fl == 0.0) then
      root = x1
      return
    else if (fh == 0.0) then
      root = x2
      return
    end if
    
    if (fl < 0.0) then
      xl = x1
      xh = x2
    else
      xh = x1
      xl = x2
    end if
    
    root = 0.5*(x1 + x2)
    diffxold = abs(x2 - x1)
    diffx = diffxold
    call funcd(root, y, x_i, f, df)
    
    do j = 1, MAXIT
      if (((root - xh)*df - f)*((root - xl)*df - f) > 0.0 .or. &
          abs(2.0*f) > abs(diffxold*df)) then
        diffxold = diffx
        diffx = 0.5*(xh - xl)
        root = xl + diffx
        if (xl == root) return 
      else
        diffxold = diffx
        diffx = f/df
        temp = root
        root = root - diffx
        if (temp == root) return 
      end if
      
      if (abs(diffx) < xacc) return 
      
      call funcd(root, y, x_i, f, df)
      
      if (f < 0.0) then
        xl = root
      else
        xh = root
      end if

    end do

    print *, 'rtsafe exceeding maximum iterations', abs(diffx), root
    stop
  end function rtsafe

  subroutine primfx(q,f,ri)
    implicit none
    real*8, intent(in) :: q(:), ri
    real*8, intent(out) :: f(:)
    real*8 :: vv, h, W, v_ef
    real*8 :: grr, alpha, beta, detgamma

    ! q = Variables primitivas
    ! f = Flujos
    ! ri = posición de la interfaz actual

    call calculate_metric(ri, grr=grr, alpha=alpha, beta=beta, detgamma=detgamma)

    vv = q(eq_vx)**2*grr
    h = 1d0 + g1 * q(eq_pr)/q(eq_de)
    W = 1d0 / sqrt(1d0 - vv)
    v_ef = alpha*q(eq_vx)- beta ! alpha*v^r - beta^r

    f(eq_de) = q(eq_de) * W * v_ef
    f(eq_pr) = q(eq_de) * h * W**2 * v_ef + q(eq_pr)*beta
    f(eq_vx) = q(eq_de) * h * W**2 * q(eq_vx) * grr * v_ef + alpha*q(eq_pr)

    f = f*sqrt(detgamma)

  end subroutine primfx

  subroutine source(q, src, ri)
  implicit none
  real*8, intent(in) :: q(:), ri
  real*8, intent(out) :: src(:)
  real*8 :: T00, Tr0, Trr, Tthth, Tphph
  real*8 :: vv, h, W, v_ef, chr(0:3,0:3,0:3)
  real*8 :: aux1, aux2
  real*8 :: grr, alpha, beta, dlnalpha, detgamma, gthth, gphph

  ! q = Variables primitivas
  ! ri = posición de la celda actual
  ! src = Términos fuente

  call calculate_stress_energy_tensor(q, ri, T00, Tr0, Trr, Tthth, Tphph)
  call calculate_christoffel_symbols(ri, chr)
  call calculate_metric(ri, grr=grr, alpha=alpha, beta=beta, dlnalpha=dlnalpha, detgamma=detgamma, gthth=gthth, gphph=gphph)

  aux1 = -alpha**2 + grr*beta**2

  src(eq_de)  = 0.0d0
  src(eq_vx)  = T00 * (aux1*chr(0,0,1) + beta*grr*chr(1,0,1)) &
              + Tr0 * (aux1*chr(0,1,1) + beta*grr*(chr(0,0,1) + chr(1,1,1)) + grr*chr(1,0,1)) &
              + Trr * (grr*chr(1,1,1) + beta*grr*chr(0,1,1)) &
              + Tthth *gthth*chr(2,1,2) + Tphph * gphph*chr(3,1,3)
  src(eq_pr)  = alpha * (-T00*chr(0,0,0) + Tr0 * (dlnalpha - 2.0d0*chr(0,0,1)) &
              - (Trr * chr(0,1,1) + Tthth * chr(0,2,2) + Tphph * chr(0,3,3)) )
  src = src*alpha*sqrt(detgamma)
  end subroutine source

end module equations

module flujos
  use equations
  use variables
  use metrics
  implicit none
  private
  public :: fluxes

contains

  subroutine fluxes(ql, qr, rl, rr, f)
    real*8, intent(in) :: ql(:), qr(:), rl, rr
    real*8, intent(out) :: f(:)
    real*8 :: ul(neq), ur(neq), fl(neq), fr(neq)
    real*8 :: ap, am, cflmin, cflmax, cfrmin, cfrmax

    ! ql = Variables primitivas a la izquierda de la interfaz
    ! qr = Variables primitivas a la derecha de la interfaz
    ! ri = posición de la interfaz actual
    ! f  = Flujos en la interfaz

    call primfx(ql, fl, rl)
    call primfx(qr, fr, rr)

    call primu(ql, ul, rl)
    call primu(qr, ur, rr)

    call soundcfp(ql, rl, cflmin, cflmax)
    call soundcfp(qr, rr, cfrmin, cfrmax)

    ap = max(0d0, cflmax, cfrmax)
    am = min(0d0, cflmin, cfrmin)

    f = (ap*fl - am*fr + ap*am*(ur - ul) ) / (ap - am)
  end subroutine fluxes

  subroutine soundcfp(q, ri, xmin, xmax)
    implicit none
    double precision, intent(in) :: q(:), ri
    double precision, intent(out) :: xmin,xmax
    double precision :: cs, cs2, aux1, aux2
    double precision :: eig1, eig2, eig3, h, vv
    double precision :: grr, alpha, beta

    ! q = Variables primitivas
    ! ri = posición de la interfaz actual
    ! xmin, xmax = velocidades características

    call calculate_metric(ri, grr=grr, alpha=alpha, beta=beta)

    h = 1d0 + g1*q(eq_pr)/q(eq_de)
    vv = q(eq_vx)**2*grr
    cs2 = adb_idx * q(eq_pr)/q(eq_de)/h
    cs = sqrt(cs2)

    aux1 = alpha/(1d0 - vv*cs2)
    aux2 = cs * sqrt( (1d0 - vv) * (alpha/grr/aux1 - q(eq_vx)**2*(1d0 - cs2)) )

    eig1 = alpha*q(eq_vx) - beta
    eig2 = aux1*(q(eq_vx)*(1d0 - cs2) + aux2) - beta
    eig3 = aux1*(q(eq_vx)*(1d0 - cs2) - aux2) - beta

    xmax = max(0d0, eig1, eig2, eig3)
    xmin = min(0d0, eig1, eig2, eig3)

  end subroutine soundcfp

end module flujos

module output
  use variables
  implicit none
  private
  public :: save_data, save_data_at_time

  contains

  subroutine save_data(filename)
    character(len=*), intent(in) :: filename
    integer :: eq, i, unit
    character(len=100) :: full_filename

    full_filename = trim(output_folder) // trim(reconstruction_method) // '/' // trim(filename)

    print *, 'writing output file:', trim(filename)
    unit = 10

    open(unit, file=full_filename, status='replace', action='write', form='formatted')

    do eq = 1, neq
      do i = 0, nx
        write(unit, *) p(eq, i)
      end do
    end do

    close(unit)

  end subroutine save_data

  subroutine save_data_at_time(time)
    real*8, intent(in) :: time
    character(len=50) :: filename

    write(filename, '(A,F0.1,".dat")') trim(output_prefix), time
    call save_data(filename)
    print *, 'Datos guardados para t =', time
  end subroutine save_data_at_time

end module output

program grhd
  use variables
  use metrics
  use initialization
  use equations
  use flujos
  use reconstruction
  use output
  use omp_lib
  implicit none

  real*8 :: start_time, end_time, integration_time
  integer :: n_steps, i, rk, k

  integration_time = 0d0
  n_steps = 0
  start_time = omp_get_wtime()

  call initialize_problem()

  !$OMP PARALLEL DO PRIVATE(i)
  do i = 0, nx
    call primu(p(:,i), up(:,i), x(i))
  end do
  !$OMP END PARALLEL DO

  call save_data_at_time(0.0d0)
  next_save_time = save_interval
  dt = CFL * dx

  do while (integration_time < final_time)
    n_steps = n_steps + 1

    if (mod(n_steps,1) == 0) then
      print *, 'Step ', n_steps, ' - Time ', integration_time, ' (', real(integration_time/final_time*100), '%)'
    end if

    if (integration_time >= next_save_time) then
      call save_data_at_time(next_save_time)
      last_save_time = next_save_time
      next_save_time = next_save_time + save_interval
    end if

    u = up

  integration_time = integration_time + dt

    do rk = 1, 3
      call reconstruct_variables(p, pL, pR)

      !$OMP PARALLEL DO PRIVATE(i)
      do i = 0, nx-1
        call fluxes(pL(:,i), pR(:,i), x(i), x(i+1), fx(:,i))
        call source(p(:,i), s(:,i), x(i))
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1, nx-1
        rhs(:,i) = -(fx(:,i) - fx(:,i-1))/ dx + s(:,i)
      end do
      !$OMP END PARALLEL DO

      ! Inyección permanente en la frontera externa
      rhs(:,nx) = 0.0d0

      if (rk == 1) then
        up = u + rhs * dt
      else if (rk == 2) then
        up = 0.75d0 * u + 0.25d0 * (up + rhs * dt)
      else
        up = (u + 2.0d0*(up + rhs * dt))/3.0d0
      end if

      call boundary_conditions(up)

      !$OMP PARALLEL DO PRIVATE(i)
      do i = 0, nx
        call uprim(up(:,i), p(:,i), x(i))
      end do
      !$OMP END PARALLEL DO

    end do
  end do

  call save_data_at_time(integration_time)

  end_time = omp_get_wtime()

  ! Imprimir el tiempo total de ejecución
  print *, 'Execution time:', end_time - start_time, 'seconds'

end program grhd