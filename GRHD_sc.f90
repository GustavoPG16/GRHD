!=============================================================================!
module variables
   implicit none
   public

   integer :: eq_pr = 1
   integer :: eq_de = 2
   integer :: eq_vx = 3
   integer :: eq_vy = 4
   integer :: eq_vz = 5
   integer, parameter :: neq = 5

   integer, parameter :: nx = 10
   integer, parameter :: ny = 10
   integer, parameter :: nz = 10
   real*8, parameter :: pi = acos(-1.0d0)

   real*8 :: x_min = 0d0
   real*8 :: y_min = -pi
   real*8 :: z_min = 0d0

   real*8 :: x_max = 1d1
   real*8 :: y_max = pi
   real*8 :: z_max = 2*pi
   
   real*8 :: final_time = 10d0

   real*8 :: u (neq, 0:nx+1, 0:ny+1, 0:nz+1)
   real*8 :: up(neq, 0:nx+1, 0:ny+1, 0:nz+1)
   real*8 :: p (neq, 0:nx+1, 0:ny+1, 0:nz+1)
   real*8 :: s (neq, 0:nx+1, 0:ny+1, 0:nz+1)

   real*8, parameter :: adb_idx = 4d0/3d0
   real*8, parameter :: M = 1d0

   character(len=50) :: fileini = 'test0.dat' ! Condiciones iniciales
   character(len=50) :: filefin = 'test1.dat' ! Valores finales

end module variables
!=============================================================================!

!=============================================================================!
module tensor_calculations

  use variables
  implicit none
  public
  
  real*8 :: g_00, g_01, g_02, g_03, g_11, g_12, g_13, g_22, g_23, g_33
  real*8 :: g00, g01, g02, g03, g11, g12, g13, g22, g23, g33
  real*8 :: alpha, beta(3), sqrt_g, sqrt_gamma, d_alpha(4), d_beta(4)
  real*8 :: Chr0_00, Chr0_01, Chr0_11, Chr0_22, Chr0_33
  real*8 :: Chr1_01, Chr1_11, Chr1_22, Chr2_12, Chr3_31, Chr3_32
  real*8 :: U0, U1, U2, U3, T00, T01, T02, T11, T12, T22, T33

contains
  subroutine tensor(p, x, y, z)
    implicit none
    real*8, intent(in) :: p(:), x, y, z
    real*8 :: vv, h, W

    ! Cálculos
    vv = sum(p(eq_vx:eq_vz)**2)  ! Módulo al cuadrado de la velocidad total
    h = 1d0 + (adb_idx)/(adb_idx - 1d0) * p(eq_pr)/p(eq_de)  ! Entalpía
    W = 1d0 / sqrt(1d0 - vv)  ! Factor de Lorentz

    g_00 = - (1 - 2*M/x)
    g_01 = 2*M/x
    g_02 = 0
    g_03 = 0
    g_11 = 1 + 2*M/x
    g_12 = 0
    g_13 = 0
    g_22 = x**2
    g_23 = 0
    g_33 = x**2 * sin(y)**2

    alpha = 1.0 / sqrt(1 + 2*M/x)
    beta(eq_vx) = (2*M/x)/(1 + 2*M/x)
    beta(eq_vy:eq_vz) = 0.0

    sqrt_g = x**2 * sin(y)
    sqrt_gamma = sqrt_g/alpha

    dr_alpha = M/(x**2*alpha**3)

    Chr_000 = 2*M**2/x**3
    Chr_001 = (x + 2*M)*M/x**3
    Chr_011 = (2*M)*(M + x)/x**3
    Chr_022 = -2*M
    Chr_033 = -(2*M)*sin(y)**2
    Chr_101 = -(2*M**2)/x**3
    Chr_111 = -Chr_001
    Chr_122 = -x + 2*M
    Chr_212 = 1/x
    Chr_331 = Chr_212
    Chr_332 = cos(y)/sin(y)

    U_0 = W / alpha
    U_1 = U_0 * (alpha * p(eq_vx) - beta(eq_vx))
    U_2 = U_0 * (alpha * p(eq_vy) - beta(eq_vy))
    U_3 = U_0 * (alpha * p(eq_vz) - beta(eq_vz))

    T_00 = p(eq_de) * h * U_0**2 + p(eq_pr) / g_00
    T_01 = p(eq_de) * h * U_0 * U_1 + p(eq_pr) * g_01
    T_02 = p(eq_de) * h * U_0 * U_2 
    T_11 = p(eq_de) * h * U_1**2 + p(eq_pr) / g_11
    T_12 = p(eq_de) * h * U_1 * U_2 
    T_22 = p(eq_de) * h * U_2**2 + p(eq_pr) / g_22
    T_33 = p(eq_de) * h * U_3**2 + p(eq_pr) / g_33

  end subroutine tensor
end module tensor_calculations
!=============================================================================!

!=============================================================================!
 module equations
  use tensor_calculations
  use variables

   implicit none
   private
   public :: primfx, source, primu, uprim

 contains

   ! find the flux in the eq direction f as a function of the primitive variables p
   subroutine primfx(p,f)
     implicit none
     real*8, intent(in)  :: p(:)
     real*8, intent(out) :: f(:)
     real*8 :: h, W, vv, v_ef(eq_vx:eq_vz)

     vv = sum(p(eq_vx:eq_vz)**2) ! Modulo al cuadrado de la velocidad total
     h = 1d0 + adb_idx/(adb_idx - 1d0) * p(eq_pr)/p(eq_de) ! Entalpia
     W = 1d0 / sqrt(1d0 - vv) ! Factor de Lorentz
     v_ef(eq_vx:eq_vz) = p(eq_vx:eq_vz) - beta(eq_vx:eq_vz)/alpha ! Variable auxiliar

     f(eq_pr) = sqrt_g * ((p(eq_de) * h * W**2 - p(eq_pr) - p(eq_de) * W) * v_ef(eq_vx) + p(eq_pr) * p(eq_vx))
     f(eq_de) = sqrt_g * p(eq_de) * W * v_ef(eq_vx)
     f(eq_vx) = sqrt_g * (p(eq_de) * h * W**2 * p(eq_vx) * g_11 * v_ef(eq_vx) + p(eq_pr))
     f(eq_vy) = sqrt_g * p(eq_de) * h * W**2 * p(eq_vy) * g_22 * v_ef(eq_vy)
     f(eq_vz) = sqrt_g * p(eq_de) * h * W**2 * p(eq_vz) * g_33 * v_ef(eq_vz)
   end subroutine primfx

  ! source terms for the equations
  subroutine source(p, s)
    implicit none
    real*8, intent(in)  :: p(:)
    real*8, intent(out) :: s(:)
    s(eq_pr) = sqrt_g * (T_01*(dr_alpha - 2*alpha*Chr_001) - alpha*(Chr_000*T_00 + Chr_011*T_11 + Chr_022*T_22 + Chr_033*T_33) )
    s(eq_de) = 0
    s(eq_vx) = sqrt_g * (T_00 * (g_00*Chr_001 + g_01*Chr_101) + T_01*(g_00*Chr_011 + g_01*(Chr_111 + Chr_001) + g_11*Chr_101) & 
        + T_11*(g_01*Chr_011 + g_11*Chr_111) + T_22*Chr_212*g_22 + T_33*Chr_331*g_33)
    s(eq_vy) = sqrt_g * ( T_02 * (Chr_022*g_00 + Chr_122*g_01) + T_12*(Chr_022*g_01 + Chr_122*g_11 + Chr_212*g_22) & 
        + T_33*Chr_332*g_33)
    s(eq_vz) = 0 
  end subroutine source

   ! find the conservative variables u from the primitive variables p
   subroutine primu(p,u)
     implicit none
     real*8, intent(in)  :: p(:)
     real*8, intent(out) :: u(:)
     real*8 :: h, W, vv
     integer :: i
     vv = p(eq_vx)**2*g_11 + p(eq_vy)**2*g_22 + p(eq_vz)**2*g_33
     h = 1d0 + adb_idx/(adb_idx - 1d0) * p(eq_pr)/p(eq_de)
     W = 1d0 / sqrt(1d0 - vv)
     u(eq_pr) = sqrt_gamma * (p(eq_de) * h * W**2 - p(eq_pr) - p(eq_de) * W)
     u(eq_de) = sqrt_gamma * p(eq_de) * W
     u(eq_vx) = sqrt_gamma * p(eq_de) * h * W**2 * p(eq_vx) * g_11
     u(eq_vy) = sqrt_gamma * p(eq_de) * h * W**2 * p(eq_vy) * g_22
     u(eq_vz) = sqrt_gamma * p(eq_de) * h * W**2 * p(eq_vz) * g_33
     print *, W
     stop
   end subroutine primu

   ! find the primitive variables p from the conservative variables u
   subroutine uprim(u,p)
     ! calculates the primitive variables as a function of the conservative variables
     double precision, intent(in)  :: u(:)
     double precision, intent(out) :: p(:)
     double precision :: W, p_guess, h, m2
     m2 = u(eq_vx)**2/g_11 + u(eq_vy)**2/g_22 + u(eq_vz)**2/g_33
     call newrap(m2,u(eq_de),u(eq_pr),p_guess,W)

     p(eq_pr) = p_guess
     p(eq_de) = u(eq_de)/(W*sqrt_gamma)
     h = 1d0 + adb_idx/(adb_idx - 1d0)*p(eq_pr)*p(eq_de)
     p(eq_vx) = u(eq_vx)/(u(eq_de)*h*W)/g_11
     p(eq_vy) = u(eq_vy)/(u(eq_de)*h*W)/g_22
     p(eq_vz) = u(eq_vz)/(u(eq_de)*h*W)/g_33
     
   end subroutine uprim

   subroutine newrap(m2,d,tau,p_guess,W)
     double precision, intent(in) :: m2, d, tau
     double precision, intent(out) :: p_guess, W
     double precision :: f, dfdp, w_prime, aux
     double precision, parameter :: eps = 1d-12
     integer :: i
     ! start with p_guess = 1, W = 1

    !  print *, m2
    !  stop

     p_guess = 1d0
     W  = 1d0
    !  aux = tau + sqrt_gamma*p_guess + d
    !  W = abs(aux/sqrt(aux**2 - m2))

    !  print *, W
    !  stop

     p_guess = (adb_idx - 1d0)*(tau + d*(1d0-W) + sqrt_gamma*(1d0-W**2))/(sqrt_gamma*W**2)

     do i = 1, 100000000
       aux = tau + sqrt_gamma*p_guess + d
       W = abs(aux/sqrt(aux**2 - m2))
       w_prime = -(sqrt_gamma * m2 * W**3)/(aux)**3

       f = p_guess - (adb_idx - 1d0)*(tau + d*(1d0-W) + sqrt_gamma*p_guess*(1d0-W**2))/(sqrt_gamma * W**2)
       dfdp = -((adb_idx - 1d0)*(d*W- 2*aux)*w_prime - W*sqrt_gamma*(adb_idx*(W**2 - 1d0) + 1d0))/(sqrt_gamma*W**3)

       p_guess = p_guess - f / dfdp                    ! Newton-raphson iteration
       if(abs(f).lt.eps) then
        print *, 'Acabe en',i,'pasos',W
        stop
        ! print *, p_guess
        return
       end if
     end do

     print *, W
     stop
   end subroutine newrap

 end module equations
!=============================================================================!

!=============================================================================!
module MOD_fluxes

   use equations
   implicit none
   private
   public :: fluxes

   contains

   subroutine fluxes(pl,pr,f)
     real*8, intent(in)  :: pl(:), pr(:)
     real*8, intent(out) :: f(:)
     real*8 :: ul(size(pl)), ur(size(pl)), fl(size(pl)), fr(size(pl))  
     real*8 :: ap, am, cflmin, cflmax, cfrmin, cfrmax
     call primfx(pl,fl)
     call primfx(pr,fr)
     call primu(pl,ul)
     call primu(pr,ur)
     call soundcfp(pl,cflmin,cflmax)
     call soundcfp(pr,cfrmin,cfrmax)
     ap = max(0d0,cflmax,cfrmax)
     am = min(0d0,cflmin,cfrmin)
     f(:) = (ap*fl(:)-am*fr(:)+ap*am*(ur(:)-ul(:)))/(ap-am)   
   end subroutine fluxes

   subroutine soundcfp(p,xmin,xmax) ! Modificar esta subrutina
    implicit none
    double precision, intent(in)  :: p(5)
    double precision, intent(out)    :: xmin,xmax
  !   double precision :: lor, gg, h, cs2, vv
  !   double precision :: c(5), sqr
  !   vv  = p(eq_vr)**2
  !   lor = sqrt(1d0 - vv)

  !   h = 1d0 + 4d0 * p(eq_pr)/p(eq_de)
  !   cs2 = gamma*p(eq_pr)/p(eq_de)/h

  !   c(1) = 1d0/cs2 - vv
  !   c(2) = p(eq_vr)*lor
  !   c(3) = p(eq_vr)**2 - 1d0
  !   sqr = sqrt(c(2)**2-c(1)*c(3))
  !   xmax = (-c(2) + sqr)/c(1)*lor + p(eq_vr)
  !   xmin = (-c(2) - sqr)/c(1)*lor + p(eq_vr)

  !   xmax = min(xmax, 1d0)
  !   xmin = max(xmin,-1d0)
   end subroutine soundcfp

end module MOD_fluxes
!=============================================================================!

!=============================================================================!
program nada
  use variables
  use tensor_calculations
  use equations
  implicit none

  real*8 :: dt, dx, dy, dz, x(nx), y(ny), z(nz), a(0:nx), b(0:ny), c(0:nz)
  real*8 :: integration_time = 0d0
  integer :: n_steps = 0, i, j, k

  ! size of the cells
  dx = abs(x_max-x_min)/ nx
  dy = abs(y_max-y_min)/ ny
  dz = abs(z_max-z_min)/ nz

  ! coordinates at the center of each cell
  do i = 1, nx
    x(i) = x_min + dx*(i-0.5d0)
  end do
  do j = 1, ny
    y(j) = y_min + dy*(j-0.5d0)
  end do
  do k = 1, nz
    z(k) = z_min + dz*(k-0.5d0)
  end do
  ! coordinates at the intercells of each cell
  do i = 0, nx
    a(i) = x_min + dx*i
  end do
  do j = 0, ny
    b(j) = y_min + dy*j
  end do
  do k = 0, nz
    c(k) = z_min + dz*k
  end do

  ! initial conditions
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        p(eq_pr,i,j,k) = 0.1d0 ! pressure
        p(eq_de,i,j,k) = 1d0   ! density
        p(eq_vx,i,j,k) = 0.4d0  ! velocity x
        p(eq_vy,i,j,k) = 0d0  ! velocity y
        p(eq_vz,i,j,k) = 0d0  ! velocity z
        !  if (x(i).gt.x(2)) then
        !   p(eq_de,i) = 1d0
        ! end if
      end do
    end do
  end do

  ! print *, 'antes',p(eq_vx, :, 1, 1)
  ! stop
  
  ! compute the conservative variables
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        call tensor(p(:,i,j,k),x(i),y(j),z(k))
        call primu(p(:,i,j,k), up(:,i,j,k))
        ! call source(p(:,i,j,k), s(:,i,j,k))
        call uprim(up(:,i,j,k), p(:,i,j,k)) ! Para verificar cuando se pasa de conservadas a primitivas
      end do
    end do
  end do

  print *, 'despues',p(eq_vx, :, 1, 1)
  stop

  ! find the timestep
  dt = 0.8*min(dx, dy,dz)
  ! print *, dt, dx, dy, dz

  do while((integration_time.lt.final_time))
    n_steps = n_steps + 1

    ! print the iteration number
    if (mod(n_steps,1).eq.0) then 
      print *, 'Paso número',n_steps,'-',real(integration_time/final_time*1d2),'%'
    end if

    ! copy the update conserved variables into the initial variables (for the next step)
    u = up

    ! pupdate the iteration time
    integration_time = integration_time + dt

    ! boundary conditions ! Duda sobre los valores de las "esquinas", es decir, fronteras que sea de x e y o casos similares
    p(:, 0, :, :) = p(:, 1, :, :)
    p(:, :, 0, :) = p(:, :, 1, :)
    p(:, :, :, 0) = p(:, :, :, 1)

    p(:, nx+1, :, :) = p(:, nx, :, :)
    p(:, :, ny+1, :) = p(:, :, ny, :)
    p(:, :, :, nz+1) = p(:, :, :, nz)

    ! compute the flux along the r direction
    ! do i = 0, nx
    !   call fluxes(p(:,i),p(:,i+1),fx(:,i), r(i))
    ! end do


  end do
  
end program nada