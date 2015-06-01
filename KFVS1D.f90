!-----------------------------------------------------------------------------
!KFVS1D - shock structure calculation using Kinetic Flux Vector Splitting Scheme
!-----------------------------------------------------------------------------
!Copyright (C) 2015 Ao Xu <http://www.linkedin.com/pub/ao-xu/30/a72/a29>
!-----------------------------------------------------------------------------

!--------------------------------------------------
!>store the global variables
!--------------------------------------------------
module global_data
    !--------------------------------------------------
    !kind selection
    !--------------------------------------------------
    integer,parameter :: RKD = 8

    !--------------------------------------------------
    !variables to control the simulation
    !--------------------------------------------------
    real(kind=RKD),parameter :: PI = 4.0*atan(1.0) !Pi
    real(kind=RKD),parameter :: SMV = tiny(real(1.0,8)) !small value to avoid 0/0
    real(kind=RKD),parameter :: UP = 1.0 !used in sign() function
    real(kind=RKD) :: cfl !global CFL number
    real(kind=RKD) :: dt !global time step
    real(kind=RKD) :: sim_time !current simulation time
    real(kind=RKD) :: max_time !maximum simulation time
    character(len=20),parameter :: HSTFILENAME = "shock.hst" !history file name
    character(len=20),parameter :: RSTFILENAME = "shock_KFVS.dat" !result file name
    integer :: iter !iteration

    !--------------------------------------------------
    !gas properties
    !--------------------------------------------------
    real(kind=RKD) :: gam !ratio of specific heat
    integer :: ck !internal degree of freedom

    !--------------------------------------------------
    !macros for a readable code
    !--------------------------------------------------
    !I/O
    integer,parameter :: HSTFILE = 20 !history file ID
    integer,parameter :: RSTFILE = 21 !result file ID
    !--------------------------------------------------
    !basic derived type
    !--------------------------------------------------
    !cell center
    type :: cell_center
        !geometry
        real(kind=RKD) :: x !cell center coordinates
        real(kind=RKD) :: length !length
        !flow field
        real(kind=RKD) :: w(3) !density, x-momentum,total energy
        real(kind=RKD) :: prim(3)
    end type cell_center

    !cell interface
    type :: cell_interface
        real(kind=RKD) :: flux(3) !mass flux, x momentum flux, energy flux
    end type cell_interface

    !--------------------------------------------------
    !flow field
    !--------------------------------------------------
    !index method
    !     ----------------
    !  (i)|      (i)     |(i+1)
    !     ----------------
    integer :: ixmin,ixmax !index range
    type(cell_center),allocatable,dimension(:) :: ctr !cell centers
    type(cell_interface),allocatable,dimension(:) :: vface !vertical interfaces

end module global_data

!--------------------------------------------------
!>define some commonly used functions/subroutines
!--------------------------------------------------
module tools
    use global_data
    implicit none
    contains
        !--------------------------------------------------
        !>convert primary variables to conservative variables
        !>@param[in] prim          :primary variables
        !>@return    get_conserved :conservative variables
        !--------------------------------------------------
        function get_conserved(prim)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD) :: get_conserved(3)

            get_conserved(1) = prim(1)
            get_conserved(2) = prim(1)*prim(2)
            get_conserved(3) = 0.5*prim(1)/prim(3)/(gam-1.0)+0.5*prim(1)*prim(2)**2
        end function get_conserved

        !--------------------------------------------------
        !>convert conservative variables to primary variables
        !>@param[in] w           :conservative variables
        !>@return    get_primary :conservative variables
        !--------------------------------------------------
        function get_primary(w)
            real(kind=RKD),intent(in) :: w(3)
            real(kind=RKD) :: get_primary(3) !primary variables

            get_primary(1) = w(1)
            get_primary(2) = w(2)/w(1)
            get_primary(3) = 0.5*w(1)/(gam-1.0)/(w(3)-0.5*w(2)**2/w(1))
        end function get_primary

        !--------------------------------------------------
        !>obtain ratio of specific heat
        !>@param[in] ck        :internal degree of freedom
        !>@return    get_gamma :ratio of specific heat
        !--------------------------------------------------
        function get_gamma(ck)
            integer,intent(in) :: ck
            real(kind=RKD) :: get_gamma

            get_gamma = float(ck+3)/float(ck+1)
        end function get_gamma

        !--------------------------------------------------
        !>obtain speed of sound
        !>@param[in] prim    :primary variables
        !>@return    get_sos :speed of sound
        !--------------------------------------------------
        function get_sos(prim)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD) :: get_sos !speed of sound

            get_sos = sqrt(0.5*gam/prim(3))
        end function get_sos

end module tools

!--------------------------------------------------
!>flux calculation
!--------------------------------------------------
module flux
    use global_data
    use tools
    implicit none

    integer,parameter :: MNUM = 6 !number of normal velocity moments
    integer,parameter :: MTUM = 4 !number of tangential velocity moments

    contains
        !--------------------------------------------------
        !>calculate flux of inner interface
        !>@param[in]    cell_L :cell left to the target interface
        !>@param[inout] face   :the target interface
        !>@param[in]    cell_R :cell right to the target interface
        !--------------------------------------------------
        subroutine calc_flux(cell_L,face,cell_R)
            type(cell_center),intent(in) :: cell_L,cell_R
            type(cell_interface),intent(inout) :: face
            real(kind=RKD) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM),Mxi(0:2) !<u^n>,<u^n>_{>0},<u^n>_{<0},<\xi^l>

!!!face%flux(1) = cell_L%prim(1)*( cell_L%prim(2)/2.0*erfc(-sqrt(cell_L%prim(3))*cell_L%prim(2)) &
!!!                                +0.5*exp(-cell_L%prim(3)*cell_L%prim(2)**2)/sqrt(PI*cell_L%prim(3)) ) &
!!!              +cell_R%prim(1)*( cell_R%prim(2)/2.0*erfc(sqrt(cell_R%prim(3))*cell_R%prim(2)) &
!!!                                -0.5*exp(-cell_R%prim(3)*cell_R%prim(2)**2)/sqrt(PI*cell_R%prim(3)) )
!!!
!!!face%flux(2) = cell_L%prim(1)*( (cell_L%prim(2)**2/2.0+1.0/4.0/cell_L%prim(3))*erfc(-sqrt(cell_L%prim(3))*cell_L%prim(2)) &
!!!                                +0.5*cell_L%prim(2)*exp(-cell_L%prim(3)*cell_L%prim(2)**2)/sqrt(PI*cell_L%prim(3)) ) &
!!!              +cell_R%prim(1)*( (cell_R%prim(2)**2/2.0+1.0/4.0/cell_R%prim(3))*erfc(sqrt(cell_R%prim(3))*cell_R%prim(2)) &
!!!                                -0.5*cell_R%prim(2)*exp(-cell_R%prim(3)*cell_R%prim(2)**2)/sqrt(PI*cell_R%prim(3)) )
!!!
!!!face%flux(3) = cell_L%prim(1)*( (cell_L%prim(2)**3/4.0+(ck+3)/8.0/cell_L%prim(3)*cell_L%prim(2)) &
!!!                                  *erfc(-sqrt(cell_L%prim(3))*cell_L%prim(2)) &
!!!                               +(cell_L%prim(2)**2/4.0+(ck+2)/8.0/cell_L%prim(3)) &
!!!                                  *exp(-cell_L%prim(3)*cell_L%prim(2)**2)/sqrt(PI*cell_L%prim(3)) ) &
!!!              +cell_R%prim(1)*( (cell_R%prim(2)**3/4.0+(ck+3)/8.0/cell_R%prim(3)*cell_R%prim(2)) &
!!!                                  *erfc(sqrt(cell_R%prim(3))*cell_R%prim(2)) &
!!!                                -(cell_R%prim(2)**2/4.0+(ck+2)/8.0/cell_R%prim(3)) &
!!!                                  *exp(-cell_R%prim(3)*cell_R%prim(2)**2)/sqrt(PI*cell_R%prim(3)) )


        call moment_u(cell_L%prim,Mu,Mxi,Mu_L,Mu_R)
        face%flux(1) = cell_L%prim(1)*Mu_L(1)
        face%flux(2) = cell_L%prim(1)*Mu_L(2)
        face%flux(3) = 0.5*cell_L%prim(1)*(Mu_L(3)+Mu_L(1)*Mxi(1))

        call moment_u(cell_R%prim,Mu,Mxi,Mu_L,Mu_R)
        face%flux(1) = face%flux(1)+cell_R%prim(1)*Mu_R(1)
        face%flux(2) = face%flux(2)+cell_R%prim(1)*Mu_R(2)
        face%flux(3) = face%flux(3)+0.5*cell_R%prim(1)*(Mu_R(3)+Mu_R(1)*Mxi(1))


        end subroutine calc_flux


        !--------------------------------------------------
        !>calculate moments of velocity
        !>@param[in] prim :primary variables
        !>@param[out] Mu        :<u^n>
        !>@param[out] Mxi       :<\xi^l>
        !>@param[out] Mu_L,Mu_R :<u^n>_{>0},<u^n>_{<0}
        !--------------------------------------------------
        subroutine moment_u(prim,Mu,Mxi,Mu_L,Mu_R)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD),intent(out) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM)
            real(kind=RKD),intent(out) :: Mxi(0:2)
            integer :: i

            !moments of normal velocity
            Mu_L(0) = 0.5*erfc(-sqrt(prim(3))*prim(2))
            Mu_L(1) = prim(2)*Mu_L(0)+0.5*exp(-prim(3)*prim(2)**2)/sqrt(PI*prim(3))
            Mu_R(0) = 0.5*erfc(sqrt(prim(3))*prim(2))
            Mu_R(1) = prim(2)*Mu_R(0)-0.5*exp(-prim(3)*prim(2)**2)/sqrt(PI*prim(3))

            do i=2,MNUM
                Mu_L(i) = prim(2)*Mu_L(i-1)+0.5*(i-1)*Mu_L(i-2)/prim(3)
                Mu_R(i) = prim(2)*Mu_R(i-1)+0.5*(i-1)*Mu_R(i-2)/prim(3)
            end do

            Mu = Mu_L+Mu_R

            !moments of \xi
            Mxi(0) = 1.0 !<\xi^0>
            Mxi(1) = 0.5*ck/prim(3) !<\xi^2>
            Mxi(2) = (ck**2+2.0*ck)/(4.0*prim(3)**2) !<\xi^4>
        end subroutine moment_u
end module flux

!--------------------------------------------------
!>UGKS solver
!--------------------------------------------------
module solver
    use global_data
    use tools
    use flux
    implicit none
    contains
        !--------------------------------------------------
        !>calculate time step
        !--------------------------------------------------
        subroutine timestep()
            real(kind=RKD) :: tmax !max 1/dt allowed
            real(kind=RKD) :: sos !speed of sound
            real(kind=RKD) :: prim(3) !primary variables
            integer :: i

            !set initial value
            tmax = 0.0

            do i=ixmin,ixmax
                !convert conservative variables to primary variables
                prim = get_primary(ctr(i)%w)

                !sound speed
                sos = get_sos(prim)

                !maximum velocity
                prim(2) = Abs(prim(2))+sos

                !maximum 1/dt allowed
                tmax = max(tmax,prim(2)/ctr(i)%length)
            end do

            !time step
            dt = cfl/tmax
        end subroutine timestep

        !--------------------------------------------------
        !>calculate the flux across the interfaces
        !--------------------------------------------------
        subroutine evolution()
            integer :: i

            do i=ixmin,ixmax+1 !with ghost cell
                call calc_flux(ctr(i-1),vface(i),ctr(i))
            end do

        end subroutine evolution

        !--------------------------------------------------
        !>update cell averaged values
        !--------------------------------------------------
        subroutine update()
            integer :: i

            do i=ixmin,ixmax
                !--------------------------------------------------
                !update W^{n+1}
                !--------------------------------------------------
                ctr(i)%w = ctr(i)%w+dt*(vface(i)%flux-vface(i+1)%flux)/ctr(i)%length !update W^{n+1}
                ctr(i)%prim = get_primary(ctr(i)%w)
            end do
        end subroutine update
end module solver

!--------------------------------------------------
!>input and output
!--------------------------------------------------
module io
    use global_data
    use tools
    implicit none
    contains
        !--------------------------------------------------
        !>main initialization subroutine
        !--------------------------------------------------
        subroutine init()
            real(kind=RKD) :: xlength !length of computational domain
            real(kind=RKD) :: xscale !cell size/mean free path

            !control
            cfl = 0.5 !CFL number
            max_time = 0.15d0 !output time interval

            !gas
            ck = 4 !internal degree of freedom
            gam = get_gamma(ck) !ratio of specific heat

            !geometry
            xlength = 1.0
            xscale = 0.001

            call init_geometry(xlength,xscale) !initialize the geometry
            call init_flow_field() !set the initial condition
        end subroutine init

        !--------------------------------------------------
        !>initialize the geometry
        !>@param[inout] xlength :domain length
        !>@param[in]    xscale  :cell size/mean free path
        !--------------------------------------------------
        subroutine init_geometry(xlength,xscale)
            real(kind=RKD),intent(inout) :: xlength
            real(kind=RKD),intent(in) :: xscale
            integer :: xnum !number of cells
            real(kind=RKD) :: dx !cell length
            integer :: i

            !adjust values
            xnum = nint(xlength/xscale)
            xlength = xnum*xscale
            write(*,*) "xnum=",xnum
            write(*,*) "xlength=",xlength

            !cell index range
            ixmin = 1
            ixmax = xnum

            !allocation
            allocate(ctr(ixmin-1:ixmax+1)) !cell center (with ghost cell)
            allocate(vface(ixmin:ixmax+1)) !vertical and horizontal cell interface

            !cell length and area
            dx = xlength/(ixmax-ixmin+1)
            write(*,*) "dx=",dx

            !cell center (with ghost cell)
            forall(i=ixmin-1:ixmax+1)
                ctr(i)%x = (i-0.5)*dx
                ctr(i)%length = dx
            end forall
        end subroutine init_geometry

        !--------------------------------------------------
        !>set the initial condition
        !>@param[in] Ma_L :Mach number in front of shock
        !--------------------------------------------------
        subroutine init_flow_field()
!!!            real(kind=RKD) :: Ma_L !Mach number in front of shock
!!!            real(kind=RKD) :: Ma_R !Mach number after shock
!!!            real(kind=RKD) :: ratio_T !T2/T1
            real(kind=RKD) :: prim_L(3), prim_R(3) !primary variables before and after shock
            real(kind=RKD) :: w_L(3), w_R(3) !conservative variables before and after shock
            integer :: i
!-------------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!! prim(1) = density
!!!!!!!!!!!!!! prim(2) = velocity
!!!!!!!!!!!!!! prim(3) = lambda = 1/temperature
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!! w(1) = density
!!!!!!!!!!!!!! w(2) = density*velocity
!!!!!!!!!!!!!! w(3) = density*energy=0.5*density*velocity**2+0.5*(ck+1)*pressure
!-------------------------------------------------------------------------------------------------------------
            !upstream condition (before shock)
            w_L(1) = 1.0d0
            w_L(2) = -2.0d0
            w_L(3) = 3.0d0

            !downstream condition (after shock)
            w_R(1) = 1.0d0
            w_R(2) = 2.0d0
            w_R(3) = 3.0d0

            prim_L = get_primary(w_L)
            prim_R = get_primary(w_R)
!-------------------------------------------------------------------------------------------------------
            !initialize field (with ghost cell)
            forall(i=ixmin-1:(ixmin+ixmax)/2)
                ctr(i)%prim = prim_L
                ctr(i)%w = w_L
            endforall

            forall(i=(ixmin+ixmax)/2+1:ixmax+1)
                ctr(i)%prim = prim_R
                ctr(i)%w = w_R
            endforall
        end subroutine init_flow_field

        !--------------------------------------------------
        !>write result
        !--------------------------------------------------
        subroutine output()
            integer :: i
            !--------------------------------------------------
            !write to file
            !--------------------------------------------------
            !open result file and write header
            open(unit=RSTFILE,file=RSTFILENAME,status="unknown",action="write")
            write(RSTFILE,*) "VARIABLES = Position, Density, Velocity, Pressure, Energy"
            do i=ixmin,ixmax
                write(RSTFILE,*) ctr(i)%x, ctr(i)%prim(1), ctr(i)%prim(2), ctr(i)%prim(1)/ctr(i)%prim(3)/2.0d0, &
                                 ctr(i)%prim(1)/ctr(i)%prim(3)/2.0d0/(gam-1)/ctr(i)%prim(1)
            enddo
            !close file
            close(RSTFILE)
        end subroutine output
end module io

!--------------------------------------------------
!>main program
!--------------------------------------------------
program main
    use global_data
    use solver
    use io
    implicit none

    !initialization
    call init()

    !set initial value
    iter = 1 !number of iteration
    sim_time = 0.0 !simulation time

    !open file and write header
    open(unit=HSTFILE,file=HSTFILENAME,status="replace",action="write") !open history file
    write(HSTFILE,*) "VARIABLES = iter, sim_time, dt" !write header

    !iteration
    do while(.true.)
        call timestep() !calculate time step
        if((sim_time+dt)>=max_time) dt=max_time-sim_time
        call evolution() !calculate flux across the interfaces
        call update() !update cell averaged value

        !write iteration situation every 10 iterations
        if (mod(iter,100)==0) then
            write(*,"(A18,I15,2E15.7)") "iter,sim_time,dt:",iter,sim_time,dt
            write(HSTFILE,"(I15,2E15.7)") iter,sim_time,dt
        end if

        !check if output
        if (sim_time>=max_time) exit

        iter = iter+1
        sim_time = sim_time+dt
    end do

    write(*,"(A18,I15,2E15.7)") "iter,sim_time,dt:",iter,sim_time,dt
    write(HSTFILE,"(I15,2E15.7)") iter,sim_time,dt

    !close history file
    close(HSTFILE)

    !output solution
    call output()
end program main

! vim: set ft=fortran tw=0:
