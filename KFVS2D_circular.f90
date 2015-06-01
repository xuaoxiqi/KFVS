!-----------------------------------------------------------------------------
!KFVS2D - solve water shallow equation using Kinetic Flux Vector Splitting Scheme
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
    real(kind=RKD) :: res(3) !residual
    real(kind=RKD) :: eps !convergence criteria
    real(kind=RKD) :: sim_time !current simulation time
    real(kind=RKD) :: max_time !maximum simulation time
    character(len=20),parameter :: HSTFILENAME = "cavity.hst" !history file name
    character(len=20),parameter :: RSTFILENAME = "cavity_KFVS.dat" !result file name
    integer :: iter !iteration
    integer :: method_output !output the solution with normalized value or not

    !--------------------------------------------------
    !gas properties
    !--------------------------------------------------
    real(kind=RKD) :: gam !ratio of specific heat
    real(kind=RKD) :: omega !temperature dependence index in HS/VHS/VSS model
    real(kind=RKD) :: pr !Prandtl number
    real(kind=RKD) :: mu_ref !viscosity coefficient in reference state
    integer :: ck !internal degree of freedom

    !--------------------------------------------------
    !macros for a readable code
    !--------------------------------------------------
    !direction
    integer,parameter :: IDIRC = 1 !i direction
    integer,parameter :: JDIRC = 2 !j direction
    !I/O
    integer,parameter :: HSTFILE = 20 !history file ID
    integer,parameter :: RSTFILE = 21 !result file ID

    !--------------------------------------------------
    !basic derived type
    !--------------------------------------------------
    !cell center
    type :: cell_center
        !geometry
        real(kind=RKD) :: x,y !cell center coordinates
        real(kind=RKD) :: length(2) !length
        !flow field
        real(kind=RKD) :: w(3) !height, height*xVelocity, height*yVelocity
        real(kind=RKD) :: prim(3)
    end type cell_center

    !cell interface
    type :: cell_interface
        !flow flux
        real(kind=RKD) :: flux(3) !mass flux, x and y momentum flux
    end type cell_interface

    !--------------------------------------------------
    !flow field
    !--------------------------------------------------
    !index method
    !     ----------------
    !  (i)|      (i)     |(i+1)
    !     ----------------
    integer :: ixmin,ixmax,iymin,iymax !index range
    type(cell_center),allocatable,dimension(:,:) :: ctr !cell centers
    type(cell_interface),allocatable,dimension(:,:) :: vface, hface !vertical and horizontal interfaces

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
            get_conserved(3) = prim(1)*prim(3)
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
            get_primary(3) = w(3)/w(1)
        end function get_primary

        !--------------------------------------------------
        !>obtain ratio of specific heat
        !>@param[in] ck        :internal degree of freedom
        !>@return    get_gamma :ratio of specific heat
        !--------------------------------------------------
        function get_gamma(ck)
            integer,intent(in) :: ck
            real(kind=RKD) :: get_gamma

            get_gamma = float(ck+4)/float(ck+2)
        end function get_gamma

end module tools

!--------------------------------------------------
!>flux calculation
!--------------------------------------------------
module flux
    use global_data
    use tools
    implicit none

    integer,parameter :: MNUM = 5 !number of normal velocity moments
    integer,parameter :: MTUM = 5 !number of tangential velocity moments

    contains
        !--------------------------------------------------
        !>calculate flux of inner interface
        !>@param[in]    cell_L :cell left to the target interface
        !>@param[inout] face   :the target interface
        !>@param[in]    cell_R :cell right to the target interface
        !--------------------------------------------------
        subroutine calc_flux(cell_L,face,cell_R,idx)
            type(cell_center),intent(in) :: cell_L,cell_R
            type(cell_interface),intent(inout) :: face
            integer,intent(in) :: idx
            real(kind=RKD) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM) !<u^n>,<u^n>_{>0},<u^n>_{<0}
            real(kind=RKD) :: Mv(0:MNUM),Mv_L(0:MNUM),Mv_R(0:MNUM) !<v^n>,<v^n>_{>0},<v^n>_{<0}

            if(idx.EQ.IDIRC) then

                call moment_u(cell_L%prim,Mu,Mu_L,Mu_R)
                call moment_v(cell_L%prim,Mv,Mv_L,Mv_R)
                face%flux(1) = cell_L%prim(1)*Mu_L(1)*Mv(0)
                face%flux(2) = cell_L%prim(1)*Mu_L(2)*Mv(0)
                face%flux(3) = cell_L%prim(1)*Mu_L(1)*Mv(1)

                call moment_u(cell_R%prim,Mu,Mu_L,Mu_R)
                call moment_v(cell_R%prim,Mv,Mv_L,Mv_R)
                face%flux(1) = face%flux(1)+cell_R%prim(1)*Mu_R(1)*Mv(0)
                face%flux(2) = face%flux(2)+cell_R%prim(1)*Mu_R(2)*Mv(0)
                face%flux(3) = face%flux(3)+cell_R%prim(1)*Mu_R(1)*Mv(1)

            elseif(idx.EQ.JDIRC) then

                call moment_u(cell_L%prim,Mu,Mu_L,Mu_R)
                call moment_v(cell_L%prim,Mv,Mv_L,Mv_R)
                face%flux(1) = cell_L%prim(1)*Mu(0)*Mv_L(1)
                face%flux(2) = cell_L%prim(1)*Mu(1)*Mv_L(1)
                face%flux(3) = cell_L%prim(1)*Mu(0)*Mv_L(2)

                call moment_u(cell_R%prim,Mu,Mu_L,Mu_R)
                call moment_v(cell_R%prim,Mv,Mv_L,Mv_R)
                face%flux(1) = face%flux(1)+cell_R%prim(1)*Mu(0)*Mv_R(1)
                face%flux(2) = face%flux(2)+cell_R%prim(1)*Mu(1)*Mv_R(1)
                face%flux(3) = face%flux(3)+cell_R%prim(1)*Mu(0)*Mv_R(2)

            endif

        end subroutine calc_flux


        !--------------------------------------------------
        !>calculate moments of velocity
        !>@param[in] prim :primary variables
        !--------------------------------------------------
        subroutine moment_u(prim,Mu,Mu_L,Mu_R)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD),intent(out) :: Mu(0:MNUM),Mu_L(0:MNUM),Mu_R(0:MNUM)
            integer :: i

            !moments of normal velocity
            Mu_L(0) = 0.5*erfc(-sqrt(1.0d0/prim(1))*prim(2))
            Mu_L(1) = prim(2)*Mu_L(0)+0.5*exp(-1.0d0/prim(1)*prim(2)**2)/sqrt(PI/prim(1))
            Mu_R(0) = 0.5*erfc(sqrt(1.0d0/prim(1))*prim(2))
            Mu_R(1) = prim(2)*Mu_R(0)-0.5*exp(-1.0d0/prim(1)*prim(2)**2)/sqrt(PI/prim(1))

            do i=2,MNUM
                Mu_L(i) = prim(2)*Mu_L(i-1)+0.5*(i-1)*Mu_L(i-2)*prim(1)
                Mu_R(i) = prim(2)*Mu_R(i-1)+0.5*(i-1)*Mu_R(i-2)*prim(1)
            end do

            Mu = Mu_L+Mu_R

        end subroutine moment_u

        subroutine moment_v(prim,Mv,Mv_L,Mv_R)
            real(kind=RKD),intent(in) :: prim(3)
            real(kind=RKD),intent(out) :: Mv(0:MNUM),Mv_L(0:MNUM),Mv_R(0:MNUM)
            integer :: i

            !moments of normal velocity
            Mv_L(0) = 0.5*erfc(-sqrt(1.0d0/prim(1))*prim(3))
            Mv_L(1) = prim(3)*Mv_L(0)+0.5*exp(-1.0d0/prim(1)*prim(3)**2)/sqrt(PI/prim(1))
            Mv_R(0) = 0.5*erfc(sqrt(1.0d0/prim(1))*prim(3))
            Mv_R(1) = prim(3)*Mv_R(0)-0.5*exp(-1.0d0/prim(1)*prim(3)**2)/sqrt(PI/prim(1))

            do i=2,MNUM
                Mv_L(i) = prim(3)*Mv_L(i-1)+0.5*(i-1)*Mv_L(i-2)*prim(1)
                Mv_R(i) = prim(3)*Mv_R(i-1)+0.5*(i-1)*Mv_R(i-2)*prim(1)
            end do

            Mv = Mv_L+Mv_R

        end subroutine moment_v
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
            real(kind=RKD) :: prim(3) !primary variables
            integer :: i,j

            !set initial value
            tmax = 0.0

            do j=iymin,iymax
                do i=ixmin,ixmax
                    !convert conservative variables to primary variables
                    prim = get_primary(ctr(i,j)%w)

                    !maximum velocity
                    prim(2) = Abs(prim(2))+Sqrt(prim(1))

                    !maximum 1/dt allowed
                    tmax = max(tmax,prim(2)/ctr(i,j)%length(1))
                end do
            enddo

            !time step
            dt = cfl/tmax
        end subroutine timestep


        !--------------------------------------------------
        !>calculate the flux across the interfaces
        !--------------------------------------------------
        subroutine evolution()
            integer :: i,j

            do j=iymin,iymax
                do i=ixmin,ixmax+1 !with ghost cell
                    call calc_flux(ctr(i-1,j),vface(i,j),ctr(i,j),IDIRC)
                end do
            enddo

            do j=iymin,iymax+1
                do i=ixmin,ixmax !with ghost cell
                    call calc_flux(ctr(i,j-1),hface(i,j),ctr(i,j),JDIRC)
                end do
            enddo

        end subroutine evolution

        !--------------------------------------------------
        !>update cell averaged values
        !--------------------------------------------------
        subroutine update()
            real(kind=RKD) :: sum_res(3),sum_avg(3)
            integer :: i,j

            !set initial value
            res = 0.0
            sum_res = 0.0
            sum_avg = 0.0

            do j=iymin,iymax
                do i=ixmin,ixmax
                    !--------------------------------------------------
                    !update W^{n+1}
                    !--------------------------------------------------
                    ctr(i,j)%w = ctr(i,j)%w+dt*(vface(i,j)%flux-vface(i+1,j)%flux)/ctr(i,j)%length(1) &
                                           +dt*(hface(i,j)%flux-hface(i,j+1)%flux)/ctr(i,j)%length(2)  !update W^{n+1}
                    ctr(i,j)%prim = get_primary(ctr(i,j)%w)
                enddo
            enddo
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
            real(kind=RKD) :: iner_gas(3), outer_gas(3) !initial condition
            integer :: xnum,ynum

            !control
            cfl = 0.01 !CFL number
            write(*,*) "CFL =",cfl
            max_time = 0.69d0 !output time interval
            write(*,*) "Max Time =",max_time

            !gas
            ck = 4 !internal degree of freedom
            gam = get_gamma(ck) !ratio of specific heat

            !geometry
            xnum = 500
            ynum = 500

            !initial condition (height,height*xVelocity,height*vVelocity)
            iner_gas = [1.0, 0.0, 0.0]
            outer_gas = [0.1, 0.0, 0.0]

            call init_geometry(xnum,ynum) !initialize the geometry
            call init_flow_field(iner_gas,outer_gas) !set the initial condition
        end subroutine init

        !--------------------------------------------------
        !>initialize the geometry
        !>@param[inout] xlength :domain length
        !>@param[in]    xscale  :cell size/mean free path
        !--------------------------------------------------
        subroutine init_geometry(xnum,ynum)
            integer :: xnum, ynum !number of cells
            real(kind=RKD) :: dx, dy !cell length
            integer :: i,j

            !adjust values
            write(*,*) "xnum=",xnum,",  ynum=",ynum

            !cell index range
            ixmin = 1
            ixmax = xnum
            iymin = 1
            iymax = ynum

            !allocation
            allocate(ctr(ixmin-1:ixmax+1,iymin-1:iymax+1)) !cell center (with ghost cell)
            allocate(vface(ixmin:ixmax+1,iymin:iymax),hface(ixmin:ixmax,iymin:iymax+1)) !vertical and horizontal cell interface

            !cell length and area
            dx = 5.0d0/(ixmax-ixmin+1)
            dy = 5.0d0/(iymax-iymin+1)
            write(*,*) "dx=",dx,",  dy=",dy

!            cell center (with ghost cell)
            forall(i=ixmin-1:ixmax+1,j=iymin-1:iymax+1)
                ctr(i,j)%x = (i-0.5)*dx
                ctr(i,j)%y = (j-0.5)*dy
                ctr(i,j)%length(1) = dx
                ctr(i,j)%length(2) = dy
            end forall


        end subroutine init_geometry

        !--------------------------------------------------
        !>set the initial condition
        !>@param[in] Ma_L :Mach number in front of shock
        !--------------------------------------------------
        subroutine init_flow_field(iner_gas,outer_gas)
            real(kind=RKD),intent(in) :: iner_gas(3),outer_gas(3)
            real(kind=RKD) :: w_iner(3),w_outer(3)
            integer :: i,j
            integer :: nxHalf,nyHalf
            integer :: radius

            nxHalf = (ixmax-ixmin)/2+ixmin
            nyHalf = (iymax-iymin)/2+iymin
            radius = 110

            !initialize field (with ghost cell)
            do j=iymin-1,iymax+1
                do i=ixmin-1,ixmax+1
                    if ( ((i-nxHalf)**2+(j-nyHalf)**2).LE.radius**2 ) then
                        ctr(i,j)%prim = iner_gas
                    else
                        ctr(i,j)%prim = outer_gas
                    endif
                enddo
            enddo


            w_iner = get_conserved(iner_gas)
            w_outer = get_conserved(outer_gas)
            do j=iymin-1,iymax+1
                do i=ixmin-1,ixmax+1
                    if ( ((i-nxHalf)**2+(j-nyHalf)**2).LE.radius**2 ) then
                        ctr(i,j)%w = w_iner
                    else
                        ctr(i,j)%w = w_outer
                    endif
                enddo
            enddo

        end subroutine init_flow_field

        !--------------------------------------------------
        !>write result
        !--------------------------------------------------

        subroutine output()
            integer :: i,j
            character(len=100) :: filename

            write(filename,*) iter
            filename = adjustl(filename)

            open(unit=01,file="results-"//trim(filename)//".plt",status="unknown")
            write(01,*) 'title="kfvs-2d"'
            write(01,*) 'VARIABLES="X" "Y" "H"'
            write(01,101) ixmax-ixmin+3,iymax-iymin+3

            do j=iymin-1,iymax+1
                do i=ixmin-1,ixmax+1
                    write(01,100) ctr(i,j)%x, ctr(i,j)%y, ctr(i,j)%w(1)
                enddo
            enddo

100 format(2x,10(e13.6,' '))
101 format('ZONE',1x,'I=',1x,i5,2x,'J=',1x,i5,1x,'F=POINT')
            close(01)
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

    call output()
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
            if (mod(iter,500)==0) then
                call output()
            endif
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
