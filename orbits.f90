! File: Orbits.f90
! Author: Ian May

program orbits

  use utility, only: fp
  use timestep, only: take_step

  implicit none
  
  integer, parameter :: nSteps = 1000  ! Number of time steps to use
  !integer, parameter :: nSteps = 2000   ! increased number of steps(a system of my choosing)
  integer :: nT                        ! Loop variable for time updates
  !!! Step 1
  real (fp) :: tFinal
  real (fp) :: dt
  !integer, parameter :: x = 2, y = 2
  real (fp), dimension(2) :: mass
  real (fp), dimension(2,2, nSteps) :: pos
  real (fp), dimension(2,2,nSteps) :: mom
  tFinal = 50
  dt = tFinal/nSteps


  ! Give all particles initial mass, position, and momentum
  call set_ics()

  ! Fill rest of the array by integrating in time
  do nT=1,nSteps-1
    call take_step(dt,mass,pos(:,:,nT),mom(:,:,nT),pos(:,:,nT + 1),mom(:,:,nT + 1))
  end do

  ! Write the arrays to a file for plotting
  call write_data()

contains

  subroutine set_ics()
    implicit none
    !!! Step 2
    ! Set particle masses
    !mass = (/1.0_fp,0.01_fp/)
    !finding solution for other initial condition (m2=0.1)
    !mass = (/1.0_fp,0.1_fp/)
    mass = (/1.0_fp,1.0_fp/)!heavier mass2=1.0


    ! First particle position and momentum
    !set first particle to lie at origin
    pos(:,1,1) = (/0.0_fp,0.0_fp/)
    mom(:,1,1) = mass(1)*(/0.0_fp,0.0_fp/)
    
    print *, "momentum = "
    print *, mom(:,1,1)
    !new velocity for first particle that cancels out with second particle
    !mom(:,1,1) = mass(1)*(/-1.0_fp,0.0_fp/)


    
    ! Second particle position and momentum
    pos(:,2,1) = (/0.0_fp,-1.0_fp/)
    mom(:,2,1) = mass(2)*(/1.0_fp,0.0_fp/)
  print *, "momentum2 = "
    print *, mom(:,2,1)
    
  end subroutine set_ics

  subroutine write_data()
    implicit none
    open(20,file = "sol.dat",status = "replace")
    do nT=1,nSteps
      write(20,*) pos(:,:,nT),mom(:,:,nT)
    end do
    close(20)
  end subroutine write_data

end program orbits
