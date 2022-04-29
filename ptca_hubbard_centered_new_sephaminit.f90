program ptca_repulsive


!!!!!!!!!!!! defining parameters for the problems!!!!!!!!!!!!!!!!
!!! in this code only m,theta and phi are transferred between the
!!! processes and the new hamitlonian is initialized using these 
!!! this should be able to save time compared to transferring hamiltonian also
  implicit none
  include "mpif.h"
  integer(8) :: i,j,ki
  integer :: my_id
  integer :: num_procs
  integer(8) :: site_clster,loc_proc
  real(8) :: tvar,rnum !! variable used to store intermediate temperature
  integer(8),parameter :: ns_unit = 3
  integer(8),parameter :: L = 10 !! system size
  integer(8),parameter :: n_sites = L * L !! number of sites in the lattice
  integer(8),parameter :: cls_sites =  8 !! cluster size
  integer(8),parameter :: ncl_by2 = 0.5*(cls_sites)+1 !! dividing cls_sites by 2
  integer(8),parameter :: n_splits = (ncl_by2)*(ncl_by2)
  integer(8),parameter :: split_sites = n_sites/n_splits
  integer(8),parameter :: cls_dim = (cls_sites)*(cls_sites) !! number of sites in the cluster
  integer(8),parameter :: n_equil  = 0 !! no of equilibrium steps
  integer(8),parameter :: n_meas  = 0 !! no of measurements
  integer(8),parameter :: meas_skip = 10 ! make measurement after this mc cycles
  integer(8),parameter :: dim_h = 6*n_sites  ! dimensionality of hamiltonian
  integer(8),parameter :: dim_clsh = 6*cls_dim ! dimensionality of cluster hamiltonian
  real(8),parameter :: temp = 0.30  !! simulation temperature
  real(8),parameter :: dtemp = 0.30 !! temperature step to lower the temperature
  real(8),parameter :: t_min = 0.01 !! minimum temperature for the simulation
  real(8),parameter :: pi = 4*atan(1.0)
  real(8),parameter :: t_hopping = 1.0
  real(8),parameter :: u_int = 2.0
  real(8),parameter :: mu = 0.5*(u_int)
  real(8),parameter :: m_max=2.0_8,m_min=0.0_8
  real :: t_strt_equil, t_end_equil
  real :: t_strt_meas , t_end_meas
  real :: delT 
  !!! this array will be initialized to -1 at the starting 
  !!! entry will be changed to 1 when that particular site is updated during mc
  integer(8),dimension(0:split_sites-1)::changed_ids,loc_ids

  !!!!!!!!!!!!!!!!!! initialize neighbour table !!!!!!!!!!!!!!!!!!!!!
  integer(8),dimension(0:n_sites-1)::sites
  integer(8),dimension(0:n_sites-1):: right,left,up,down,right_up,left_down

!!!!!!!!!!!!!!!!!!!!!!!!!!!! array that will hold all the sites in the cluster!!!
  integer(8), dimension(0:n_sites-1,0:cls_dim-1)::cl_st ! sites in the cluster at site j
  integer(8),dimension(0:n_splits-1,0:split_sites-1) :: sites_array !! array to store information about the split sites

  !!!!!!!!!!!!!!!!!! variational parameters of the monte carlo procedure !!!!!!!!
  real(8),dimension(0:3,0:n_sites-1):: m,m_loc   !! m
  real(8),dimension(0:3,0:n_sites-1):: theta,loc_theta !! theta
  real(8),dimension(0:3,0:n_sites-1):: phi,loc_phi  !! phi

!!!!!!!!!!!!!!! full hamiltonian and cluster hamiltonian !!!!!!!!!!!!
  complex(8),dimension(0:dim_h-1,0:dim_h-1) :: hamiltonian
  complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1) :: hamil_cls

!!!!!!!!!!!!!!!  parameters for the lapack !!!!!!!!!!!!!!!!!!!!!!!!!
integer(8),parameter :: lwork  = (2*dim_clsh)+(dim_clsh**2)
integer(8),parameter :: lrwork = 2*(dim_clsh**2)+(5*(dim_clsh)+1)
integer(8),parameter :: liwork = (5*dim_clsh)+3
integer(8) :: info = 10

complex(8),dimension(lwork)::work
real(8),dimension(lrwork) :: rwork
integer(8),dimension(liwork) :: iwork

!!!!!!!!!!!!!!!!!! files names !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=200):: fname


!!!!!!!!!!!!!!!!!! parameters for the parallelization !!!!!!!!!!!!!!!

integer :: ierr
integer, dimension(MPI_STATUS_SIZE)::status

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

    !! calling the subroutine to initialize the neighbours
    call neighbour_table(right,left,up,down,L,n_sites,sites,right_up,left_down)
   
    !! subroutinee to initialize the mc variables(m,theta,phi) at each site
    if (my_id == 0) then
        call mcvar_init(ns_unit,n_sites,m,theta,phi,pi,m_min,m_max)
    end if
    
    !! synchronize all the processes and broadcast m configuration
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(m,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    !! synchronize all the processes and broadcast theta configuration
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(theta,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    !! synchronize all the processes and broadcast theta configuration
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(phi,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


    !! subroutine to initialize the full lattice hamiltonian (2L*L,2*L*L) by master
    
    call ham_init(right,left,up,down,right_up,left_down,L,n_sites,hamiltonian,t_hopping,&
                            mu,u_int,m,theta,phi,dim_h,ns_unit)

    !! wait for all the process to finish this
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !! subroutine to initialize the arrays with sites splits in 4 subgroup 
    call cluster_sites(cls_sites,L,n_sites,cl_st,cls_dim)

    !! subroutine to split the lattice based on the cluster dimensions
    call  lattice_splt(split_sites,n_splits,n_sites,ncl_by2,sites_array,L)
  
      !! initialize the temperature and we will loop over this
    tvar  = temp
    !!! temperature loop over all the temperatures
    do while (tvar > t_min)
    print *,tvar,"started"
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !! time when the equilibration started
    call cpu_time(t_strt_equil)
    !!! Equlibration cycle
    do i = 0, n_equil, 1
       !! loop over all the splits 
!       print *,m,my_id 
        do j=0,n_splits-1,1
          !! intializing changed vars to -1 and broadcast it to all the processes
          if (my_id==0) then
            do loc_proc=0,split_sites-1,1
                changed_ids(loc_proc) = -1
              enddo
          end if
          
          !! all process should wait here till root send the data
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_BCAST(changed_ids,split_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          
          !! loop over all the sites within the partition
          do ki=my_id,split_sites-1,num_procs !uncomment this one to parallelize
            
            site_clster = sites_array(j,ki)
            changed_ids(ki) = site_clster
            !print *,'before sweep',m(site_clster),my_id,site_clster
            !!  initialize cluster hamiltonian
            call cluster_ham(site_clster,L,n_sites,cls_sites, &
                          hamil_cls,cls_dim,t_hopping,hamiltonian,dim_h,dim_clsh,cl_st,ns_unit)

            !!  try to update the mc variables at the given site
            call  mc_sweep(cls_sites,hamil_cls,dim_h,dim_clsh,n_sites,m,theta,phi,site_clster,&
                 cls_dim,hamiltonian,mu,u_int,pi,work,lwork,rwork,lrwork,iwork,liwork,info,tvar,cl_st,m_min,m_max)
            !print *,'after sweep',m(site_clster),my_id,site_clster

!          print *,'bb',m(site_clster),my_id,site_clster
          end do
          
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!          print *,'bb',m(j*(split_sites-1):(j+1)*(split_sites-1)),my_id
          
          
          !!! transfer the new m,theta,phi,hamiltonian between all the processors
          if (my_id==0) then
            !! loop over all the process and recieve the data into the root processes
            do loc_proc=1,num_procs-1,1
              call MPI_RECV(loc_ids,split_sites,MPI_DOUBLE_PRECISION,loc_proc,11,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(m_loc,n_sites,MPI_DOUBLE_PRECISION,loc_proc,12,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(loc_theta,n_sites,MPI_DOUBLE_PRECISION,loc_proc,13,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(loc_phi,n_sites,MPI_DOUBLE_PRECISION,loc_proc,14,MPI_COMM_WORLD,status,ierr)

              !! loop over the loc_ids array and get the site index that is updated
              do ki=0,split_sites-1,1
                if (loc_ids(ki)>=0) then
                  m(loc_ids(ki)) = m_loc(loc_ids(ki))
                  theta(loc_ids(ki)) = loc_theta(loc_ids(ki))
                  phi(loc_ids(ki)) = loc_phi(loc_ids(ki))
                endif
              enddo
            end do
          else
              loc_ids = changed_ids
              m_loc = m 
              loc_theta = theta
              loc_phi = phi
              call MPI_SEND(loc_ids,split_sites,MPI_DOUBLE_PRECISION,0,11,MPI_COMM_WORLD,ierr)
              call MPI_SEND(m_loc,n_sites,MPI_DOUBLE_PRECISION,0,12,MPI_COMM_WORLD,ierr)
              call MPI_SEND(loc_theta,n_sites,MPI_DOUBLE_PRECISION,0,13,MPI_COMM_WORLD,ierr)
              call MPI_SEND(loc_phi,n_sites,MPI_DOUBLE_PRECISION,0,14,MPI_COMM_WORLD,ierr)
          end if

          !! synchronize all the processes       
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !! broadcast the updated m vlaues from root
        call MPI_BCAST(m,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !! broadcast the updated theta vlaues from root
        call MPI_BCAST(theta,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        !! broadcast the updated phi vlaues from root
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(phi,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !! initializing the most updated hamiltonian using the updated
        !! monte carlo configurations of m,theta and phi
        call ham_init(right,left,up,down,right_up,left_down,L,n_sites,hamiltonian,t_hopping,&
                            mu,u_int,m,theta,phi,dim_h,ns_unit)
        
        !do loc_proc=0,10,1
        !print *,'ab',m(j*(split_sites-1):(j+1)*(split_sites-1)),my_id
        !print *,'ab',theta(j*(split_sites-1):(j+1)*(split_sites-1)),my_id
        !print *, phi(j*(split_sites-1):(j+1)*(split_sites-1)),my_id
        !end do

        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !print *,'ab',m,my_id
      
      end do
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !! time when the equilibration cycle finishes
      call cpu_time(t_end_equil)
      delT  = t_end_equil-t_strt_equil
      open(21,file='total_equilibration_time_L10_cl8',action='write',position='append')
      if (my_id==0) then
          write(21,*) my_id, tvar,delT 
          do i=1,num_procs-1,1
              call MPI_RECV(delT,1,MPI_DOUBLE_PRECISION,i,69,MPI_COMM_WORLD,status,ierr)
              write(21,*) i, tvar,delT 
          end do
      else 
              call MPI_SEND(delT,1,MPI_DOUBLE_PRECISION,0,69,MPI_COMM_WORLD,ierr)
      end if
      close(21)

      !!! measurement cycle
      call cpu_time(t_strt_meas)
      do i = 1, n_meas, 1
        !print *,'measurement loop with temp',tvar
        !! loop over all partition of the lattice
        do j=0,n_splits-1,1
  
          !! intializing changed vars to -1 and broadcast it to all the processes
          if (my_id==0) then
            do loc_proc=0,split_sites-1,1
                changed_ids(loc_proc) = -1
              enddo
          end if
          
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_BCAST(changed_ids,split_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !! loop over all the sites within the lattice
         do ki=my_id,split_sites-1,num_procs !uncomment this one to parallelize
            site_clster = sites_array(j,ki)
            changed_ids(ki) = site_clster

            !!    initialize cluster hamiltonian
            call cluster_ham(site_clster,L,n_sites,cls_sites, &
                                hamil_cls,cls_dim,t_hopping,hamiltonian,dim_h,dim_clsh,cl_st,ns_unit)

            !!     try to update the mc variables at the given site
            call  mc_sweep(cls_sites,hamil_cls,dim_h,dim_clsh,n_sites,m,theta,phi,site_clster,&
                    cls_dim,hamiltonian,mu,u_int,pi,work,lwork,rwork,lrwork,iwork,liwork,info,tvar,cl_st,m_min,m_max)
            !! loop over the sites in each non-interacting split of the lattice
            end do
        
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         
          !!! transfer the new m,theta,phi,hamiltonian between all the processors
          if (my_id==0) then
            !!! loop over all the processors and recieve the data from each one of them
            do loc_proc=1,num_procs-1,1
              call MPI_RECV(loc_ids,split_sites,MPI_DOUBLE_PRECISION,loc_proc,28,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(m_loc,n_sites,MPI_DOUBLE_PRECISION,loc_proc,38,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(loc_theta,n_sites,MPI_DOUBLE_PRECISION,loc_proc,48,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(loc_phi,n_sites,MPI_DOUBLE_PRECISION,loc_proc,58,MPI_COMM_WORLD,status,ierr)
              
              !!! set the variables arrays in master using values from other slave processes
              do ki=0,split_sites-1,1
              if (loc_ids(ki)>=0) then
                m(loc_ids(ki)) = m_loc(loc_ids(ki))
                theta(loc_ids(ki)) = loc_theta(loc_ids(ki))
                phi(loc_ids(ki)) = loc_phi(loc_ids(ki))               
               endif
              enddo
            end do

          !!! send the information about the site that is changed and the observables that are changed to the master
          else 
              loc_ids = changed_ids
              m_loc = m 
              loc_theta = theta
              loc_phi = phi
              call MPI_SEND(loc_ids,split_sites,MPI_DOUBLE_PRECISION,0,28,MPI_COMM_WORLD,ierr)
              call MPI_SEND(m_loc,n_sites,MPI_DOUBLE_PRECISION,0,38,MPI_COMM_WORLD,ierr)
              call MPI_SEND(loc_theta,n_sites,MPI_DOUBLE_PRECISION,0,48,MPI_COMM_WORLD,ierr)
              call MPI_SEND(loc_phi,n_sites,MPI_DOUBLE_PRECISION,0,58,MPI_COMM_WORLD,ierr)
            end if         
          
          !!! synchronize all the processes 
          !!! send the updated m configurations
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_BCAST(m,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          
          !!! send the updated theta configurations
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_BCAST(theta,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          
          !!! send the updated phi configurations
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_BCAST(phi,n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
          
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

          !! initializing the most updated hamiltonian using the updated
          !! monte carlo configurations of m,theta and phi
          call ham_init(right,left,up,down,L,n_sites,hamiltonian,t_hopping,&
                            mu,u_int,m,theta,phi,dim_h)
        
        !!! loop end for all the splits
        end do    
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !! printing the results in an output file after every meas_skip mc cycles
        if (my_id==0) then 
        if (mod(i,meas_skip)==0) then
          
          call print_f(fname,u_int,tvar,L,cls_sites)
          print *,i,tvar,fname

          open(16,file=fname,action='write',position='append')
          do j=0,n_sites-1,1
            write(16,20) i,j,m(j),theta(j),phi(j)
            20  format(I4,2X,I4,2X,ES22.8,2X,ES22.8,2X,ES22.8)
          end do
          close(16)
        end if
      end if 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !!! end of the measurement loop
      end do

      call cpu_time(t_end_meas)
      delT = t_end_meas-t_strt_meas 
      print *,'measurement time elapsed', delT,my_id
      open(22,file='total_measurement_time_L10_cl8',action='write',position='append')
      if (my_id==0) then
          write(22,*) my_id, tvar, delT 
          do i=1,num_procs-1,1
              call MPI_RECV(delT,1,MPI_DOUBLE_PRECISION,i,69,MPI_COMM_WORLD,status,ierr)
              write(22,*) i, tvar,delT 
          end do
      else 
              call MPI_SEND(delT,1,MPI_DOUBLE_PRECISION,0,69,MPI_COMM_WORLD,ierr)
      end if
      close(22)
      !print *,tvar,"finished"
      !!! lower the temperature of the system
      tvar = tvar-dtemp
    end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call MPI_FINALIZE(ierr)

end program ptca_repulsive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! setting up the neighbour table!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine neighbour_table(right,left,up,down,L,n_sites,sites,right_up,left_down)
implicit none
    integer(8) :: i,xi,yi,ri,li,ui,di
    integer(8) :: L
    integer(8) :: n_sites
    integer(8),dimension(0:n_sites-1):: right,left,up,down,right_up,left_down
    integer(8),dimension(0:n_sites-1)::sites

    !!! sites in the lattice starting from 0-->(n_sites-1)
    do i = 0, n_sites-1, 1
      sites(i) = i
    end do

    do i = 0, n_sites-1, 1
      yi  =  i/L
      xi  = mod(i,L)
      ri = mod(xi+1,L)
      li = mod(xi-1+L,L)
      ui = mod(yi+1,L)
      di = mod(yi-1+L,L)
      right(i) = ri + (yi * L)
      left(i) = li + (yi * L)
      up(i) = xi + (ui * L)
      down(i) = xi + (di * L)
      left_down(i) = li + (di*L)
      right_up(i) = ri + (ui*L)
      !print *,i,right(i),left(i),up(i),down(i)
    end do

end subroutine neighbour_table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! ---------------------------------------------------------------------------!!
!! -------------------- initialize monte carlo variables----------------------!!
!! ---------------------------------------------------------------------------!!
subroutine mcvar_init(ns_unit,n_sites,m,theta,phi,pi,m_min,m_max)
implicit none
  integer(8) :: n_sites,i,ns_unit,j
  
  real(8) :: pi
  real(8) :: m_min,m_max
  real(8):: rand3,rand1,rand2
  real(8) :: rand_int3,rand_int1,rand_int2
  real(8),dimension(0:n_sites-1) :: m,theta,phi  

  do i = 0, n_sites-1, 1
    do j=0,ns_unit-1,1
    call random_number(rand1)
    call random_number(rand2)
    call random_number(rand3)

    rand_int1=int(rand1*5000)
    m(j,i)=((((m_min)**3)+(((m_max)**3)-((m_min)**3)))*(rand_int1/5000.0))**(1.0_8/3.0_8)

    rand_int2=int(rand2*1000)
    theta(j,i)=acos(-1.0_8+(rand_int2/500.0))

    rand_int3=int(rand3*2000)
    phi(j,i)=(2.0_8*pi)*(rand_int3/2000.0)
    end do

    
  end do
end subroutine mcvar_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!initializing the hamiltonian!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ham_init(right,left,up,down,right_up,left_down,L,n_sites,hamiltonian,t_hopping,mu, &
                      u_int,m,theta,phi,dim_h,ns_unit)
implicit none
  integer :: i,id0,id1,id2
  integer :: s0l,s0ld 
  integer :: s1r,s1d
  integer :: s2u,s2ur
  integer(8) :: L,ns_unit
  integer(8) :: n_sites
  integer(8) :: ri,li,ui,di,rui,ldi
  integer(8) :: dim_h
  real(8) :: t_hopping,mu,u_int
  integer(8),dimension(0:n_sites-1) :: right,left,up,down,right_up,left_down

  real(8) ::mx,my,mz
  real(8),dimension(0:ns_unit,0:n_sites-1)::m,theta,phi
  complex(8),dimension(0:dim_h-1,0:dim_h-1) :: hamiltonian
  hamiltonian(:,:) = cmplx(0.0,0.0)

  !! setting hopping dependent terms
  call haminitinUnitcell(hamiltonian,ns_unit,n_sites,dim_h)

  do i = 0,n_sites-1, 1
    !! sites in the unit cell
    id0 = (ns_unit*i); id1 = (ns_unit*i)+1; id2 = (ns_unit*i)+2

    !! neighbouring unit cell 
    ri = right(i); li = left(i); ui = up(i); di = down(i); rui = right_up(i); ldi = left_down(i)

    !! neighbour of site 0
    s0l = (ns_unit*li) + 1 ; s0ld = (ns_unit*ldi)+2

    !! neighbour of site 1
    s1r = (ns_unit*ri) ; s1d = (ns_unit*di + 2)

    !! neighbour of site 2
    s2u = (ns_unit*ui) + 1; s2ur = (ns_unit*rui)


    !!! mx for three sites in the unit cell
    mx0 = m(0,i)  * cos(phi(0,i)) * sin(theta(0,i))
    mx1 = m(1,i)  * cos(phi(1,i)) * sin(theta(1,i))
    mx2 = m(2,i)  * cos(phi(2,i)) * sin(theta(2,i))

    !!! my for three sites in the unit cell
    my0 = m(0,i) * sin(phi(0,i)) * sin(theta(0,i))
    my1 = m(1,i) * sin(phi(1,i)) * sin(theta(1,i))
    my2 = m(2,i) * sin(phi(1,i)) * sin(theta(2,i))
    
    !! mz for three sites in the unit cell
    mz0 = m(0,i) *  cos(theta(0,i));mz1 = m(1,i) *  cos(theta(1,i));mz2 = m(2,i) *  cos(theta(2,i))
    
    !! for up spins
    hamiltonian(id0,id0) =  -(mu-0.5*u_int) - (0.5*u_int)*mz0
    hamiltonian(id1,id1) =  -(mu-0.5*u_int) - (0.5*u_int)*mz1
    hamiltonian(id2,id2) =  -(mu-0.5*u_int) - (0.5*u_int)*mz2

    !! for down spins
    hamiltonian(id0+n_sites,id0+n_sites) =  -(mu-0.5*u_int) + (0.5*u_int)*mz0
    hamiltonian(id1+n_sites,id1+n_sites) =  -(mu-0.5*u_int) + (0.5*u_int)*mz1
    hamiltonian(id2+n_sites,id2+n_sites) =  -(mu-0.5*u_int) + (0.5*u_int)*mz2

    ! setting the updn and dnup components for site 0 in the unit cell
    hamiltonian(id0,id0+n_sites) = -(0.5*u_int)*cmplx(mx0,-my0)
    hamiltonian(id0+n_sites,id0) = -(0.5*u_int)*cmplx(mx0,my0)


    ! setting the updn and dnup components for site 1 in the unit cell
    hamiltonian(id1,id1+n_sites) = -(0.5*u_int)*cmplx(mx1,-my1)
    hamiltonian(id1+n_sites,id1) = -(0.5*u_int)*cmplx(mx1,my1)

    ! setting the updn and dnup components for site 2 in the unit cell
    hamiltonian(id2,id2+n_sites) = -(0.5*u_int)*cmplx(mx2,-my2)
    hamiltonian(id2+n_sites,id2) = -(0.5*u_int)*cmplx(mx2,my2)

  end do
  !print *,hamiltonian
end subroutine ham_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! this function will set the hopping element. It will set the hopping between the sites
!!! within the unit cell

subroutine haminitinUnitcell(hamiltonian,ns_unit,n_sites,dim_h,right,left,up,down,t_hopping)
implicit none  
  integer :: id0 ,id1,id2 ,id0n , id1n , id2n
  integer(8) :: ns_unit , n_sites,dim_h
  integer(8),dimension(0:n_sites-1)::right,left,up,down
  real(8) :: t_hopping 
  complex(8),dimension(0:dim_h-1,0:dim_h-1) :: hamiltonian
  
  
  do i = 0,n_sites-1, 1
    id0 = (ns_unit*i) !! first site in the unit cell i
    id1 = (ns_unit*i)+ 1 !! second site in the unit cell i
    id2 = (ns_unit*i)+ 2  !! third site in the unit cell i
    id0n = id0+n_sites  !! index for the spin down in unit cell i position 0
    id1n = id1+n_sites  !! index for the spin down in unit cell i position 1
    id2n = id2+n_sites  !! index for the spin down in unit cell i position 2
    
    !!! hopping from site 0-->1 & 0--->2 for up and down spin
    hamiltonian(id0,id1) = -t_hopping; hamiltonian(id0n,id1n) = -t_hopping
    hamiltonian(id0,id2) = -t_hopping; hamiltonian(id0n,id2n) = -t_hopping

    ! conjugate hopping 1--->0 & 1--->2 for up and down spin
    hamiltonian(id1,id0) = -t_hopping; hamiltonian(id1n,id0n) = -t_hopping
    hamiltonian(id1,id2) = -t_hopping; hamiltonian(id1n,id2n) = -t_hopping

    ! conjugate hopping from 2-->0 and 2--->1 for up and down spin
    hamiltonian(id2,id0) = -t_hopping; hamiltonian(id2n,id0n) = -t_hopping
    hamiltonian(id2,id1) = -t_hopping; hamiltonian(id2n,id1n) = -t_hopping
  end do
call haminitOutunitcell(right,left,up,down,hamiltonian,ns_unit,n_sites,t_hopping,site,dim_h)
end subroutine haminitinUnitcell
!!-----------------------------------------------------------------------------!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! initializing the hamiltonian for the hopping for sites within the outsite
subroutine haminitOutunitcell(right,left,up,down,hamiltonian,ns_unit,n_sites,t_hopping,dim_h)
  integer(8) :: ns_unit,n_sites,site,dim_h
  integer :: si 
  integer :: si0, s0l, s0ld, si0n, s0ln, s0ldn ! sites related to 0 site in the unit cell 
  integer :: si1, s1r, s1d, si1n, s1rn, s1dn ! sites related to site 1 in the unit cell
  integer :: si2 , s2u, s2ur , si2n , s2un , s2urn ! sites related to site 2 in the unit cell
  real(8) :: t_hopping
  integer(8),dimension(0:n_sites-1) :: right,left,up,down
  complex(8),dimension(0:dim_h-1,0:dim_h-1) :: hamiltonian
  

  do si=0,n_sites-1,1  
    !! 1st site in the unit cell and site it's connected to 2nd and 3rd site of other unit cell
    si0 = ns_unit*si 
    s0l = left(si) * ns_unit + 1 ; s0ld = left(down(si))*ns_unit + 2
    si0n = si0+n_sites ; s0ln = s0l+n_sites ; s0ldn = s0ld+n_sites

    !! 2nd site in the unit cell is connected to 1st and 3rd site in the other unit cell
    si1 = (ns_unit*si)+1 
    s1r = right(si) * ns_unit ; s1d = ns_unit*down(si) + 2
    si1n = si1+n_sites; s1rn = s1r+n_sites;s1dn = s1d+n_sites

    !! 3rd site in the unit cell is connected to 1st site of unit cell in up and 0 site 
    !! of unit cell in up right direction
    si2 = (ns_unit*si)+2
    s2u = (up(site)*ns_unit)+1 ; s2ur = right(up(i))*ns_unit 
    si2n = si2+n_sites; s2un  = s2u + n_sites ; s2urn = s2ur+n_sites

    !! matrix element between 0 site and neighbour for up and down spins
    hamiltonian(si0,s0l) = -t_hopping ; hamiltonian(si0n,s0ln) = -t_hopping
    hamiltonian(s0l,si0) = -t_hopping ; hamiltonian(s0ln,si0n) = -t_hopping

    hamiltonian(si0,s0ld) = -t_hopping ; hamiltonian(si0n,s0ldn) = -t_hopping
    hamiltonian(s0ld,si0) = -t_hopping ; hamiltonian(s0ldn,si0n) = -t_hopping
  
    !! matrix element between 1 site and neighbour for up and down spins
    hamiltonian(si1,s1r) = -t_hopping;  hamiltonian(si1n,s1rn) = -t_hopping
    hamiltonian(s1r,si1) = -t_hopping; hamiltonian(s1rn,si1n) = -t_hopping

    hamiltonian(si1,s1d) = -t_hopping ; hamiltonian(si1n,s1dn) = -t_hopping
    hamiltonian(s1d,si1) = -t_hopping ; hamiltonian(s1dn,si1n) = -t_hopping

    !! matrix element between 2nd site and neighbour for up and down spins
    hamiltonian(si2,s2u)  = -t_hopping ; hamiltonian(si2n,s2un) = -t_hopping
    hamiltonian(s2u,si2)  = -t_hopping ; hamiltonian(s2un,si2n) = -t_hopping

    hamiltonian(si2,s2ur) = -t_hopping ; hamiltonian(si2n,s2urn) = -t_hopping
    hamiltonian(s2ur,si2) = -t_hopping ; hamiltonian(s2urn,si2n) = -t_hopping
  
  end do

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!initializing the cluster hamiltonian!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_ham(site_clster,L,n_sites,cls_sites, hamil_cls,cls_dim,&
                        t_hopping,hamiltonian,dim_h,dim_clsh,cl_st,ns_unit)
implicit none
  integer(8) :: cls_dim,L, cls_sites,n_sites,ns_unit
  integer(8) :: i,si,x,y,xip,yip
  integer :: sri,sui ! for the nn hopping
  integer(8) :: site_clster ! site that has the variable that will be changed
  integer(8) :: dim_h,dim_clsh
  real(8) :: t_hopping ! hopping strength

  !! arrays with all sites and the sites in the cluster
  integer(8),dimension(0:n_sites-1,0:cls_dim-1) :: cl_st

  complex(8),dimension(0:dim_h-1,0:dim_h-1) :: hamiltonian ! hamiltonian
  complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1)::hamil_cls ! cluster hamiltonian

  hamil_cls(:,:)=cmplx(0.0,0.0)
  
  !! constructing the cluster hamiltonian
  do i = 0, cls_dim-1, 1
      x = mod(i,cls_sites) ! x index in the  site cluster x --> [0,1,cls_sites-1]
      y = i/cls_sites       ! y index in the  site cluster y --> [0,1,cls_sites-1]
      si = x+(cls_sites*y)  ! site produced in the  site cluster si-->[0,cls_dim-1]


      xip = mod(x+1,cls_sites)  ! x+1 in the cluster with pbc (with cls_sites)
      yip = mod(y+1,cls_sites) ! y+1 in the cluster with pbc

      xim = mod(x-1+cls_sites,cls_sites) !! x-1 in the cluster (pbc)
      yim = mod(y-1+cls_sites,cls_sites) !! y-1 in the cluster (pbc)

      sri = xip + (y*cls_sites) ! right site of the cluster
      sui = x + (yip*cls_sites) ! up side of the cluster

      sli = xim + (y*cls_sites)  !! left site in the cluster
      sdi = x + (yim*cls_sites) !! down site in the cluster
      
      sild = xim + (yim*cls_sites) !! site in the left down of si
      siru = xip + (yip*cls_sites) !! site in the right up of si

      !print *,si,cl_st(site_clster,si),site_clster
      
      !!! diagonal part of the cluster hamiltonian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! setting up the up,down and down up terms for the 1st site in the unit cell.
      hamil_cls(si*ns_unit,(si*ns_unit)+cls_dim) = hamiltonian(cl_st(site_clster,si)*ns_unit,cl_st(site_clster,si)*ns_unit+n_sites)
      hamil_cls((si*ns_unit)+cls_dim,si*ns_unit) = hamiltonian(cl_st(site_clster,si)*ns_unit+n_sites,cl_st(site_clster,si)*ns_unit)

      !!! setting up the up,down and down up terms for the 2nd site in the unit cell
      hamil_cls(si*ns_unit+1,(si*ns_unit)+1+cls_dim) = hamiltonian(cl_st(site_clster,si)*ns_unit+1,cl_st(site_clster,si)*ns_unit+1+n_sites)
      hamil_cls((si*ns_unit)+cls_dim+1,si*ns_unit+1) = hamiltonian(cl_st(site_clster,si)*ns_unit+1+n_sites,cl_st(site_clster,si)*ns_unit+1)

      !!! setting up the up,down and down up terms for the 3rd site in the unit cell
      hamil_cls(si*ns_unit+2,(si*ns_unit)+cls_dim+2) = hamiltonian(cl_st(site_clster,si)*ns_unit+2,cl_st(site_clster,si)*ns_unit+n_sites+2)
      hamil_cls((si*ns_unit)+cls_dim+2,si*ns_unit+2) = hamiltonian(cl_st(site_clster,si)*ns_unit+n_sites+2,cl_st(site_clster,si)*ns_unit+2)


      !!! diagonal part of the hamiltonian for upup and dndn for site 0,1,2
      hamil_cls(ns_unit*si,ns_unit*si)=hamiltonian(cl_st(site_clster,si),cl_st(site_clster,si))
      hamil_cls(ns_unit*si+cls_dim,ns_unit*si+cls_dim)=hamiltonian(cl_st(site_clster,si)+n_sites,cl_st(site_clster,si)+n_sites)

      hamil_cls(ns_unit*si+1,ns_unit*si+1)=hamiltonian(cl_st(site_clster,si),cl_st(site_clster,si))
      hamil_cls(ns_unit*si+cls_dim,si+cls_dim)=hamiltonian(cl_st(site_clster,si)+n_sites,cl_st(site_clster,si)+n_sites)

      hamil_cls(ns_unit*si+2,ns_unit*si+2)=hamiltonian(cl_st(site_clster,si),cl_st(site_clster,si))
      hamil_cls(si+cls_dim,si+cls_dim)=hamiltonian(cl_st(site_clster,si)+n_sites,cl_st(site_clster,si)+n_sites)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      hamclsinitOutunitcell(right,left,up,down,hamil_cls,site,ns_unit,n_sites,t_hopping,site,dim_clsh)
      !print *,si,si+cls_dim!,cl_st(site_clster,si),cl_st(site_clster,si)+n_sites
    end do
    

!! uncomment to see the site and neighbour
!    print *,'ri',right_cls
!    print *,'ui',up_cls

end subroutine cluster_ham
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hamclsinitinUnitcell()

      !! 0--->1 & 0--->2 for spin up in the unit cell at si
      hamil_cls(ns_unit*si,(ns_unit*si)+1) = -t_hopping
      hamil_cls(ns_unit*si,(ns_unit*si)+2) = -t_hopping
      
      !! 1--->0 & 1---> 2 for spin up in the unit cell at si
      hamil_cls((ns_unit*si)+1,(ns_unit*si)) = -t_hopping
      hamil_cls((ns_unit*si)+1,(ns_unit*si)+2) = -t_hopping

      !! 2--->0 & 2---> 1 for spin up in the unit cell at si
      hamil_cls((ns_unit*si)+2,(ns_unit*si)) = -t_hopping
      hamil_cls((ns_unit*si)+2,(ns_unit*si)+1) = -t_hopping
      

      !! 0--->1 & 0--->2 for spin down in the unit cell at si
      hamil_cls(ns_unit*si+cls_dim,(ns_unit*si)+1+cls_dim) = -t_hopping
      hamil_cls(ns_unit*si+cls_dim,(ns_unit*si)+2+cls_dim) = -t_hopping
      
      !! 1--->0 & 1---> 2 for spin down in the unit cell at si
      hamil_cls((ns_unit*si)+1+cls_dim,(ns_unit*si)+cls_dim) = -t_hopping
      hamil_cls((ns_unit*si)+1+cls_dim,(ns_unit*si)+2+cls_dim) = -t_hopping

      !! 2--->0 & 2---> 1 for spin down in the unit cell at si
      hamil_cls((ns_unit*si)+2+cls_dim,(ns_unit*si)+cls_dim) = -t_hopping
      hamil_cls((ns_unit*si)+2+cls_dim,(ns_unit*si)+1+cls_dim) = -t_hopping

  end subroutine hamclsinitinUnitcell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! this subroutine will initialize the hamiltonian for the bonds outside the unit cell
subroutine hamclsinitOutunitcell(right,left,up,down,hamil_cls,site,ns_unit,n_sites,t_hopping,site,dim_clsh)
  integer(8) :: ns_unit,n_sites,site,dim_clsh
  real(8) :: t_hopping
  integer(8),dimension(0:ns_unit-1,0:n_sites-1) :: right,left,up,down
  complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1)
  !! 1st site in the unit cell and site it's connected to 2nd and 3rd site of other unit cell
  si0 = ns_unit*site 
  s0l = left(si) * ns_unit + 1 ; s0ld = left(down(si))*ns_unit + 2

  !! 2nd site in the unit cell is connected to 1st and 3rd site in the other unit cell
  si1 = (ns_unit*site)+1 
  s1r = right(si) * ns_unit ; s1d = ns_unit*down(si) + 2

  !! 3rd site in the unit cell is connected to 1st site of unit cell in up and 0 site 
  !! of unit cell in up right direction
  si2 = (ns_unit*site)+2
  s2u = (up(site)*ns_unit)+1 ; s2ur = right(up(i))*ns_unit 

  !! matrix element between 0 site and neighbour
  hamil_cls(si0,s0l) = -t_hopping ; hamil_cls(si0+n_sites,s0l+n_sites) = -t_hopping
  hamil_cls(s0l,si0) = -t_hopping ; hamil_cls(s0l+n_sites,si0+n_sites) = -t_hopping

  hamil_cls(si0,s0ld) = -t_hopping ; hamil_cls(si0+n_sites,s0ld+n_sites) = -t_hopping
  hamil_cls(s0ld,si0) = -t_hopping ; hamil_cls(s0ld+n_sites,si0+n_sites) = -t_hopping
  
  !! matrix element between 1 site and neighbour
  hamil_cls(si1,s1r) = -t_hopping;  hamil_cls(si1+n_sites,s1r+n_sites) = -t_hopping
  hamil_cls(s1r,si1) = -t_hopping;hamil_cls(s1r+n_sites,si1+n_sites) = -t_hopping

  hamil_cls(si1,s1d) = -t_hopping ; hamil_cls(si1+n_sites,s1d) = -t_hopping
  hamil_cls(s1d,si1) = -t_hopping ; hamil_cls(s1d,si1+n_sites) = -t_hopping

  !! matrix element between 2nd site and neighbour
  hamil_cls(si2,s2u)  = -t_hopping ; hamil_cls(si2+n_sites,s2u+n_sites) = -t_hopping
  hamil_cls(s2u+n_sites,si2)  = -t_hopping ; hamil_cls(s2u+n_sites,si2+n_sites) = -t_hopping

  hamil_cls(si2,s2ur) = -t_hopping ; hamil_cls(si2+n_sites,s2ur+n_sites) = -t_hopping
  hamil_cls(s2ur,si2) = -t_hopping ; hamil_cls(s2ur+n_sites,si2+n_sites) = -t_hopping


end subroutine hamclsinitOutunitcell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! generate cluster sites!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_sites(cls_sites,L,n_sites,cl_st,cls_dim)
implicit none
    integer(8) :: i,j,sil
    integer(8) :: xi,yi,xn,yn,n_sites,cls_dim
    integer(8) :: L,cls_sites
    integer(8) :: x_init,y_init
    integer(8),dimension(0:n_sites-1,0:cls_dim-1) :: cl_st

    !! generating all the sites in the cluster
    do j=0,n_sites-1,1
      i = 0

      !! x and y of the sites
      x_init = int(mod(j,L))
      y_init = int(j/L)

    !! y index of the site -cls/2>y>cls/2
    !! x index of the site -cls/2>x>cls/2

    do yi = -int(0.5*cls_sites),int(0.5*cls_sites)-1,1
      do xi = -int(0.5*cls_sites),int(cls_sites*0.5)-1, 1
        !! finding the x neighbour of the site
        if (xi>=0) then
            xn = mod(xi + x_init,L)
        else if(xi<0) then
            xn = mod(x_init-abs(xi)+L,L)
          end if

        !! finding the y neighbour of the site
          if (yi>=0) then
              yn = mod(yi + y_init,L)
          else if(yi<0) then
              yn = mod(y_init-abs(yi)+L,L)
            end if

        !! generating the site
        sil = xn+(yn*L)
        !! storing in the array
        cl_st(j,i) = sil
        !! incrementing the array index
        
        i=i+1
      end do
   end do
  !print *,cl_st(j,:)
  end do
 
end subroutine cluster_sites
!! ---------------------------------------------------------------------------!!


!!-----------------------------------------------------------------------------!!
!!---------------------------splitting lattice to------------------------------!!
!!-----------------------------------------------------------------------------!!
subroutine lattice_splt(split_sites,n_splits,n_sites,ncl_by2,sites_array,L)
implicit none
  !! given the x coordinate(i) and y coordinate(j) 
  !! this subroutine can categorize the lattice into 
  !! certain number of  groups based on the size of the cluster
 
  integer(8) :: n_sites !! number of sites in the lattice 
  integer(8) :: xi,yi  !! loop variables
  integer(8) :: i ,j
  integer(8) :: splt_var ! --> 0 , # of partitions of lattice
  integer(8) :: ki !! used in the array as an index
  integer(8) :: n_splits !! number of splits we want to create based on the cluster dim
  integer(8) :: split_sites
  integer(8) :: ncl_by2
  integer (8) :: l_minus_1 
  integer (8) :: L
  integer(8),dimension(0:n_splits-1,0:split_sites-1) :: sites_array !! array to store the sites and splits
  l_minus_1 = L - 1
  !print *,'SS',split_sites,'ns:',n_splits
  !! looping over the allowed values of x and y in the split lattice
  do j = 0 ,ncl_by2-1, 1
    do i = 0,ncl_by2-1,1
      !print *,"i:",i,"j:",j
      !! array index for the multidimensional array
      !! this determines the row index and the column is the site
      splt_var = ((ncl_by2)*j)+i 
      ki = 0
   !   print *,splt_var
      do yi = 0,l_minus_1,1
        do xi = 0,l_minus_1,1         
          !print *,i,j,splt_var,ncl_by2
          !print *,xi,yi,ki,splt_var
          if((mod(xi,(ncl_by2)) == i).and.(mod(yi,(ncl_by2)) == j)) then
           sites_array(splt_var,ki)=xi+(L*yi)
            ki = ki+1
          end if
         
        end do
      end do
    end do
  end do
  
end subroutine lattice_splt
!!-----------------------------------------------------------------------------!!


!! ---------------------------------------------------------------------------!!
!! -------------------- monte carlo sweep--------------------------------------!!
!! ---------------------------------------------------------------------------!!
subroutine mc_sweep(cls_sites,hamil_cls,dim_h,dim_clsh,n_sites,m,theta,phi,site_clster,&
     cls_dim,hamiltonian,mu,u_int,pi,work,lwork,rwork,lrwork,iwork,liwork,info,tvar,cl_st,m_min,m_max)
!implicit none
  integer(8) :: dim_h,dim_clsh,n_sites,site_clster,cls_dim
  integer(8) :: loc_site !! local site that is being flipped
  integer(8) :: info
  integer(8) :: lwork,lrwork,liwork
  integer(8) :: cls_sites
  real(8) :: mu, u_int
  real(8) :: pi,beta,e_u,e_v,enr_loc
  real(8)  :: tvar
  real(8) :: m_min,m_max,mc_prob
  integer(8), dimension(0:n_sites-1,0:cls_dim-1)::cl_st ! sites in the cluster at site j
 
  !!! variables for the monte carlo procedure
  real(8) :: tempmx0,tempmy0,tempmz0 ! to store initial value of mx,my,mz
  real(8) :: tempmx1,tempmy1,tempmz1 ! to store initial value of mx,my,mz
  real(8) :: tempmx2,tempmy2,tempmz2 ! to store initial value of mx,my,mz

  real(8) :: rand1, rand2 , rand3 ! random number to generate,m,theta,phi
  real(8) :: rand_int1,rand_int2, rand_int3
  real(8) :: delE,tempm0,temptheta0,tempphi0 
  real(8) :: tempm1,temptheta1,tempphi1
  real(8) :: tempm2,temptheta2,tempphi2

  complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1) :: hamil_cls
  complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1) :: temp_clsham
  complex(8),dimension(0:dim_h-1,0:dim_h-1) :: hamiltonian
  real(8),dimension(0:ns_unit,0:n_sites-1) :: m,theta,phi
  real(8),dimension(0:ns_unit,0:n_sites-1) :: loc_m
  real(8),dimension(0:ns_unit,0:dim_clsh-1) :: egval
  
  complex(8),dimension(lwork)::work
  real(8),dimension(lrwork) :: rwork
  integer(8),dimension(liwork)::iwork

  loc_m(:) = 0
  beta = 1./tvar
  ! position of global site in the cluster (center)
  loc_site = int((0.5*cls_sites)+cls_sites*(0.5*cls_sites))
!  print *,'loc',loc_site,site_clster
  !! generating uniform random numbers (new mc variables)
  !print *,'lcs',i,cls_sites,loc_site,site_clster
  
    
  call random_number(rand1)
  call random_number(rand2)
  call random_number(rand3)
  !! temp m
  rand_int1=rand1*5000
  temp0=((((m_min)**3)+(((m_max)**3)-((m_min)**3)))*(rand_int1/5000.0))**(1.0_8/3.0_8)
  
  !! temp theta
  rand_int2=rand2*1000
  temptheta0=acos(-1.0_8+(rand_int2/500.0))
  
  !! temp phi
  rand_int3=rand3*2000
  tempphi0=(2.0_8*pi)*(rand_int3/2000.0)
  

  !! value of mx,my,mz
  tempmx0 = tempm0 * cos(tempphi0) * sin(temptheta0)
  tempmy0 = tempm0 * sin(tempphi0) * sin(temptheta0)
  tempmz0 = tempm0 * cos(temptheta0)
    
  !!! monte carlo variable for the second site in the unit cell
  call random_number(rand1)
  call random_number(rand2)
  call random_number(rand3)
  !! temp m
  rand_int1=rand1*5000
  tempm1=((((m_min)**3)+(((m_max)**3)-((m_min)**3)))*(rand_int1/5000.0))**(1.0_8/3.0_8)
  
  !! temp theta
  rand_int2=rand2*1000
  temptheta1=acos(-1.0_8+(rand_int2/500.0))
  
  !! temp phi
  rand_int3=rand3*2000
  tempphi1=(2.0_8*pi)*(rand_int3/2000.0)
  

  !! value of mx,my,mz
  tempmx1 = tempm1 * cos(tempphi1) * sin(temptheta1)
  tempmy1 = tempm1 * sin(tempphi1) * sin(temptheta1)
  tempmz1 = tempm1 * cos(temptheta1)
    
  !!! monte carlo variables for the third site in the unit cell
  call random_number(rand1)
  call random_number(rand2)
  call random_number(rand3)
  !! temp m
  rand_int1=rand1*5000
  tempm2=((((m_min)**3)+(((m_max)**3)-((m_min)**3)))*(rand_int1/5000.0))**(1.0_8/3.0_8)
  
  !! temp theta
  rand_int2=rand2*1000
  temptheta2=acos(-1.0_8+(rand_int2/500.0))
  
  !! temp phi
  rand_int3=rand3*2000
  tempphi2=(2.0_8*pi)*(rand_int3/2000.0)
  

  !! value of mx,my,mz
  tempmx2 = tempm2 * cos(tempphi2) * sin(temptheta2)
  tempmy2 = tempm2 * sin(tempphi2) * sin(temptheta2)
  tempmz2 = tempm2 * cos(temptheta2)


!  print *,site_clster,tempmx,tempmy,tempmz,loc_site

  info = 10
  !print *,'tca_site',loc_site,'lattice_site',site_clster
  !! diagonalize old cluster hamiltonian but before that copy the copy of the hamiltonian
  !! and m configurations

  temp_clsham(:,:)  =  hamil_cls(:,:)
  loc_m(:) = m(:)
  
  egval(:)=0
  call zheevd('V','U', dim_clsh, temp_clsham, dim_clsh, egval, work, lwork, &
                                          rwork, lrwork, iwork, liwork, info)
  call enr_calc(egval,dim_clsh,loc_m,n_sites,cl_st,cls_dim,enr_loc,site_clster,tvar,u_int)
  e_u = enr_loc

  !! updating the diagonal part of the matrix, upup & dndn for the 1st site
  hamil_cls(loc_site*ns_unit,loc_site*ns_unit) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz0
  hamil_cls((loc_site*ns_unit)+cls_dim,(loc_site*ns_unit)+cls_dim) = -(mu-0.5*u_int) + (0.5*u_int)*tempmz0

  !! updating the blocks updn & dnup for the 1st site
  hamil_cls(loc_site*ns_unit,(loc_site*ns_unit)+cls_dim) = -(0.5*u_int)*cmplx(tempmx0,-tempmy0)
  hamil_cls((loc_site*ns_unit)+cls_dim,(ns_unit*loc_site)) = -(0.5*u_int)*cmplx(tempmx0,tempmy0)


  !! updating the diagonal part of the matrix, upup & dndn for the 2nd site
  hamil_cls(loc_site*ns_unit+1,loc_site*ns_unit+1) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz1
  hamil_cls((loc_site*ns_unit)+cls_dim+1,(loc_site*ns_unit)+1+cls_dim) = -(mu-0.5*u_int) + (0.5*u_int)*tempmz1

  !! updating the blocks updn & dnup for the 1st site
  hamil_cls(loc_site*ns_unit+1,(loc_site*ns_unit)+cls_dim+1) = -(0.5*u_int)*cmplx(tempmx1,-tempmy1)
  hamil_cls((loc_site*ns_unit)+1+cls_dim,(ns_unit*loc_site)+1) = -(0.5*u_int)*cmplx(tempmx1,tempmy1)


  !! updating the diagonal part of the matrix, upup & dndn for the 3rd site
  hamil_cls(loc_site*ns_unit+2,loc_site*ns_unit+2) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz2
  hamil_cls((loc_site*ns_unit)+cls_dim+2,(loc_site*ns_unit)+2+cls_dim) = -(mu-0.5*u_int) + (0.5*u_int)*tempmz2

  !! updating the blocks updn & dnup for the 1st site
  hamil_cls(loc_site*ns_unit+2,(loc_site*ns_unit)+cls_dim+2) = -(0.5*u_int)*cmplx(tempmx2,-tempmy2)
  hamil_cls((loc_site*ns_unit)+2+cls_dim,(ns_unit*loc_site)+2) = -(0.5*u_int)*cmplx(tempmx2,tempmy2)


  info = 10
  loc_m(site_clster*ns_unit) = tempm0
  loc_m((site_clster*ns_unit)+1) = tempm1
  loc_m((site_clster*ns_unit)+2) = tempm2 
  
  !! diagonalize the updated cluster hamiltonian
  
  egval(:)=0
  call zheevd('V','U', dim_clsh, hamil_cls, dim_clsh, egval, work, lwork, &
                                           rwork, lrwork, iwork, liwork, info)
  

  call enr_calc(egval,dim_clsh,loc_m,n_sites,cl_st,cls_dim,enr_loc,site_clster,tvar,u_int)
  e_v = enr_loc
  
  delE = e_v - e_u
  !print *,site_clster
  if ( delE <0.0 ) then
    !! setting up the mc variables
      m(site_clster*ns_unit) = tempm0 !! update the mc variable if energy is reduced  for the 1st site in the unit cell
      theta(site_clster*ns_unit) = temptheta0
      phi(site_clster*ns_unit) = tempphi0

      m(site_clster*ns_unit+1) = tempm1 !! update the mc variable if energy is reduced  for the 2nd site in the unit cell
      theta(site_clster*ns_unit+1) = temptheta1
      phi(site_clster*ns_unit+1) = tempphi1

      m(site_clster*ns_unit+2) = tempm2 !! update the mc variable if energy is reduced  for the 3rd site in the unit cell
      theta(site_clster*ns_unit+2) = temptheta2
      phi(site_clster*ns_unit+2) = tempphi2

      !!!! interaction term for the 1st site in the unit cell
      hamiltonian(site_clster*ns_unit,site_clster*ns_unit) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz0
      hamiltonian((site_clster*ns_unit)+n_sites,(site_clster*ns_unit)+n_sites) = -(mu - 0.5*u_int) +  (0.5*u_int)*tempmz0

      hamiltonian(site_clster*ns_unit,(site_clster*ns_unit)+n_sites) =   -(0.5*u_int)*cmplx(tempmx0,-tempmy0)
      hamiltonian((site_clster*ns_unit)+n_sites,site_clster*ns_unit) =  -(0.5*u_int)*cmplx(tempmx0,tempmy0)
  


      !!!! interaction term for the 2nd site in the unit cell
      hamiltonian((site_clster*ns_unit)+1,(site_clster*ns_unit)+1) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz1
      hamiltonian((site_clster*ns_unit)+n_sites,(site_clster*ns_unit)+n_sites) = -(mu - 0.5*u_int) +  (0.5*u_int)*tempmz1

      hamiltonian((site_clster*ns_unit)+1,(site_clster*ns_unit)+n_sites+1) =   -(0.5*u_int)*cmplx(tempmx1,-tempmy1)
      hamiltonian((site_clster*ns_unit)+1+n_sites,(site_clster*ns_unit)+1) =  -(0.5*u_int)*cmplx(tempmx1,tempmy1)
  

      !!!! interaction term for the 3nr site in the unit cell
      hamiltonian((site_clster*ns_unit)+2,(site_clster*ns_unit)+2) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz2
      hamiltonian((site_clster*ns_unit)+2+n_sites,(site_clster*ns_unit)+2+n_sites) = -(mu - 0.5*u_int) +  (0.5*u_int)*tempmz2

      hamiltonian((site_clster*ns_unit)+2,(site_clster*ns_unit)+n_sites+2) =   -(0.5*u_int)*cmplx(tempmx2,-tempmy2)
      hamiltonian((site_clster*ns_unit)+2+n_sites,(site_clster*ns_unit)+2) =  -(0.5*u_int)*cmplx(tempmx2,tempmy2)
  
  
  else
     call random_number(mc_prob)
     if(mc_prob < exp(-beta*delE)) then
        
      m(site_clster*ns_unit) = tempm0 !! update the mc variable if energy is reduced  for the 1st site in the unit cell
      theta(site_clster*ns_unit) = temptheta0
      phi(site_clster*ns_unit) = tempphi0

      m(site_clster*ns_unit+1) = tempm1 !! update the mc variable if energy is reduced  for the 2nd site in the unit cell
      theta(site_clster*ns_unit+1) = temptheta1
      phi(site_clster*ns_unit+1) = tempphi1

      m(site_clster*ns_unit+2) = tempm2 !! update the mc variable if energy is reduced  for the 3rd site in the unit cell
      theta(site_clster*ns_unit+2) = temptheta2
      phi(site_clster*ns_unit+2) = tempphi2

      !!!! interaction term for the 1st site in the unit cell
      hamiltonian(site_clster*ns_unit,site_clster*ns_unit) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz0
      hamiltonian((site_clster*ns_unit)+n_sites,(site_clster*ns_unit)+n_sites) = -(mu - 0.5*u_int) +  (0.5*u_int)*tempmz0

      hamiltonian(site_clster*ns_unit,(site_clster*ns_unit)+n_sites) =   -(0.5*u_int)*cmplx(tempmx0,-tempmy0)
      hamiltonian((site_clster*ns_unit)+n_sites,site_clster*ns_unit) =  -(0.5*u_int)*cmplx(tempmx0,tempmy0)
  
      !!!! interaction term for the 2nd site in the unit cell
      hamiltonian((site_clster*ns_unit)+1,(site_clster*ns_unit)+1) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz1
      hamiltonian((site_clster*ns_unit)+n_sites,(site_clster*ns_unit)+n_sites) = -(mu - 0.5*u_int) +  (0.5*u_int)*tempmz1

      hamiltonian((site_clster*ns_unit)+1,(site_clster*ns_unit)+n_sites+1) =   -(0.5*u_int)*cmplx(tempmx1,-tempmy1)
      hamiltonian((site_clster*ns_unit)+1+n_sites,(site_clster*ns_unit)+1) =  -(0.5*u_int)*cmplx(tempmx1,tempmy1)
  
      !!!! interaction term for the 3nr site in the unit cell
      hamiltonian((site_clster*ns_unit)+2,(site_clster*ns_unit)+2) = -(mu-0.5*u_int) - (0.5*u_int)*tempmz2
      hamiltonian((site_clster*ns_unit)+2+n_sites,(site_clster*ns_unit)+2+n_sites) = -(mu - 0.5*u_int) +  (0.5*u_int)*tempmz2

      hamiltonian((site_clster*ns_unit)+2,(site_clster*ns_unit)+n_sites+2) =   -(0.5*u_int)*cmplx(tempmx2,-tempmy2)
      hamiltonian((site_clster*ns_unit)+2+n_sites,(site_clster*ns_unit)+2) =  -(0.5*u_int)*cmplx(tempmx2,tempmy2)
      end if
  end if
end subroutine mc_sweep
!! ---------------------------------------------------------------------------!!



!! ---------------------------------------------------------------------------!!
!! ------------------------measurement of energy------------------------------!!
!! ---------------------------------------------------------------------------!!
subroutine enr_calc(egval,dim_clsh,loc_m,n_sites,cl_st,cls_dim,enr_loc,site_clster,tvar,u_int)
implicit none
  integer(8) :: i,n_sites,
  integer(8) :: si,ns_unit
  integer(8) :: dim_clsh,cls_dim,site_clster
  integer(8),dimension(0:n_sites-1,0:cls_dim-1)::cl_st
  real(8),dimension(0:dim_clsh-1):: egval
  real(8),dimension(0:ns_unit,0:n_sites-1)::loc_m
  real(8)::m_i,u_int
  real(8):: fac_be ! product of beta and e_i
  real(8) :: sum_e,sum_cl,tvar,beta,enr_loc,val


  beta = 1./tvar
  sum_e = 0.0
  sum_cl = 0.0
  
  !! calculating sum over all the log(1+beta E)
  do i = 0,dim_clsh-1, 1
    fac_be = -(beta*egval(i))
    if ( fac_be>40.0 ) then
        val = fac_be
    else 
        val = log(1.0+exp(fac_be))
    end if
    sum_e=sum_e+val
   end do
   sum_e = (-1.0)*(sum_e)/beta

  !! classical energies (U/4)sum m_{i}*m_{i}
  do i=0,cls_dim-1,1
    si = cl_st(site_clster,i)
    m_i = loc_m(0,si) + loc_m(1,si)+loc_m(2,si)
    sum_cl = sum_cl + (m_i)*(m_i)
  end do
    enr_loc = sum_e + (sum_cl)*(0.25*u_int)
end subroutine enr_calc

!! ---------------------------------------------------------------------------!!



!!-----------------------------------------------------------------------------!!
!!-------------------- printing -----------------------------------------------!!
!!-----------------------------------------------------------------------------!!
 subroutine print_f(fname,u_int,tvar,L,cls_sites)
   integer(8)::L,cls_sites
   real(8):: u_int,tvar
   character(len=50) :: format_L
   character(len=50) :: format_U
   character(len=50) :: format_temp="(F8.6)"
   character(len=50) :: format_cls="(I1)"
   character(len=20):: str_1
   character(len=20):: str_2
   character(len=20):: str_3
   character(len=20)::str_4   
   character(len=200):: fname



   !!! set format based on L
   if(L<10)then
     format_L="(I1)"
   else if(L>=10 .and. L<100) then
     format_L="(I2)"
   else if(L>100) then
     format_L="(I3)"
   end if

   
   !!! set format based on U
   if (u_int < 10.0) then
     format_U="(F5.3)"
   else if (u_int >= 10.0) then
     format_U="(F6.3)"
   end if

   

   write(str_1,format_L)L
   write(str_2,format_temp)tvar
   write(str_3,format_U)u_int
   write(str_4,format_cls)cls_sites

   fname=trim('sepconfigurations_L')//trim(str_1)//trim('_temp')//trim(str_2)//trim('_Uint')&
                     //trim(str_3)//trim('_cluster')//trim(str_4)
!   print*,fname

 end subroutine
