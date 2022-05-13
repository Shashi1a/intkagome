!!!!!!!!!!!! defining parameters for the problems!!!!!!!!!!!!!!!!
!!! in this code only m,theta and phi are transferred between the
!!! processes and the new hamitlonian is initialized using these 
!!! this should be able to save time compared to transferring hamiltonian also

program ptca_repulsive
use varmodule 
  implicit none
  include "mpif.h"
  integer(8) :: i,j,ki,loc_si1,loc_si2,loc_si3,suc0,suc1,suc2
  integer :: my_id,num_procs !! process id and number of processes
  integer(8) :: site_clster,loc_proc !! local site in the cluster
  real(8) :: rnum                 !! variable used to store intermediate temperature
  real(8) :: mu_init,sum_mu=0.0,mu_avg=0.0 !! mu calculations
  real(8) :: uninitial !! variable to store initial value of the interaction strength
  !!!! for time calcualtions
  real :: t_strt_equil, t_end_equil
  real :: t_strt_meas , t_end_meas
  real :: delT , mu_initial 

  complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1) :: copy_ham !! copy of cluster hamiltonian

!!!!!!!!!!!!!!!!!! parameters for the parallelization !!!!!!!!!!!!!!!
  integer :: ierr
  integer, dimension(MPI_STATUS_SIZE)::status

   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

    !! calling the subroutine to initialize the neighbours
    call neighbour_table()

    !! subroutinee to initialize the mc variables(m,theta,phi) at each site
    if (my_id == 0) then
        call mcvar_init()
    end if

    !! synchronize all the processes and broadcast m configuration
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(m,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    !! synchronize all the processes and broadcast theta configuration
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(theta,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    !! synchronize all the processes and broadcast theta configuration
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(phi,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(charge_confs,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    !! subroutine to initialize the full lattice hamiltonian (2L*L,2*L*L) by all processes
    !! since we have broadcasted monte carlo variables to all processes all of them will 
    !! have the same hamiltonian.
    
    !! initialize the hamiltonian and its copy to  zero
    hamiltonian(:,:) = cmplx(0.0,0.0)
    ham_noint = cmplx(0.0,0.0)

    !! set the non interacting part of the hamiltonian
    call  haminitinUnitcell()
    
    !! copy the non interacting part of the hamiltonian
    ham_noint(:,:) = hamiltonian(:,:)

    !!! set the interacting part of the hamiltonian
    call ham_init()

    !! initialize all the cluster hamiltonian matrix elements to zero 
    hamil_cls(:,:) = cmplx(0.0,0.0) 
    hamcls_noint(:,:) = cmplx(0.0,0.0)

    !! get the non interacting cluster hamiltonian
    call hamclsinitinUnitcell()
    
    !! copy the non interacting cluster hamiltonian
    hamcls_noint(:,:) = hamil_cls(:,:)

    !call zheevd('V','U', dim_clsh, hamil_cls, dim_clsh, egval, work, lwork, &
    !                                       rwork, lrwork, iwork, liwork, info)
               
    !print *,egval(0),egval(dim_clsh-1)
 
    !call zheevd('V','U', dim_h, hamiltonian, dim_h, egval_fl, work_full, lwork_full, &
    !                                       rwork_full, lrwork_full, iwork_full, liwork_full, info)
    !print *,egval_fl(0),egval_fl(dim_h-1),egval_fl(int(0.5*dim_h)),egval_fl(int(0.5*dim_h)-1)

    !hamiltonian(:,:) = cmplx(0.0,0.0)
    !hamiltonian(:,:) = ham_noint(:,:)
    !mu = (egval_fl(int(0.5*dim_h))+egval_fl(int(0.5*dim_h)-1))*0.5
    !mu = 0.0
    !!! set the interacting and diagonal part of the hamiltonian
    !call ham_init()

    !call zheevd('V','U', dim_h, hamiltonian, dim_h, egval_fl, work_full, lwork_full, &
    !                                       rwork_full, lrwork_full, iwork_full, liwork_full, info)

    !print *,'mus',egval_fl(0),egval_fl(dim_h-1),egval_fl(int(0.5*dim_h)),egval_fl(int(0.5*dim_h)-1)


    !open(89,file='eigenval_L16_nomufix.dat',action='write',position='append')            
    !do i=0,dim_h-1,1
    !  write(89,*) egval_fl(i)
    !end do
    !close(89)                
    

    !! wait for all the process to finish this
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !! subroutine to initialize the arrays with sites splits in 4 subgroup 
    call cluster_sites()

    !do j=0,15,1
    !  print * ,my_id, j ,cl_st(8,j)
    !end do
    !! subroutine to split the lattice based on the cluster dimensions
    call lattice_splt()
  
    !! initialize the temperature and we will loop over this
    tvar  = temp
   
    !!! temperature loop over all the temperatures
    do while (tvar > t_min)
      print *,tvar,"started"
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !! time when the equilibration started
      call cpu_time(t_strt_equil)
      call print_f()

      sum_mu = 0.0
      !!! Equlibration cycle
      do i = 0, n_equil, 1
      
        !!! for first 20 steps calculate the mu and use if for rest of the iterations
        !!! form a cluster centered around site 0
        if (i<mu_cnf) then
            open(81,file=fname_mu,action='write',position='append')
            !! form a cluster centered around site 0
            !site_clster = 0
            !mu = 0.0 !! set mu to zero for this calculation
            !hamil_cls(:,:) = cmplx(0.0,0.0) !! set all matrix element to zero
            !hamil_cls(:,:) = hamcls_noint(:,:) !! set the non interacting part of the hamiltonian
            !call cluster_ham(site_clster)  !! set the diagonal elements
            !info  = 10
            !call zheevd('V','U', dim_clsh, hamil_cls, dim_clsh, egval, work, lwork, &
            !                               rwork, lrwork, iwork, liwork, info)
            
            hamiltonian(:,:) = cmplx(0.0,0.0)
            hamiltonian(:,:) = ham_noint(:,:)
            call ham_init()
            info = 10
            call zheevd('V','U', dim_h, hamiltonian, dim_h, egval_fl, work_full, lwork_full, &
                                           rwork_full, lrwork_full, iwork_full, liwork_full, info)
            
            mu_init = 0.5*(egval_fl(int(0.5*dim_h)-1)+egval_fl(int(0.5*dim_h)))
            mu = mu_init
            sum_mu  = sum_mu + mu_init
            
            !print *,my_id,i,mu,egval_fl(dim_h-1),egval_fl(0),egval_fl(int(0.5*dim_h)-1)
            write(81,*) tvar,i,mu_init
            close(81)
        else  
            mu_avg = sum_mu/mu_cnf
            mu = mu_avg
            !write(81,*) tvar,i,mu 
            !print *,i,my_id,mu,mu_avg,int(0.5*dim_clsh)
        end if
        
        !!! loop over all the splits 
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

           
          !! loop over all the sites within each partition
          do ki=my_id,split_sites-1,num_procs !uncomment this one to parallelize
            site_clster = sites_array(j,ki)
            changed_ids(ki) = site_clster
            
            !  print *,'before sweep',my_id,site_clster
            !! initialize cluster hamiltonian 
            hamil_cls(:,:) = cmplx(0.0,0.0) !! set all elements to zero
            hamil_cls(:,:) = hamcls_noint(:,:) !! copy the non interacting hamiltonian
            
            !! initialize the diagonal part
            call cluster_ham(site_clster)
            
            !!  try to update the mc variables at the given site
            call mc_sweep(site_clster)
          end do
          
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          
          !!! transfer the new m,theta,phi,hamiltonian between all the processors
          if (my_id==0) then
            !! loop over all the process and recieve the data into the root processes
            do loc_proc=1,num_procs-1,1
              call MPI_RECV(loc_ids,split_sites,MPI_DOUBLE_PRECISION,loc_proc,11,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(m_loc,ns_unit*n_sites,MPI_DOUBLE_PRECISION,loc_proc,12,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(loc_theta,ns_unit*n_sites,MPI_DOUBLE_PRECISION,loc_proc,13,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(loc_phi,ns_unit*n_sites,MPI_DOUBLE_PRECISION,loc_proc,14,MPI_COMM_WORLD,status,ierr)

              !! loop over the loc_ids array and get the site index that is updated
              do ki=0,split_sites-1,1
                if (loc_ids(ki)>=0) then
                  loc_si1 = 0+(loc_ids(ki)*ns_unit)
                  loc_si2 = loc_si1+1
                  loc_si3 = 2+loc_si1
                  m(loc_si1) = m_loc(loc_si1)
                  m(loc_si2) = m_loc(loc_si2)
                  m(loc_si3) = m_loc(loc_si3)
                  theta(loc_si1) = loc_theta(loc_si1)
                  theta(loc_si2) = loc_theta(loc_si2)
                  theta(loc_si3) = loc_theta(loc_si3)
                  phi(loc_si1) = loc_phi(loc_si1)
                  phi(loc_si2) = loc_phi(loc_si2)
                  phi(loc_si3) = loc_phi(loc_si3)
                endif
              enddo
            end do
          else
              loc_ids = changed_ids
              m_loc = m 
              loc_theta = theta
              loc_phi = phi
              call MPI_SEND(loc_ids,split_sites,MPI_DOUBLE_PRECISION,0,11,MPI_COMM_WORLD,ierr)
              call MPI_SEND(m_loc,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,12,MPI_COMM_WORLD,ierr)
              call MPI_SEND(loc_theta,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,13,MPI_COMM_WORLD,ierr)
              call MPI_SEND(loc_phi,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,14,MPI_COMM_WORLD,ierr)
          end if

          !! synchronize all the processes       
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !! broadcast the updated m vlaues from root and wait for all the proccess
        call MPI_BCAST(m,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        !! broadcast the updated theta vlaues from root and wait for all the proccess
        call MPI_BCAST(theta,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        !! broadcast the updated phi vlaues from root
        call MPI_BCAST(phi,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !! initializing the most updated hamiltonian using the updated
        !! monte carlo configurations of m,theta and phi
        hamiltonian(:,:) = cmplx(0.0,0.0) !! initialize the hamiltonian to zero
        hamiltonian(:,:) = ham_noint(:,:) !! copy the non interacting part
        call ham_init()  !! intializing the interacting part

      end do
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
      
      end do
      !close(81)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      !! time when the equilibration cycle finishes
      call cpu_time(t_end_equil)
      delT  = t_end_equil-t_strt_equil
      open(21,file=fname_eqt,action='write',position='append')
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
      do i = 0,n_meas-1, 1
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
            hamil_cls(:,:) = cmplx(0.0,0.0)
            hamil_cls(:,:) = hamcls_noint(:,:)
            call cluster_ham(site_clster)
             
            !!     try to update the mc variables at the given site
            call  mc_sweep(site_clster)

            !! loop over the sites in each non-interacting split of the lattice
            end do
        
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         
          !!! transfer the new m,theta,phi,hamiltonian between all the processors
          if (my_id==0) then
            !!! loop over all the processors and recieve the data from each one of them
            do loc_proc=1,num_procs-1,1
              call MPI_RECV(loc_ids,split_sites,MPI_DOUBLE_PRECISION,loc_proc,28,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(m_loc,ns_unit*n_sites,MPI_DOUBLE_PRECISION,loc_proc,38,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(loc_theta,ns_unit*n_sites,MPI_DOUBLE_PRECISION,loc_proc,48,MPI_COMM_WORLD,status,ierr)
              call MPI_RECV(loc_phi,ns_unit*n_sites,MPI_DOUBLE_PRECISION,loc_proc,58,MPI_COMM_WORLD,status,ierr)
              
              !!! set the variables arrays in master using values from other slave processes
              do ki=0,split_sites-1,1
              if (loc_ids(ki)>=0) then
                loc_si1 = 0+(loc_ids(ki)*ns_unit)
                  loc_si2 = loc_si1+1
                  loc_si3 = 2+loc_si1
                  m(loc_si1) = m_loc(loc_si1)
                  m(loc_si2) = m_loc(loc_si2)
                  m(loc_si3) = m_loc(loc_si3)
                  theta(loc_si1) = loc_theta(loc_si1)
                  theta(loc_si2) = loc_theta(loc_si2)
                  theta(loc_si3) = loc_theta(loc_si3)
                  phi(loc_si1) = loc_phi(loc_si1)
                  phi(loc_si2) = loc_phi(loc_si2)
                  phi(loc_si3) = loc_phi(loc_si3)            
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
              call MPI_SEND(m_loc,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,38,MPI_COMM_WORLD,ierr)
              call MPI_SEND(loc_theta,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,48,MPI_COMM_WORLD,ierr)
              call MPI_SEND(loc_phi,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,58,MPI_COMM_WORLD,ierr)
          end if         
          
          !!! synchronize all the processes 
          !!! send the updated m configurations
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_BCAST(m,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          
          !!! send the updated theta configurations
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_BCAST(theta,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          
          !!! send the updated phi configurations
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_BCAST(phi,ns_unit*n_sites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

          !! initializing the most updated hamiltonian using the updated
          !! monte carlo configurations of m,theta and phi 
          hamiltonian(:,:) = cmplx(0.0,0.0)
          hamiltonian(:,:) = ham_noint(:,:)
          call ham_init()

        !!! loop end for all the splits
        end do    
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        !! printing the results in an output file after every meas_skip mc cycles
        if (my_id==0) then 
        if (mod(i,meas_skip)==0) then
          
          print *,i,tvar,fname_conf

          open(16,file=fname_conf,action='write',position='append')
          do j=0,n_sites-1,1
            suc0 = 0 + (j*ns_unit) ;suc1 = suc0+1 ;suc2 = suc0+2 ;  
            write(16,20) i,j,m(suc0),m(suc1),m(suc2),theta(suc0),theta(suc1),theta(suc2),&
            phi(suc0),phi(suc1),phi(suc2)
            20  format(I4,2X,I4,2X,ES22.8,2X,ES22.8,2X,ES22.8,2X,ES22.8,2X,ES22.8,2X,ES22.8,2X,ES22.8,2X,ES22.8,2X,ES22.8)
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
      open(22,file=fname_meast,action='write',position='append')
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
      print *,tvar,"finished"
      !!! lower the temperature of the system
      tvar = tvar-dtemp
    end do

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call MPI_FINALIZE(ierr)

end program ptca_repulsive
!!-----------------------------------------------------------------------------!!
!!-----------------------------------------------------------------------------!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! setting up the neighbour table!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine neighbour_table()
use varmodule
implicit none
    integer(8) :: i,xi,yi,ri,li,ui,di

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
    end do

end subroutine neighbour_table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! ---------------------------------------------------------------------------!!
!! -------------------- initialize monte carlo variables----------------------!!
!! ---------------------------------------------------------------------------!!
subroutine mcvar_init()
use varmodule
implicit none
  
  integer(8) :: i,j,si
  real(8):: rand1,rand2,rand3,rand_int3,rand_int1,rand_int2
 
  !! initialize the initial configuration of the monte carlo variables
  !! loop over the sites in the lattice
  do i=0,n_sites-1,1

    !! loop over the sites in the unit cell
    do j=0,ns_unit-1,1
    call random_number(rand1)
    call random_number(rand2)
    call random_number(rand3)
    
    si = (i*ns_unit) + j
    rand_int1=int(rand1*5000)
    m(si)=((((m_min)**3)+(((m_max)**3)-((m_min)**3)))*(rand_int1/5000.0))**(1.0_8/3.0_8)

    rand_int2=int(rand2*1000)
    theta(si)=acos(-1.0_8+(rand_int2/500.0))

    rand_int3=int(rand3*2000)
    phi(si)=(2.0_8*pi)*(rand_int3/2000.0)
    
    charge_confs(si) = 1.5*u_int
    end do
  end do
end subroutine mcvar_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!initializing the hamiltonian!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ham_init()
use varmodule
implicit none
  
  integer :: i,id0,id1,id2,id0n,id1n,id2n
  real(8) ::mx0,mx1,mx2
  real(8) ::my0,my1,my2
  real(8) ::mz0,mz1,mz2
  
  !! setting hopping dependent part of hamiltonian
  !call  haminitinUnitcell()

  do i = 0,n_sites-1, 1
    !! sites in the unit cell
    id0 = (ns_unit*i) 
    id0n = id0 + (ns_unit*n_sites)
    id1 = id0+1 
    id1n = id1 + (ns_unit*n_sites)
    id2 = id0+2
    id2n = id2 + (ns_unit*n_sites)


    !!! mx for three sites in the unit cell
    mx0 = m(id0)  * cos(phi(id0)) * sin(theta(id0))
    mx1 = m(id1)  * cos(phi(id1)) * sin(theta(id1))
    mx2 = m(id2)  * cos(phi(id2)) * sin(theta(id2))

    !!! my for three sites in the unit cell
    my0 = m(id0)  * sin(phi(id0)) * sin(theta(id0))
    my1 = m(id1)  * sin(phi(id1)) * sin(theta(id1))
    my2 = m(id2)  * sin(phi(id2)) * sin(theta(id2))

    !! mz for three sites in the unit cell
    mz0 = m(id0)  * cos(theta(id0))
    mz1 = m(id1)  * cos(theta(id1))
    mz2 = m(id2)  * cos(theta(id2))
    
    !! diagonal element for site 0 for up and dn spins
    hamiltonian(id0,id0) =  -(-charge_confs(id0)) - (0.5*u_int)*mz0
    hamiltonian(id0n,id0n) = -(-charge_confs(id0)) + (0.5*u_int)*mz0

    !! diagonal element for site 1 for up and dn spins
    hamiltonian(id1,id1) =  -(-charge_confs(id1)) - (0.5*u_int)*mz1
    hamiltonian(id1n,id1n) = -(-charge_confs(id1)) + (0.5*u_int)*mz1

    !! diagonal element for site20 for up and dn spins
    hamiltonian(id2,id2) =  -(-charge_confs(id2)) - (0.5*u_int)*mz2
    hamiltonian(id2n,id2n) = -(-charge_confs(id2)) + (0.5*u_int)*mz2

    ! setting the updn and dnup components for site 0 in the unit cell
    hamiltonian(id0,id0n) = -(0.5*u_int)*cmplx(mx0,-my0)
    hamiltonian(id0n,id0) = -(0.5*u_int)*cmplx(mx0,my0)


    ! setting the updn and dnup components for site 1 in the unit cell
    hamiltonian(id1,id1n) = -(0.5*u_int)*cmplx(mx1,-my1)
    hamiltonian(id1n,id1) = -(0.5*u_int)*cmplx(mx1,my1)

    ! setting the updn and dnup components for site 2 in the unit cell
    hamiltonian(id2,id2n) = -(0.5*u_int)*cmplx(mx2,-my2)
    hamiltonian(id2n,id2) = -(0.5*u_int)*cmplx(mx2,my2)

    
    !print*,'h1',id0,id1,id2
    !print*,i,'m',m(0+ns_unit*i),m(1+i*ns_unit),m(2+ns_unit*i),n_sites,ns_unit
    !print*,'h2',id0n,id1n,id2n,theta(0,i),theta(1,32),theta(2,32)
    !print*,i,'th',theta(0+ns_unit*i),theta(1+i*ns_unit),theta(2+ns_unit*i)
    !print*,i,'phi',phi(0+ns_unit*i),phi(1+i*ns_unit),phi(2+ns_unit*i)
    
  end do
  
end subroutine ham_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! this function will set the hopping element. It will set the hopping between the sites
!!! within the unit cell

subroutine haminitinUnitcell()
use varmodule
implicit none  
  integer :: id0 ,id1,id2 ,id0n , id1n , id2n,i

  do i = 0,n_sites-1, 1
    id0 = (ns_unit*i) !! first site in the unit cell i
    id1 = (ns_unit*i)+ 1 !! second site in the unit cell i
    id2 = (ns_unit*i)+ 2  !! third site in the unit cell i
    id0n = id0+(ns_unit*n_sites)  !! index for the spin down in unit cell i position 0
    id1n = id1+(ns_unit*n_sites)  !! index for the spin down in unit cell i position 1
    id2n = id2+(ns_unit*n_sites)  !! index for the spin down in unit cell i position 2
    
    !!! hopping from site 0-->1 & 0--->2 for up and down spin
    hamiltonian(id0,id1) = -t_hopping; hamiltonian(id0n,id1n) = -t_hopping
    hamiltonian(id0,id2) = -t_hopping; hamiltonian(id0n,id2n) = -t_hopping

    ! conjugate hopping 1--->0 & 1--->2 for up and down spin
    hamiltonian(id1,id0) = -t_hopping; hamiltonian(id1n,id0n) = -t_hopping
    hamiltonian(id1,id2) = -t_hopping; hamiltonian(id1n,id2n) = -t_hopping

    ! conjugate hopping from 2-->0 and 2--->1 for up and down spin
    hamiltonian(id2,id0) = -t_hopping; hamiltonian(id2n,id0n) = -t_hopping
    hamiltonian(id2,id1) = -t_hopping; hamiltonian(id2n,id1n) = -t_hopping
    !print *,i,id0,id1,id2,id0n,id1n,id2n
  
  end do
!call haminitOutunitcell(right,left,up,down,hamiltonian)
call haminitOutunitcell()

end subroutine haminitinUnitcell
!!-----------------------------------------------------------------------------!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! initializing the hamiltonian for the hopping for sites within the outsite
subroutine haminitOutunitcell()

use varmodule
  implicit none
  integer :: si 
  integer :: si0, s0l, s0ld, si0n, s0ln, s0ldn ! sites related to 0 site in the unit cell 
  integer :: si1, s1r, s1d, si1n, s1rn, s1dn ! sites related to site 1 in the unit cell
  integer :: si2 , s2u, s2ur , si2n , s2un , s2urn ! sites related to site 2 in the unit cell

  
  do si=0,n_sites-1,1  
    !! 1st site in the unit cell and site it's connected to 2nd and 3rd site of other unit cell
    si0 = ns_unit*si 
    s0l = left(si) * ns_unit + 1 ; s0ld = left(down(si))*ns_unit + 2
    si0n = si0+(ns_unit*n_sites) ; s0ln = s0l+(ns_unit*n_sites) ; s0ldn = s0ld+(ns_unit*n_sites)

    !! 2nd site in the unit cell is connected to 1st and 3rd site in the other unit cell
    si1 = (ns_unit*si)+1 
    s1r = right(si) * ns_unit ; s1d = ns_unit*down(si) + 2
    si1n = si1+(ns_unit*n_sites); s1rn = s1r+(ns_unit*n_sites);s1dn = s1d+(ns_unit*n_sites)

    !! 3rd site in the unit cell is connected to 1st site of unit cell in up and 0 site 
    !! of unit cell in up right direction
    si2 = (ns_unit*si)+2
    s2u = (up(si)*ns_unit)+1 ; s2ur = right(up(si))*ns_unit 
    si2n = si2+(ns_unit*n_sites); s2un  = s2u + (ns_unit*n_sites) ; s2urn = s2ur+(ns_unit*n_sites)

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
!    print *,si,si0,si1,si2
!    print *,si,si0,si0n,s0ln
  end do

end subroutine haminitOutunitcell



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!initializing the cluster hamiltonian!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_ham(site_clster)
use varmodule
implicit none
  integer(8) :: i,si,x,y,xip,yip,xim,yim,sri,sui
  integer(8) :: site_clster ! site that has the variable that will be changed
  integer(8) :: sli,sdi,sild,siru
  integer(8) :: sic, sicn , clsi,clsin
  integer(8) :: sic1 , sicn1 , clsi1 , clsin1 
  integer(8) :: sic2 , sicn2 , clsi2 , clsin2 

  !! constructing the cluster hamiltonian
  do i = 0,cls_dim-1, 1
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

      !!! diagonal part of the cluster hamiltonian !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! setting up the up,down and down up terms for the 1st site in the unit cell.

      sic = si*ns_unit ; sicn = sic + (ns_unit*cls_dim)  !! index of a site 0 in the cluster for up and dn spin
      clsi = cl_st(site_clster,si)*ns_unit ; clsin = clsi +(ns_unit*n_sites) !! index of site 0 in the lattice for up and dn spin


      sic1 = sic + 1 ; sicn1 = sicn+1 !! index of site1 in the cluster for up and down spin
      clsi1 = clsi + 1; clsin1 = clsin+1 !! index of site1 in the lattice for up and down spin

      sic2 = sic1 + 1 ; sicn2 = sicn1+1   !! index of site2 in the clsuter for up and down spin
      clsi2 = clsi1 + 1; clsin2 = clsin1+1 !! index of site2 in the lattice for up and down spin

      !! setting the term between up,down for site 0 in the cluster from the hamiltonian
      hamil_cls(sic,sicn) = hamiltonian(clsi,clsin)  !! term that is mx+imy
      hamil_cls(sicn,sic) = hamiltonian(clsin,clsi)  !! term that is mx-imy

      !!! setting up term between up,down for site 1 in the cluster from hamiltonian
      hamil_cls(sic1,sicn1) = hamiltonian(clsi1,clsin1) 
      hamil_cls(sicn1,sic1) = hamiltonian(clsin1,clsi1) 

      !!! setting up the term between up,down for site 2 in the cluster from hamiltonian
      hamil_cls(sic2,sicn2) = hamiltonian(clsi2,clsin2)
      hamil_cls(sicn2,sic2) = hamiltonian(clsin2,clsi2)


      !!! diagonal part of the hamiltonian for upup and dndn for site 0
      hamil_cls(sic,sic)=hamiltonian(clsi,clsi) - mu  !! term that is mz
      hamil_cls(sicn,sicn)=hamiltonian(clsin,clsin) - mu !! term that is -mz

      !!! diagonal part of the hamiltonian for upup and dndn for site 1
      hamil_cls(sic1,sic1)=hamiltonian(clsi1,clsi1)  - mu 
      hamil_cls(sicn1,sicn1)=hamiltonian(clsin1,clsin1) - mu

      !!! diagonal part of the hamiltonian for upup and dndn for site 2
      hamil_cls(sic2,sic2)=hamiltonian(clsi2,clsi2) - mu
      hamil_cls(sicn2,sicn2)=hamiltonian(clsin2,clsin2) - mu
      
      
        !print *,si,site_clster,sic,sic1,sic2,clsi,clsi1,clsi2
        !print *,'h1',hamiltonian(clsi,clsi),hamiltonian(clsi1,clsi1),hamiltonian(clsi2,clsi2)
        !print *,'h2',hamiltonian(clsi,clsin),hamiltonian(clsi1,clsin1),hamiltonian(clsi2,clsin2)  
        !print *,'h2s',clsi,clsin,clsi1,clsin1,clsi2,clsin2
      
        !print *,clsi2,clsi2,clsin2,clsin2
      !print *,si,cl_st(site_clster,si),site_clster
      
    end do
    
end subroutine cluster_ham


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hamclsinitinUnitcell()
use varmodule
  implicit none
  
  integer :: i,x,y,si
  integer(8) :: nsi0,nsi1,nsi2
  integer(8) :: nsi0c,nsi1c,nsi2c
  
  do i = 0, cls_dim-1, 1
      x = mod(i,cls_sites) ! x index in the  site cluster x --> [0,1,cls_sites-1]
      y = i/cls_sites       ! y index in the  site cluster y --> [0,1,cls_sites-1]
      si = x+(cls_sites*y)  ! site produced in the  site cluster si-->[0,cls_dim-1]


      nsi0 = ns_unit*si ; nsi1 = nsi0+1 ; nsi2 = nsi0+2
      nsi0c = nsi0+(ns_unit*cls_dim) ; nsi1c = nsi1+(ns_unit*cls_dim) ;  nsi2c = nsi2+(ns_unit*cls_dim)
      !print *,nsi0,nsi1,nsi2,nsi0c,nsi1c,nsi2c
      
      !! 0--->1 & 0--->2 for spin up in the unit cell at si
      hamil_cls(nsi0,nsi1) = -t_hopping
      hamil_cls(nsi0,nsi2) = -t_hopping
      
      !! 1--->0 & 1---> 2 for spin up in the unit cell at si
      hamil_cls(nsi1,nsi0) = -t_hopping
      hamil_cls(nsi1,nsi2) = -t_hopping

      !! 2--->0 & 2---> 1 for spin up in the unit cell at si
      hamil_cls(nsi2,nsi0) = -t_hopping
      hamil_cls(nsi2,nsi1) = -t_hopping
      

      !! 0--->1 & 0--->2 for spin down in the unit cell at si
      hamil_cls(nsi0c,nsi1c) = -t_hopping
      hamil_cls(nsi0c,nsi2c) = -t_hopping
      
      !! 1--->0 & 1---> 2 for spin down in the unit cell at si
      hamil_cls(nsi1c,nsi0c) = -t_hopping
      hamil_cls(nsi1c,nsi2c) = -t_hopping

      !! 2--->0 & 2---> 1 for spin down in the unit cell at si
      hamil_cls(nsi2c,nsi0c) = -t_hopping
      hamil_cls(nsi2c,nsi1c) = -t_hopping
  end do
  call hamclsinitOutunitcell()

end subroutine hamclsinitinUnitcell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! this subroutine will initialize the hamiltonian for the bonds outside the unit cell
subroutine hamclsinitOutunitcell()
use varmodule
  implicit none
  integer :: si,i,x,y,sl,sr,su,sd ,sru,sld
  integer(8) :: si0 , s0l , s0ld , si0n , s0ln , s0ldn
  integer(8) :: si1, s1r , s1d , si1n , s1rn , s1dn
  integer(8) :: si2 , s2u ,s2ur , si2n , s2un , s2urn
  
  do i = 0, cls_dim-1, 1
      x = mod(i,cls_sites) ! x index in the  site cluster x --> [0,1,cls_sites-1]
      y = i/cls_sites       ! y index in the  site cluster y --> [0,1,cls_sites-1]
      si = x+(cls_sites*y)  ! site produced in the  site cluster si-->[0,cls_dim-1]
      sl = mod(x-1+cls_sites,cls_sites) + (y*cls_sites) !! left neighbour of the site si
      sr = mod(x+1,cls_sites) + (y*cls_sites)       !! right neighbour of the site si
      su = (mod(y+1,cls_sites)*cls_sites)+x !! up neighbour of the site si
      sd = (mod(y-1+cls_sites,cls_sites)*cls_sites) + x !! down neighbour of the site si
      sru = mod(x+1,cls_sites) + (mod(y+1,cls_sites)*cls_sites) !! right neighbour of site si
      sld = mod(x-1+cls_sites,cls_sites) + (mod(y-1+cls_sites,cls_sites)*cls_sites) !! left down neighbour of site si

      !! 1st site in the unit cell and site it's connected to 2nd and 3rd site of other unit cell
      si0 = ns_unit*si ; si0n = si0 + (ns_unit*cls_dim)
      s0l = (sl*ns_unit)+1 ; s0ln = s0l + (ns_unit*cls_dim)
      s0ld = (sld*ns_unit )+2  ; s0ldn = s0ld + (ns_unit*cls_dim)
      
      
      !! 2nd site in the unit cell is connected to 1st and 3rd site in the other unit cell
      si1 = si0+1 ; si1n = si1+(ns_unit*cls_dim)
      s1r = (sr*ns_unit) ; s1rn = s1r + (ns_unit*cls_dim)
      s1d = (sd*ns_unit)+ 2 ; s1dn = s1d + (ns_unit*cls_dim)

      !! 3rd site in the unit cell is connected to 1st site of unit cell in up and 0 site 
      !! of unit cell in up right direction
      si2 = si0+2 ; si2n = si2+(ns_unit*cls_dim)
      s2u = (su*ns_unit)+1 ; s2un = s2u + (ns_unit*cls_dim)
      s2ur = (sru*ns_unit) ; s2urn = s2ur + (ns_unit*cls_dim)

      !! matrix element between 0 site and neighbour for up and down
      hamil_cls(si0,s0l) = -t_hopping ; hamil_cls(si0n,s0ln) = -t_hopping
      hamil_cls(s0l,si0) = -t_hopping ; hamil_cls(s0ln,si0n) = -t_hopping

      hamil_cls(si0,s0ld) = -t_hopping ; hamil_cls(si0n,s0ldn) = -t_hopping
      hamil_cls(s0ld,si0) = -t_hopping ; hamil_cls(s0ldn,si0n) = -t_hopping
  
      !! matrix element between 1 site and neighbour
      hamil_cls(si1,s1r) = -t_hopping;  hamil_cls(si1n,s1rn) = -t_hopping
      hamil_cls(s1r,si1) = -t_hopping;hamil_cls(s1rn,si1n) = -t_hopping

      hamil_cls(si1,s1d) = -t_hopping ; hamil_cls(si1n,s1dn) = -t_hopping
      hamil_cls(s1d,si1) = -t_hopping ; hamil_cls(s1dn,si1n) = -t_hopping

      !! matrix element between 2nd site and neighbour
      hamil_cls(si2,s2u)  = -t_hopping ; hamil_cls(si2n,s2un) = -t_hopping
      hamil_cls(s2u,si2)  = -t_hopping ; hamil_cls(s2un,si2n) = -t_hopping

      hamil_cls(si2,s2ur) = -t_hopping ; hamil_cls(si2n,s2urn) = -t_hopping
      hamil_cls(s2ur,si2) = -t_hopping ; hamil_cls(s2urn,si2n) = -t_hopping
  
   end do 

end subroutine hamclsinitOutunitcell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! generate cluster sites!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_sites()
use varmodule
implicit none
    integer(8) :: i,j,sil
    integer(8) :: xi,yi,xn,yn
    integer(8) :: x_init,y_init

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
subroutine lattice_splt()
use varmodule
implicit none
  !! given the x coordinate(i) and y coordinate(j) 
  !! this subroutine can categorize the lattice into 
  !! certain number of  groups based on the size of the cluster
 
  integer(8) :: xi,yi  !! loop variables
  integer(8) :: i ,j
  integer(8) :: splt_var ! --> 0 , # of partitions of lattice
  integer(8) :: ki !! used in the array as an index
  integer (8) :: l_minus_1 

  l_minus_1 = L - 1

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
subroutine mc_sweep(site_clster)
use varmodule
implicit none
  integer(8) :: loc_site,site_clster!! local site that is being flipped

  
  ! position of global site in the cluster (center)
  loc_site = int((0.5*cls_sites)+cls_sites*(0.5*cls_sites))

  !!! updating first site in the unit cell
  call mcsweepUnitCell(loc_site,site_clster)
end subroutine mc_sweep
!! ---------------------------------------------------------------------------!!

!! ---------------------------------------------------------------------------!!
!! --------------  monte carlo sweep for the unit cell------------------------!!
!! ---------------------------------------------------------------------------!!

subroutine mcsweepUnitCell(loc_site,site_clster)
use varmodule
implicit none
  integer(8) :: loc_site,site_clster !! index of site in the cluster
  integer(8) :: sil0,silc0,silt0,siltn0
  integer(8) :: sil1,silc1,silt1,siltn1
  integer(8) :: sil2,silc2,silt2,siltn2

  real(8) :: beta,e_u,e_v,enr_loc,mc_prob

  !!! variables for the monte carlo procedure
  real(8) :: mx0,my0,mz0 ! mx,my,mz for site 0 in the unit cell 
  real(8) :: mx1,my1,mz1 ! mx,my,mz for site 1 in the unit cell 
  real(8) :: mx2,my2,mz2 ! mx,my,mz for site 2 in the unit cell 
  
  real(8) :: rand1, rand2 , rand3 ! random number to generate,m,theta,phi
  real(8) :: rand_int1,rand_int2, rand_int3
  real(8) :: delE
  real(8) :: tempm0,temptheta0,tempphi0 !! store temp value of m,theta,phi for site 0
  real(8) :: tempm1,temptheta1,tempphi1 
  real(8) :: tempm2,temptheta2,tempphi2 

  complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1) :: temp_clsham,temp2_clsham !! to store a copy of the cluster hamiltonian 
  real(8),dimension(0:ns_unit*n_sites-1) :: loc_m  !! local m array

  !!!! generate variables for mx0,my0,mz0
  call random_number(rand1)
  call random_number(rand2)
  call random_number(rand3)

  !! temp m
  rand_int1=rand1*5000
  tempm0=((((m_min)**3)+(((m_max)**3)-((m_min)**3)))*(rand_int1/5000.0))**(1.0_8/3.0_8)
  
  !! temp theta
  rand_int2=rand2*1000
  temptheta0=acos(-1.0_8+(rand_int2/500.0))
  
  !! temp phi
  rand_int3=rand3*2000
  tempphi0=(2.0_8*pi)*(rand_int3/2000.0)

  mx0 = tempm0 * cos(tempphi0) * sin(temptheta0)
  my0 = tempm0 * sin(tempphi0) * sin(temptheta0)
  mz0 = tempm0 * cos(temptheta0)
  
  !!!! generate variables for mx1,my1,mz1 
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

  mx1 = tempm1 * cos(tempphi1) * sin(temptheta1)
  my1 = tempm1 * sin(tempphi1) * sin(temptheta1)
  mz1 = tempm1 * cos(temptheta1)
  
  !!!! generate variables for mx2,my2,mz2
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

  mx2 = tempm2 * cos(tempphi2) * sin(temptheta2)
  my2 = tempm2 * sin(tempphi2) * sin(temptheta2)
  mz2 = tempm2 * cos(temptheta2)
  
  temp_clsham(:,:)  = hamil_cls(:,:)
  loc_m =  m
  
  !!! call diagonalization subroutine
  call zheevd('V','U', dim_clsh, temp_clsham, dim_clsh, egval, work, lwork, &
                                           rwork, lrwork, iwork, liwork, info)
  
  !!! call subroutine to calculate the energy
  call enr_calc(loc_m,enr_loc,site_clster)
  e_u = enr_loc
   
  !! copying the cluster hamiltonian 
  temp_clsham(:,:) = cmplx(0.0,0.0)
  temp_clsham(:,:) = hamil_cls(:,:)
  
  !!! position of sites in in a given unit cell
  sil0 = 0 + (loc_site * ns_unit);  sil1 = sil0 + 1  ;  sil2 = sil0 + 2
  silc0 = sil0 + (ns_unit*cls_dim);  silc1 = silc0+1;  silc2 = silc0 + 2

  !!! mapping site in the cluster to the lattice
  silt0 = 0+(site_clster*ns_unit) 
  silt1 = silt0 + 1 ; silt2 = silt0 + 2
  siltn0 = silt0 + (ns_unit*n_sites) !! positin of the cluster center for the down spin
  siltn1 = siltn0 + 1 ; siltn2 = siltn0 + 2

  !!! updating cluster hamiltonian for site 0
  temp_clsham(sil0,sil0) =  -(mu-charge_confs(silt0)) - (0.5*u_int)*mz0
  temp_clsham(silc0,silc0) = -(mu-charge_confs(silt0)) + (0.5*u_int)*mz0
  temp_clsham(sil0,silc0) =  -(0.5*u_int)*cmplx(mx0,-my0)
  temp_clsham(silc0,sil0) = -(0.5*u_int)*cmplx(mx0,my0)
  
  !!! updating cluster hamiltonian for site 1
  temp_clsham(sil1,sil1) =  -(mu-charge_confs(silt1)) - (0.5*u_int)*mz1
  temp_clsham(silc1,silc1) = -(mu-charge_confs(silt1)) + (0.5*u_int)*mz1
  temp_clsham(sil1,silc1) =  -(0.5*u_int)*cmplx(mx1,-my1)
  temp_clsham(silc1,sil1) = -(0.5*u_int)*cmplx(mx1,my1)
  
 !!! updating cluster hamiltonian for site 2
  temp_clsham(sil2,sil2) =  -(mu-charge_confs(silt2)) - (0.5*u_int)*mz2
  temp_clsham(silc2,silc2) = -(mu-charge_confs(silt2)) + (0.5*u_int)*mz2
  temp_clsham(sil2,silc2) =  -(0.5*u_int)*cmplx(mx2,-my2)
  temp_clsham(silc2,sil2) = -(0.5*u_int)*cmplx(mx2,my2)

  
  loc_m(silt0) = tempm0
  loc_m(silt1) = tempm1
  loc_m(silt2) = tempm2
  
  !!! call diagonalization subroutine
  call zheevd('V','U', dim_clsh, temp_clsham, dim_clsh, egval, work, lwork, &
                                           rwork, lrwork, iwork, liwork, info)
  
  !print *,site_clster,(site_clster*ns_unit)+si,egval(int(0.5*dim_clsh)+1),egval(int(0.5*dim_clsh)),egval(0),mu
  !print *,ns_unit,si,si+(ns_unit*site_clster),site_clster
  !!! call subroutine to calculate the energy
  call enr_calc(loc_m,enr_loc,site_clster)
  e_v = enr_loc
  
  !print *,'si:',si,site_clster,loc_site,si,sil,silc,silt,siltn
  
  !! calculate the energy difference
  delE = e_v - e_u
  beta = 1./tvar
  if (delE < 0.0) then
    m(silt0) = tempm0  !! update m 
    theta(silt0) = temptheta0 !! update theta
    phi(silt0) = tempphi0  !! update phi
    m(silt1) = tempm1  !! update m 
    theta(silt1) = temptheta1 !! update theta
    phi(silt1) = tempphi1  !! update phi
    m(silt2) = tempm2  !! update m 
    theta(silt2) = temptheta2 !! update theta
    phi(silt2) = tempphi2  !! update phi
      
    !!! updating cluster hamiltonian for site 0
    hamil_cls(sil0,sil0) =  -(mu-charge_confs(silt0)) - (0.5*u_int)*mz0
    hamil_cls(silc0,silc0) = -(mu-charge_confs(silt0)) + (0.5*u_int)*mz0
    hamil_cls(sil0,silc0) =  -(0.5*u_int)*cmplx(mx0,-my0)
    hamil_cls(silc0,sil0) = -(0.5*u_int)*cmplx(mx0,my0)
    
    !!! updating cluster hamiltonian for site 1
    hamil_cls(sil1,sil1) =  -(mu-charge_confs(silt1)) - (0.5*u_int)*mz1
    hamil_cls(silc1,silc1) = -(mu-charge_confs(silt1)) + (0.5*u_int)*mz1
    hamil_cls(sil1,silc1) =  -(0.5*u_int)*cmplx(mx1,-my1)
    hamil_cls(silc1,sil1) = -(0.5*u_int)*cmplx(mx1,my1)
    
    !!! updating cluster hamiltonian for site 2
    hamil_cls(sil2,sil2) =  -(mu-charge_confs(silt2)) - (0.5*u_int)*mz2
    hamil_cls(silc2,silc2) = -(mu-charge_confs(silt2)) + (0.5*u_int)*mz2
    hamil_cls(sil2,silc2) =  -(0.5*u_int)*cmplx(mx2,-my2)
    hamil_cls(silc2,sil2) = -(0.5*u_int)*cmplx(mx2,my2)

    !! updating the full hamiltonian with the new monte carlo variables  
    hamiltonian(silt0,silt0) = -(-charge_confs(silt0)) - (0.5*u_int)*mz0
    hamiltonian(silt1,silt1) = -(-charge_confs(silt0)) + (0.5*u_int)*mz0
    hamiltonian(silt0,siltn0) = -(0.5*u_int)*cmplx(mx0,-my0)
    hamiltonian(siltn0,silt0) = -(0.5*u_int)*cmplx(mx0,my0)

    hamiltonian(silt1,silt1) = -(-charge_confs(silt1)) - (0.5*u_int)*mz1
    hamiltonian(siltn1,siltn1) = -(-charge_confs(silt1)) + (0.5*u_int)*mz1
    hamiltonian(silt1,siltn1) = -(0.5*u_int)*cmplx(mx1,-my1)
    hamiltonian(siltn1,silt1) = -(0.5*u_int)*cmplx(mx1,my1)

    hamiltonian(silt2,silt2) = -(-charge_confs(silt2)) - (0.5*u_int)*mz2
    hamiltonian(siltn2,siltn2) = -(-charge_confs(silt2)) + (0.5*u_int)*mz2
    hamiltonian(silt2,siltn2) = -(0.5*u_int)*cmplx(mx2,-my2)
    hamiltonian(siltn2,silt2) = -(0.5*u_int)*cmplx(mx2,my2)



  else 
    !! if the energy is not lowered update the site with the probability exp(-beta*delE)
    call random_number(mc_prob)
    
    if (mc_prob < exp(-beta*delE)) then
        m(silt0) = tempm0  !! update m 
        theta(silt0) = temptheta0 !! update theta
        phi(silt0) = tempphi0  !! update phi
        m(silt1) = tempm1  !! update m 
        theta(silt1) = temptheta1 !! update theta
        phi(silt1) = tempphi1  !! update phi
        m(silt2) = tempm2  !! update m 
        theta(silt2) = temptheta2 !! update theta
        phi(silt2) = tempphi2  !! update phi
    
        !!! updating cluster hamiltonian for site 0
        hamil_cls(sil0,sil0) =  -(mu-charge_confs(silt0)) - (0.5*u_int)*mz0
        hamil_cls(silc0,silc0) = -(mu-charge_confs(silt0)) + (0.5*u_int)*mz0
        hamil_cls(sil0,silc0) =  -(0.5*u_int)*cmplx(mx0,-my0)
        hamil_cls(silc0,sil0) = -(0.5*u_int)*cmplx(mx0,my0)
        
        !!! updating cluster hamiltonian for site 1
        hamil_cls(sil1,sil1) =  -(mu-charge_confs(silt1)) - (0.5*u_int)*mz1
        hamil_cls(silc1,silc1) = -(mu-charge_confs(silt1)) + (0.5*u_int)*mz1
        hamil_cls(sil1,silc1) =  -(0.5*u_int)*cmplx(mx1,-my1)
        hamil_cls(silc1,sil1) = -(0.5*u_int)*cmplx(mx1,my1)
        
        !!! updating cluster hamiltonian for site 2
        hamil_cls(sil2,sil2) =  -(mu-charge_confs(silt2)) - (0.5*u_int)*mz2
        hamil_cls(silc2,silc2) = -(mu-charge_confs(silt2)) + (0.5*u_int)*mz2
        hamil_cls(sil2,silc2) =  -(0.5*u_int)*cmplx(mx2,-my2)
        hamil_cls(silc2,sil2) = -(0.5*u_int)*cmplx(mx2,my2)

        !! updating the full hamiltonian with the new monte carlo variables  
        hamiltonian(silt0,silt0) = -(-charge_confs(silt0)) - (0.5*u_int)*mz0
        hamiltonian(silt1,silt1) = -(-charge_confs(silt0)) + (0.5*u_int)*mz0
        hamiltonian(silt0,siltn0) = -(0.5*u_int)*cmplx(mx0,-my0)
        hamiltonian(siltn0,silt0) = -(0.5*u_int)*cmplx(mx0,my0)

        hamiltonian(silt1,silt1) = -(-charge_confs(silt1)) - (0.5*u_int)*mz1
        hamiltonian(siltn1,siltn1) = -(-charge_confs(silt1)) + (0.5*u_int)*mz1
        hamiltonian(silt1,siltn1) = -(0.5*u_int)*cmplx(mx1,-my1)
        hamiltonian(siltn1,silt1) = -(0.5*u_int)*cmplx(mx1,my1)

        hamiltonian(silt2,silt2) = -(-charge_confs(silt2)) - (0.5*u_int)*mz2
        hamiltonian(siltn2,siltn2) = -(-charge_confs(silt2)) + (0.5*u_int)*mz2
        hamiltonian(silt2,siltn2) = -(0.5*u_int)*cmplx(mx2,-my2)
        hamiltonian(siltn2,silt2) = -(0.5*u_int)*cmplx(mx2,my2)

    end if
  end if
end subroutine mcsweepUnitCell
!! ---------------------------------------------------------------------------!!




!! ---------------------------------------------------------------------------!!
!! ------------------------measurement of energy------------------------------!!
!! ---------------------------------------------------------------------------!!
subroutine enr_calc(loc_m,enr_loc,site_clster)
use varmodule
implicit none
  integer(8) :: i,si,site_clster
  real(8),dimension(0:ns_unit*n_sites-1)::loc_m
  real(8)::m_i
  real(8):: fac_be ! product of beta and e_i
  real(8) :: sum_e,sum_cl,beta,enr_loc,val

  
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
    m_i = loc_m(0+(si*ns_unit)) + loc_m(1+(si*ns_unit))+loc_m(2+(si*ns_unit))
    sum_cl = sum_cl + (m_i)*(m_i)
  end do
    enr_loc = sum_e + (sum_cl)*(0.25*u_int)
end subroutine enr_calc

!! ---------------------------------------------------------------------------!!

!!-----------------------------------------------------------------------------!!
!!-------------------- printing -----------------------------------------------!!
!!-----------------------------------------------------------------------------!!
subroutine print_f()
use varmodule
implicit none

   character(len=50) :: format_L
   character(len=50) :: format_U
   character(len=50) :: format_temp="(F8.6)"
   character(len=50) :: format_cls="(I1)"
   character(len=20):: str_1
   character(len=20):: str_2
   character(len=20):: str_3
   character(len=20)::str_4   

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

  fname_conf=trim('mcconfigurations_L')//trim(str_1)//trim('_temp')//trim(str_2)//trim('_Uint')&
                     //trim(str_3)//trim('_cluster')//trim(str_4)//'.dat'

  fname_eqt=trim('totalEquiltime_L')//trim(str_1)//trim('_temp')//trim(str_2)//trim('_Uint')&
                     //trim(str_3)//trim('_cluster')//trim(str_4)//'.dat'

  fname_meast=trim('totalmeastime_L')//trim(str_1)//trim('_temp')//trim(str_2)//trim('_Uint')&
                     //trim(str_3)//trim('_cluster')//trim(str_4)//'.dat'

  fname_mu=trim('mu_L')//trim(str_1)//trim('_temp')//trim(str_2)//trim('_Uint')&
                     //trim(str_3)//trim('_cluster')//trim(str_4)//'.dat'

!   print*,fname

 end subroutine print_f
!!-----------------------------------------------------------------------------!!
