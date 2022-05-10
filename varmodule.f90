!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!-------------- module to store the variables used during the monte carlo procedure -----------!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module  varmodule
    implicit none
    integer(8),parameter :: ns_unit = 3
    integer(8),parameter :: L = 20 !! system size
    integer(8),parameter :: n_sites = L * L !! number of sites in the lattice
    integer(8),parameter :: cls_sites =  6 !! cluster size
    integer(8),parameter :: ncl_by2 = 0.5*(cls_sites)+1 !! dividing cls_sites by 2
    integer(8),parameter :: n_splits = (ncl_by2)*(ncl_by2)
    integer(8),parameter :: split_sites = n_sites/n_splits
    integer(8),parameter :: cls_dim = (cls_sites)*(cls_sites) !! number of sites in the cluster
    integer(8),parameter :: n_equil  = -1 !! no of equilibrium steps
    integer(8),parameter :: n_meas  = -1 !! no of measurements
    integer(8),parameter :: meas_skip = 10 ! make measurement after this mc cycles
    integer(8),parameter :: dim_h = 6*n_sites  ! dimensionality of hamiltonian
    integer(8),parameter :: dim_clsh = 6*cls_dim ! dimensionality of cluster hamiltonian
    integer(8),parameter :: mu_cnf = 25 !! monte carlo steps to set mu
    real(8),parameter :: pi = 4*atan(1.0)
    real(8),parameter :: t_hopping = 1.0
    real(8),parameter :: m_max = 2.0_8
    real(8),parameter :: m_min = 0.0_8
    real(8) :: mu  = 0.0 !! chemical potential
    real(8),parameter :: temp = 0.30  !! simulation temperature
    real(8),parameter :: dtemp = 0.01 !! temperature step to lower the temperature
    real(8),parameter :: t_min = 0.28 !! minimum temperature for the simulation
    real(8)  :: tvar !! variable used for varying tempreature
    real(8) :: u_int = 0.0 !! interaction strength

    !!! this array will be initialized to -1 at the starting 
    !!! entry will be changed to 1 when that particular site is updated during mc
    integer(8),dimension(0:split_sites-1)::changed_ids,loc_ids
    real(8),dimension(0:dim_h-1) :: egval_fl !! eigenvalues for the full hamiltonian
    real(8),dimension(0:dim_clsh-1) :: egval !! eigenvalues for the cluster hamiltonian
    
    !!!!!!!!!!!!!!!!!! initialize neighbour table !!!!!!!!!!!!!!!!!!!!!
    integer(8),dimension(0:n_sites-1):: right,left,up,down,right_up,left_down

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! array that will hold all the sites in the cluster!!!
    integer(8), dimension(0:n_sites-1,0:cls_dim-1)::cl_st ! sites in the cluster at site j
    integer(8),dimension(0:n_splits-1,0:split_sites-1) :: sites_array !! array to store information about the split sites

    !!!!!!!!!!!!!!!!!! variational parameters of the monte carlo procedure !!!!!!!!
    real(8),dimension(0:(ns_unit*n_sites)-1):: m,m_loc !!  m for master, m_loc = local
    real(8),dimension(0:(ns_unit*n_sites)-1):: theta,loc_theta !! theta
    real(8),dimension(0:(ns_unit*n_sites)-1):: phi,loc_phi  !! phi
    real(8),dimension(0:(ns_unit*n_sites)-1):: charge_confs !! to store charge configurations
    
    !!!!!!!!!!!!!!! full hamiltonian and cluster hamiltonian !!!!!!!!!!!!
    complex(8),dimension(0:dim_h-1,0:dim_h-1) :: hamiltonian !! full hamiltonian
    complex(8),dimension(0:dim_h-1,0:dim_h-1) :: ham_noint  !! non interacting hamiltonian
    complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1) :: hamil_cls !! cluster hamiltonian
    complex(8),dimension(0:dim_clsh-1,0:dim_clsh-1) :: hamcls_noint !! non interacting cluster hamiltonian
    
    !!!!!!!!!!!!!!! for lapack subroutines for cluster hamiltonian !!!!!!!!!!!!!!!
    integer(8),parameter :: lwork  = (2*dim_clsh)+(dim_clsh**2)
    integer(8),parameter :: lrwork = 2*(dim_clsh**2)+(5*(dim_clsh)+1)
    integer(8),parameter :: liwork = (5*dim_clsh)+3
    integer(8) :: info = 10
   
    complex(8),dimension(lwork)::work
    real(8),dimension(lrwork) :: rwork
    integer(8),dimension(liwork) :: iwork

    !!!!!!!!!!!!!!! for lapack subroutines for cluster hamiltonian !!!!!!!!!!!!!!!
    integer(8),parameter :: lwork_full  = (2*dim_h)+(dim_h**2)
    integer(8),parameter :: lrwork_full = 2*(dim_h**2)+(5*(dim_h)+1)
    integer(8),parameter :: liwork_full = (5*dim_h)+3
    
   
    complex(8),dimension(lwork_full)::work_full
    real(8),dimension(lrwork_full) :: rwork_full
    integer(8),dimension(liwork_full) :: iwork_full
      
end module varmodule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!