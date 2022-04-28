! gfortran -llapack -lblas lieb_lattice_heating_sample1.f90 

program lieb

implicit none

real(8),parameter :: pi=3.1415926535897932_8
real(8),parameter :: U=6.0_8
real(8),parameter :: m_max=2.0_8,m_min=0.0_8

integer(8) :: a,b,i,j,l1,x,y,z,mc_sweep,temp_sweep,steps,skip,info,count_filling
integer(8),parameter :: sites=3,L=18,n=sites*(L**2),m=2*n
integer(8) :: count1,config,count2,steps_mu_fix

complex(8),dimension(m,m) :: hamil
real(8),dimension(L,L,sites) :: m_config,theta_config,phi_config,charge_config

integer(8),parameter :: L_cl=4,n_cl=sites*(L_cl**2),m_cl=2*n_cl
integer(8),parameter :: lwork_cl=2*m_cl+m_cl**2,liwork_cl=5*m_cl+3,lrwork_cl=2*(m_cl**2)+5*m_cl+1
integer(8),parameter :: lwork=2*m+m**2,liwork=5*m+3,lrwork=2*(m**2)+5*m+1

real(8) :: dummy,dummy1,dummy2,rand_int,rand_int1,rand_int2,rand_int3,rand_int4,rand_int5
real(8) :: rand,rand1,rand2,rand3,rand4,rand5,rand6,mu,t1,t2,tprime,temp,beta,quad_energy,quad_energy1
real(8) :: sum1,sum2,value1,value2,engy1,engy2,Eelec1,Eelec2,F1,F2,delta_F,prb,mu_initial,mu_final
real(8) :: sum_mu_initial

complex(8),dimension(m_cl,m_cl) :: hamil_cl
real(8),dimension(L_cl,L_cl,sites) :: m_config_cl,theta_config_cl,phi_config_cl,charge_config_cl
complex(8),dimension(lwork_cl) :: work_cl
real(8),dimension(lrwork_cl) :: rwork_cl
real(8),dimension(m_cl) :: Eigval_cl
integer(8),dimension(liwork_cl) :: iwork_cl

complex(8),dimension(lwork) :: work
real(8),dimension(lrwork) :: rwork
real(8),dimension(m) :: Eigval
integer(8),dimension(liwork) :: iwork


                                  open(15,file='configurations_mag_lieb_kagome_L18_U6.0_tprime1.0.dat',status='UNKNOWN')
                                  open(16,file='temperature_mag_lieb_kagome_L18_U6.0_tprime1.0.dat',status='UNKNOWN')
                                  open(20,file='mu_optimized_mag_lieb_kagome_L18_U6.0_tprime1.0.dat',status='UNKNOWN')
                                     
!------------------------------------------------------------------------------------------
!Parameters
!------------------------------------------------------------------------------------------
t1=1.0_8
t2=1.0_8
tprime=1.0_8
!------------------------------------------------------------------------------------------

                                    do i=1,m_cl
                                       do j=1,m_cl
                                          hamil_cl=cmplx(0.0,0.0)
                                       enddo
                                    enddo   

!-----------------------------------------------------------------------------------------              
!Initial configurations                                                                         
!-----------------------------------------------------------------------------------------
                                    do i=1,L
                                       do j=1,L
                                          do l1=1,sites
                                             
                                             call random_number(rand)
                                             call random_number(rand1)
                                             call random_number(rand2)

                                             rand_int=int(rand*5000)
                                             m_config(i,j,l1)=((((m_min)**3)+(((m_max)**3)-((m_min)**3)))*&
                                                  (rand_int/5000.0))**(1.0_8/3.0_8)
                                             rand_int1=int(rand1*1000)
                                             theta_config(i,j,l1)=acos(-1.0_8+(rand_int1/500.0))
                                             rand_int2=int(rand2*2000)
                                             phi_config(i,j,l1)=(2.0_8*pi)*(rand_int2/2000.0)
                                             charge_config(i,j,l1)=1.5_8*U        !3*U/2
                                          enddo
                                       enddo
                                    enddo   
!-----------------------------------------------------------------------------------------              
!Temperature                                                          
!-----------------------------------------------------------------------------------------
                                    
                                    temp=0.21_8
                                    
                                    do temp_sweep=1,100
                                       config=0.0_8                                     
                                       temp=temp-0.01                                                    
                                     if(temp<0.001_8)then
                                       exit
                                     endif
                                         
                                       beta=1/temp
                                       steps=5000
                                       skip=3000
                                       count1=10

                                       steps_mu_fix=500
                                       count2=25

                                       sum_mu_initial=0.0_8
                                       count_filling=1

!-----------------------------------------------------------------------------------------              
!Monte Carlo Sweeps                                                                   
!----------------------------------------------------------------------------------------- 
                                       
                                       do mc_sweep=1,steps

                                          do x=1,L          !Site-loop
                                             do y=1,L        !Site-loop
                                                do z=1,sites !Site-loop

                          if((mc_sweep<steps_mu_fix).and.mod((mc_sweep+count2), count2)==0)then

                           call hamiltonian_lattice(hamil,sites,L,n,m,m_config,theta_config,phi_config,charge_config,&
                                t1,t2,tprime,U)    
                            info =10
                            call zheevd('V','U', m, hamil, m, Eigval, work, lwork,rwork,lrwork,iwork,liwork, info)
                                  
                            mu_initial=(Eigval(0.5*m+1)+Eigval(0.5*m))/2.0_8
                            sum_mu_initial = sum_mu_initial + mu_initial
                            mu=mu_initial

                            if(mc_sweep == steps_mu_fix)then
                            mu_final = sum_mu_initial/20
                            mu=mu_final
                         endif                       
                      endif
                            

!-----------------------------------------------------------------------------------------
!Cluster mapping
!-----------------------------------------------------------------------------------------

                                                   do i=1,L_cl
                                                      do j=1,L_cl
                                                        do l1=1,sites
                                                           a=(x+i-0.5_8*L_cl-1)
                                                           b=(y+j-0.5_8*L_cl-1)

                                                           if(a.lt.(1)) a = a + L
                                                           if(a.gt.(L)) a = a - L
                                                           if(b.lt.(1)) b = b + L
                                                           if(b.gt.(L)) b = b - L

                                                           m_config_cl(i,j,l1)=m_config(a,b,l1)
                                                           theta_config_cl(i,j,l1)=theta_config(a,b,l1)
                                                           phi_config_cl(i,j,l1)=phi_config(a,b,l1)
                                                           charge_config_cl(i,j,l1)=charge_config(a,b,l1)                          

                                                         enddo
                                                       enddo
                                                     enddo  
                                                   
!-----------------------------------------------------------------------------------------              
!Classical cost                                                                                
!----------------------------------------------------------------------------------------- 

                                                   quad_energy=0.0_8
                                                   
                                                   do i=1,L_cl
                                                      do j=1,L_cl
                                                         do l1=1,sites
                                                            quad_energy=quad_energy+real(m_config_cl(i,j,l1)**2)
                                                         enddo
                                                      enddo
                                                   enddo   

                                 call hamiltonian(hamil_cl,sites,L_cl,n_cl,m_cl,m_config_cl,theta_config_cl,phi_config_cl,&
                                      charge_config_cl,t1,t2,tprime,mu,U)
                                           info =10
                                         call zheevd('V','U', m_cl, hamil_cl, m_cl, Eigval_cl, work_cl, lwork_cl,&
                                                        rwork_cl, lrwork_cl, iwork_cl, liwork_cl, info)
                                                   
                                                   sum1=0.0_8
                                                   do i=1,m_cl
                                                      value1=(-beta*Eigval_cl(i))
                                                      if(value1>40)then
                                                         engy1=value1
                                                      else
                                                         engy1=log(1+exp(value1))
                                                      endif
                                                      sum1=sum1+engy1
                                                   enddo
                                                   Eelec1=(-(sum1)/beta)
                                                                                                      
!----------------------------------------------------------------------------------------
!Free energy before update
!----------------------------------------------------------------------------------------
                                                   F1=Eelec1+0.25_8*U*quad_energy

                                                   dummy=m_config(x,y,z)
                                                   dummy1=theta_config(x,y,z)
                                                   dummy2=phi_config(x,y,z)
                                                   
!----------------------------------------------------------------------------------------
!Update
!----------------------------------------------------------------------------------------
                                                  
                                                   call random_number(rand3)
                                                   call random_number(rand4)
                                                   call random_number(rand5)

                                                   rand_int3=int(rand3*5000)
                                                   m_config(x,y,z)=((((m_min)**3)+(((m_max)**3)-&
                                                        ((m_min)**3)))*(rand_int3/5000.0))**(1.0_8/3.0_8)
                                                   rand_int4=int(rand4*1000)
                                                   theta_config(x,y,z)=acos(-1.0_8 +(rand_int4/500.0))
                                                   rand_int5=int(rand5*2000)
                                                   phi_config(x,y,z)=(2.0_8*pi)*(rand_int5/2000.0)
                                                   charge_config(x,y,z)=1.5_8*U     !3*U/2
!----------------------------------------------------------------------------------------------
!Cluster mapping
!----------------------------------------------------------------------------------------------

                                                 do i=1,L_cl
                                                      do j=1,L_cl
                                                         do l1=1,sites

                                                            a=(x+i-0.5_8*L_cl-1)
                                                            b=(y+j-0.5_8*L_cl-1)

                                                            if(a.lt.(1)) a = a + L
                                                            if(a.gt.(L)) a = a - L
                                                            if(b.lt.(1)) b = b + L
                                                            if(b.gt.(L)) b = b - L

                                                           m_config_cl(i,j,l1)=m_config(a,b,l1)
                                                           theta_config_cl(i,j,l1)=theta_config(a,b,l1)
                                                           phi_config_cl(i,j,l1)=phi_config(a,b,l1)
                                                           charge_config_cl(i,j,l1)=charge_config(a,b,l1)
                                                         enddo
                                                      enddo
                                                   enddo
!------------------------------------------------------------------------------------------------
!Classical cost
!------------------------------------------------------------------------------------------------
                                                   
                                                   quad_energy1=0.0_8
                                                   do i=1,L_cl
                                                      do j=1,L_cl
                                                         do l1=1,sites
                                                          quad_energy1=quad_energy1+real(m_config_cl(i,j,l1)**2)
                                                         enddo
                                                      enddo
                                                   enddo   

                                 call hamiltonian(hamil_cl,sites,L_cl,n_cl,m_cl,m_config_cl,theta_config_cl,phi_config_cl,&
                                      charge_config_cl,t1,t2,tprime,mu,U)
                                                   info =10
                                 call zheevd('V','U', m_cl, hamil_cl, m_cl, Eigval_cl, work_cl, lwork_cl,&
                                                        rwork_cl, lrwork_cl, iwork_cl, liwork_cl, info)
                                                   
                                                  sum2=0.0_8
                                                   do i=1,m_cl
                                                      value2=(-beta*Eigval_cl(i))
                                                      if(value2>40)then
                                                         engy2=value2
                                                      else
                                                         engy2=log(1+exp(value2))
                                                      endif
                                                      sum2=sum2+engy2
                                                   enddo
                                                   Eelec2=(-(sum2)/beta)
                                                   
!-------------------------------------------------------------------------------------------------------
!Free energy after update
!-------------------------------------------------------------------------------------------------------
                                                   F2=Eelec2+0.25_8*U*quad_energy1
!-------------------------------------------------------------------------------------------------------
!Metropolis
!-------------------------------------------------------------------------------------------------------
                                                   if(F2<F1)then
                                                      F1=F2
                                                   else
   
                                                   delta_F=F2-F1                                                 
                                                   prb=exp(-delta_F*beta)            
                                                   call random_number(rand6)
                                                   if(rand6>prb)then
                                                   m_config(x,y,z)=dummy
                                                   theta_config(x,y,z)=dummy1
                                                   phi_config(x,y,z)=dummy2                                     
      
                                                   else
                                                      F1=F2
                                                   endif
                                                 endif 
                                                enddo       !Site-loop
                                             enddo          !Site-loop
                                          enddo             !Site-loop
                                         

                                          if((mc_sweep>skip) .and. mod((mc_sweep+count1), count1)==0)then
                                             
                                             do i=1,L
                                                do j=1,L
                                                   do l1=1,sites
                                                      write(15,11)i,j,l1,m_config(i,j,l1),theta_config(i,j,l1),phi_config(i,j,l1),&
                                                           charge_config(i,j,l1)
11                                                    format(I4,2X,I4,2X,I4,2X,F15.10,2X,F15.10,2X,F15.10,2X,F15.10)
                                                   enddo
                                                enddo
                                             enddo                                               
                                             config=config+1
                                           endif   

                                       enddo    !MC-loop   

                            call hamiltonian_lattice(hamil,sites,L,n,m,m_config,theta_config,phi_config,charge_config,&
                                t1,t2,tprime,U)
                            info =10
                            call zheevd('V','U', m, hamil, m, Eigval, work, lwork,rwork,lrwork,iwork,liwork, info)
                                      do i=1,m
                                         if(Eigval(i).le.mu)then
                                        count_filling=count_filling+1
                                     endif
                                  enddo
                                     write(20,*)temp,mu,count_filling

                                       write(16,*)temp
                                    enddo       !Temperature-loop

end program lieb
!--------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------
subroutine hamiltonian(hamil_cl,sites,L_cl,n_cl,m_cl,m_config_cl,theta_config_cl,phi_config_cl,&
     charge_config_cl,t1,t2,tprime,mu,U)

implicit none

integer(8) :: i,j,l1,l2
integer(8) :: sites,L_cl,n_cl,m_cl,K0_cl,K1_cl,K2_cl
real(8) :: mu,t1,t2,tprime,U

complex(8),dimension(m_cl,m_cl) :: hamil_cl
real(8),dimension(L_cl,L_cl,sites) :: m_config_cl,theta_config_cl,phi_config_cl,charge_config_cl
real(8),dimension(L_cl) :: P1_cl,P2_cl

!Initialization
                              do i=1,m_cl
                                 do j=1,m_cl
                                    hamil_cl(i,j)=cmplx(0.0,0.0)
                                 enddo
                              enddo                                
!Periodic boundary condition                                                              
                              P2_cl(1)=L_cl
                               do i=1,L_cl                                  
                                  P1_cl(i)=i+1
                                  P1_cl(L_cl)=1
                               enddo   
                               do i=1,L_cl
                                  P2_cl(i)=i
                                  P2_cl(L_cl)=0
                               enddo   
!Hopping elements                                 
                               do i=1,L_cl
                                  do j=1,L_cl
                                                                        
                                     K0_cl=(sites*((i-1)*L_cl+j))-2
                                                                         
                                     hamil_cl(K0_cl,(K0_cl+1))=-t1
                                     hamil_cl(K0_cl,(K0_cl+2))=-t1
                                     hamil_cl((K0_cl+1),K0_cl)=-t1
                                     hamil_cl((K0_cl+2),K0_cl)=-t1

                                     hamil_cl((K0_cl+n_cl),(K0_cl+n_cl+1))=-t1
                                     hamil_cl((K0_cl+n_cl),(K0_cl+n_cl+2))=-t1
                                     hamil_cl((K0_cl+n_cl+1),(K0_cl+n_cl))=-t1
                                     hamil_cl((K0_cl+n_cl+2),(K0_cl+n_cl))=-t1

                                       do l1=1,sites                   !Additional index for Lieb
                                          K0_cl=(sites*((i-1)*L_cl+j))-2+(l1-1)
                                          hamil_cl(K0_cl,K0_cl)=-(mu-charge_config_cl(i,j,l1))-0.5*U*&
                                               m_config_cl(i,j,l1)*cos(theta_config_cl(i,j,l1))
                                          hamil_cl((K0_cl+n_cl),(K0_cl+n_cl))=-(mu-charge_config_cl(i,j,l1))+0.5*U*&
                                               m_config_cl(i,j,l1)*cos(theta_config_cl(i,j,l1))
!Pairing field elements                                          
                                          hamil_cl(K0_cl,(K0_cl+n_cl))=-0.5*U*m_config_cl(i,j,l1)*&
                                               sin(theta_config_cl(i,j,l1))*cmplx(cos(phi_config_cl(i,j,l1)),&
                                               -sin(phi_config_cl(i,j,l1)))
                                          hamil_cl((K0_cl+n_cl),K0_cl)=-0.5*U*m_config_cl(i,j,l1)*&
                                               sin(theta_config_cl(i,j,l1))*cmplx(cos(phi_config_cl(i,j,l1)),&
                                               sin(phi_config_cl(i,j,l1)))
!Additional hoppings due to Lieb sites
                                          do l2=1,sites
                                             K1_cl=(sites*(P2_cl(i)*L_cl+j))-2+(l2-1)
                                             K2_cl=(sites*((i-1)*L_cl+P1_cl(j)))-2+(l2-1)

                                             if((l1==2).and.(l2==1))then
                                                hamil_cl(K0_cl,K1_cl)=-t2
                                                hamil_cl(K1_cl,K0_cl)=-t2
                                                hamil_cl((K0_cl+n_cl),(K1_cl+n_cl))=-t2
                                                hamil_cl((K1_cl+n_cl),(K0_cl+n_cl))=-t2
                                             endif   
                                             if((l1==3).and.(l2==1))then
                                                hamil_cl(K0_cl,K2_cl)=-t2
                                                hamil_cl(K2_cl,K0_cl)=-t2
                                                hamil_cl((K0_cl+n_cl),(K2_cl+n_cl))=-t2
                                                hamil_cl((K2_cl+n_cl),(K0_cl+n_cl))=-t2
                                             endif                                                
                                             if((l1==2).and.(l2==3))then
                                                hamil_cl(K0_cl,K1_cl)=-tprime
                                                hamil_cl(K1_cl,K0_cl)=-tprime
                                                hamil_cl((K0_cl+n_cl),(K1_cl+n_cl))=-tprime
                                                hamil_cl((K1_cl+n_cl),(K0_cl+n_cl))=-tprime
                                             endif
                                             if((l1==3).and.(l2==2))then
                                                  hamil_cl(K0_cl,K2_cl)=-tprime
                                                  hamil_cl(K2_cl,K0_cl)=-tprime
                                                  hamil_cl((K0_cl+n_cl),(K2_cl+n_cl))=-tprime
                                                  hamil_cl((K2_cl+n_cl),(K0_cl+n_cl))=-tprime
                                             endif
                                          enddo   
                                       enddo   
                                   enddo
                               enddo   
                                                                                                       
end subroutine hamiltonian
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------

subroutine hamiltonian_lattice(hamil,sites,L,n,m,m_config,theta_config,phi_config,&
     charge_config,t1,t2,tprime,U)

implicit none

integer(8) :: i,j,l1,l2
integer(8) :: sites,L,n,m,K0,K1,K2
real(8) :: t1,t2,tprime,U

complex(8),dimension(m,m) :: hamil
real(8),dimension(L,L,sites) :: m_config,theta_config,phi_config,charge_config
real(8),dimension(L) :: P1,P2

!Initialization
                              do i=1,m
                                 do j=1,m
                                    hamil(i,j)=cmplx(0.0,0.0)
                                 enddo
                              enddo                                
!Periodic boundary condition                   
                               P2(1)=L
                               do i=1,L                                  
                                  P1(i)=i+1
                                  P1(L)=1
                               enddo   
                               do i=1,L
                                  P2(i)=i
                                  P2(L)=0
                               enddo   
!Hopping elements                                 
                               do i=1,L
                                  do j=1,L
                                                                        
                                     K0=(sites*((i-1)*L+j))-2
                                                                         
                                     hamil(K0,(K0+1))=-t1
                                     hamil(K0,(K0+2))=-t1
                                     hamil((K0+1),K0)=-t1
                                     hamil((K0+2),K0)=-t1

                                     hamil((K0+n),(K0+n+1))=-t1
                                     hamil((K0+n),(K0+n+2))=-t1
                                     hamil((K0+n+1),(K0+n))=-t1
                                     hamil((K0+n+2),(K0+n))=-t1

                                       do l1=1,sites                   !Additional index for Lieb
                                          K0=(sites*((i-1)*L+j))-2+(l1-1)
                                          hamil(K0,K0)=-(-charge_config(i,j,l1))-0.5*U*m_config(i,j,l1)*&
                                               cos(theta_config(i,j,l1))
                                          hamil((K0+n),(K0+n))=-(-charge_config(i,j,l1))+0.5*U*m_config(i,j,l1)*&
                                               cos(theta_config(i,j,l1))
!Pairing field elements                                          
                                          hamil(K0,(K0+n))=-0.5*U*m_config(i,j,l1)*sin(theta_config(i,j,l1))*&
                                               cmplx(cos(phi_config(i,j,l1)),-sin(phi_config(i,j,l1)))
                                          hamil((K0+n),K0)=-0.5*U*m_config(i,j,l1)*sin(theta_config(i,j,l1))*&
                                               cmplx(cos(phi_config(i,j,l1)),sin(phi_config(i,j,l1)))
!Additional hoppings due to Lieb sites
                                          do l2=1,sites
                                             K1=(sites*(P2(i)*L+j))-2+(l2-1)
                                             K2=(sites*((i-1)*L+P1(j)))-2+(l2-1)

                                             if((l1==2).and.(l2==1))then
                                                hamil(K0,K1)=-t2
                                                hamil(K1,K0)=-t2
                                                hamil((K0+n),(K1+n))=-t2
                                                hamil((K1+n),(K0+n))=-t2
                                             endif   
                                             if((l1==3).and.(l2==1))then
                                                hamil(K0,K2)=-t2
                                                hamil(K2,K0)=-t2
                                                hamil((K0+n),(K2+n))=-t2
                                                hamil((K2+n),(K0+n))=-t2                                                
                                             endif                      
                                             if((l1==2).and.(l2==3))then
                                                hamil(K0,K1)=-tprime
                                                hamil(K1,K0)=-tprime
                                                hamil((K0+n),(K1+n))=-tprime
                                                hamil((K1+n),(K0+n))=-tprime
                                             endif
                                             if((l1==3).and.(l2==2))then
                                                  hamil(K0,K2)=-tprime
                                                  hamil(K2,K0)=-tprime
                                                  hamil((K0+n),(K2+n))=-tprime
                                                  hamil((K2+n),(K0+n))=-tprime
                                             endif
                                          enddo   
                                       enddo   
                                   enddo
                               enddo   
                                                                                                       
end subroutine hamiltonian_lattice
