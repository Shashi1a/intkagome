
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as lng
import glob
import itertools
import seaborn as sns




### class to define the Lattice
class Lattice:
    def __init__(self, l_cell, n_cell,s_cell, hopping, n_sites,dim_h,mu,u_int):
        '''
        This function will initialize the object and set the
        parameters:
        parameters:
            self(object): object
            l_cell(int):  number of unit cell
            n_cell(int): total number of unit cells
            hopping(float): nn nearest neighbour hopping
            n_sites(int): number of sites in the lattice
            dim_h(int): dimensionality of the hamiltonian matrix
            mu(float): chemical potential
            u_int(float): interaction strength
        '''
        self.l_cell = l_cell
        self.n_cell = n_cell
        self.s_cell = s_cell
        self.hopping = hopping
        self.n_sites = n_sites
        self.dim_h = dim_h
        self.mu = mu
        self.u_int = u_int

    ## this function will create the square lattice for
    ## the unit cells
    def neighbour(self):
        '''
        This function will create the neighbour table with
        periodic boundary conditions
        parameters:
            None
        return:
            None
        '''
        def sip(x): return (x+1) % self.l_cell
        def sim(x): return (x-1+self.l_cell) % self.l_cell
        ri, li, ui, di, rui, dli = [], [], [], [], [], []
        for j in range(self.n_cell):
            yi, xi = divmod(j, self.l_cell)
            ri.append(sip(xi)+yi*self.l_cell)
            li.append(sim(xi)+yi*self.l_cell)
            ui.append(xi+sip(yi)*self.l_cell)
            di.append(xi+sim(yi)*self.l_cell)
            rui.append(sip(xi)+sip(yi)*self.l_cell)
            dli.append(sim(xi)+sim(yi)*self.l_cell)
        self.right = np.array(ri, dtype='int')
        self.left = np.array(li, dtype='int')
        self.up = np.array(ui, dtype='int')
        self.down = np.array(di, dtype='int')
        self.right_up = np.array(rui, dtype='int')
        self.left_down = np.array(dli, dtype='int')



## defining the class for the hamiltonian 
## it uses lattice as the parent class
class Hamiltonian(Lattice):
    ## initializing the object
    def __init(self,l_cell,n_cell,s_cell,hopping,n_sites,dim_h,mu,u_int):
        
        ## using init of the parent class
        super(Hamiltonian, self).__init__(l_cell, n_cell, s_cell, hopping,n_sites, dim_h, mu, u_int)
                                          
    ## construct the hamiltonian matrix
    def hamunitcell(self):
        '''
        This function will initialize the elements of the matrix only within the unit cell.
        Hopping outside the unit cells is not allowed
        0-->1, 0-->2, 1-->0,1--->2, 2--->0, 2-->1
        parameters:
            self(object): An instance of the class
        return :
            hammat(matrix): the incomplete hamiltonian matrix
        '''
        hammat = np.zeros((self.dim_h,self.dim_h),dtype='complex')

        for i in range(self.n_cell):
            si0 = int(self.s_cell*i) ; si0n = int(si0 + (self.s_cell*self.n_cell))
            si1 = int(si0 + 1) ; si1n = int(si1 + (self.s_cell*self.n_cell))
            si2 = int(si0+2) ; si2n = int(si2 + (self.s_cell*self.n_cell))
            
            #print(f'si0: {si0} si0n: {si0n}')
            ## setting the bonds inside the unit cell
            hammat[si0, si1] = -self.hopping # up spin
            hammat[si0n, si1n] = -self.hopping # down spin
            hammat[si0, si2] = -self.hopping  # up spin
            hammat[si0n, si2n] = -self.hopping # down spin

            hammat[si1,si0] = -self.hopping # up spin
            hammat[si1n, si0n] = -self.hopping # down spin
            hammat[si1n,si2n] = -self.hopping # up spin
            hammat[si1n, si2n] = -self.hopping # down spin

            hammat[si2, si0] = -self.hopping # up spin
            hammat[si2n, si0n] = -self.hopping # down spin
            hammat[si2, si1] = -self.hopping # up spin
            hammat[si2n, si1n] = -self.hopping  # down spin

        return hammat

    ### initialize hopping outside the unit cells
    def hamOutUnitCell(self,hammat):
        '''
        This function will initialize hopping outside the unit cells
        parameters:
            self(object): Instance of the class
        return:
            None
        '''
        print(self.n_cell,self.n_sites)
        for i in range(self.n_cell):
            
            ## sites in a given unit cell
            s0,s1,s2 = int(self.s_cell*i),int(self.s_cell*i)+1,int(self.s_cell*i)+2
            s0n,s1n,s2n = s0+int(self.s_cell*self.n_cell),s1+int(self.s_cell*self.n_cell) , s2+int(self.s_cell*self.n_cell) 
            
            ### right,left,up,down neighbours of the unitcell
            ri, li, ui, di = int(self.right[i]), int(self.left[i]), int(self.up[i]), int(self.down[i])
            
            ## allowed neighbours of the site 0 of the unit cell
            s0l = int(self.s_cell*li+1); s0ld = int(self.s_cell*self.left_down[i]+2)
            s0ln = int(s0l+(self.n_cell*self.s_cell));s0ldn = int(s0ld + (self.n_cell*self.s_cell))

            ## allowed neighbours of the site 1 of the unit cell
            s1r = int(self.s_cell*ri); s1d = int(self.s_cell*di+2)
            s1rn = int(s1r + (self.s_cell*self.n_cell)) ; s1dn = int(s1d + (self.s_cell*self.n_cell))

            ## allowed neighbours of the site 2 of the unit cell
            s2u = int(self.s_cell*ui+1); s2ur = int(self.s_cell*self.right_up[i])
            s2un = int(s2u +(self.s_cell * self.n_cell)) ; s2urn = int(s2ur + (self.s_cell*self.n_cell))

            ### setting up the bonds outside the unitcell
            ## for the 1st site in the unitcell

            hammat[s0,s0l] = -self.hopping; hammat[s0n,s0ln] = -self.hopping
            hammat[s0l,s0] = -self.hopping ; hammat[s0ln,s0n] = -self.hopping

            hammat[s0,s0ld] = -self.hopping ; hammat[s0n,s0ldn] = -self.hopping
            hammat[s0ld,s0] = -self.hopping ; hammat[s0ldn,s0n] = -self.hopping

            ## for the 2nd site in the unit cell
            hammat[s1,s1r] = -self.hopping ; hammat[s1n,s1rn] = -self.hopping
            hammat[s1r,s1] = -self.hopping ; hammat[s1rn,s1n] = -self.hopping

            hammat[s1,s1d] = -self.hopping ; hammat[s1n,s1dn] = -self.hopping
            hammat[s1d,s1] = -self.hopping ; hammat[s1dn,s1n]  = -self.hopping

            ## for the 3rd site in the unit cell
            hammat[s2,s2u] = -self.hopping  ; hammat[s2n,s2un] = -self.hopping
            hammat[s2u,s2] = -self.hopping ; hammat[s2un,s2n]  = -self.hopping

            hammat[s2,s2ur] = -self.hopping ; hammat[s2n,s2urn] = -self.hopping
            hammat[s2ur,s2] = -self.hopping ; hammat[s2urn,s2n] = -self.hopping

        self.ham = hammat

    ### write the full hamiltonian
    def haminit(self):
        '''
        This function will construct the full hamiltonian in two steps:
        1. In first step it will sep up the hamiltonian for sites in the unit cell
        2. In the second step the hopping outside the unit cell will be set up.
        parameters:
            self(object): The class instance
        return:
            None
        '''
        mat = self.hamunitcell()
        self.hamOutUnitCell(mat)

    ### diagonalize the hamiltonian
    def diag(self):
        '''
        This function will diagonalize the hamiltonian
        parameters:
            self(object): Instance of the class
        return:
            evals (float): Eigenvalues of the hamiltonian
            evecs(float): eigenvectors of the hamiltonian
        '''
        evals,evecs = lng.eigh(self.ham)
        return evals,evecs


### class to obtain the files in the folder that we want to use
class Files:
    """
    This is a class method to obtain the files inside the folder.
    Since we are using the class method we don't have to create an 
    instance of the class to be able to use this method
    """
    @classmethod
    def getFiles(cls,file_path,pattern):
        """
        This function will get the files inside the golder given the file_path and the
        pattern we want to match
        parameter:
            cls (class)
        return:
            file_dict(dictionary): dictionary containing the temperature as the key
                                    and value as the name of the file.    
        """
        file_list = glob.glob(file_path+pattern)
        file_list = file_list
        n_files = len(file_list)

        ## temporary variable to store the key and value pair
        ki,vi  = [],[] # list to store the key and values
        for j in file_list:
            ki.append(j.split('temp')[1].split('_Uint')[0])
            vi.append(j)
        file_dict = dict(zip(np.array(ki,dtype='float'),vi))
        return file_dict



### class definition to incorporate the interacting part to the system
class Interaction(Hamiltonian):
    ###Constructor for the object
    def __init__(self,l_cell,n_cell,s_cell,hopping,n_sites,dim_h,mu,u_int):
        """
        Constructor for the class.
        parameters:
            self(object): class instance
            l_cell(integer): linear dimension of the lattice
            n_cell(integer): number of  unit cells in the lattice (L*L)
            s_cell(integer): number of sites in the unit cell
            hopping(float): hopping strength
            n_sites(integer): number of sites in the system (n_cell*s_cell)
            dim_h(integer) : dimensionaity of the hamiltonian considering the spin up and spin down particles
            u_int(float): interaction strength
            file_path(string): location where the file is stores
        return:
            None
        """

        super(Interaction, self).__init__(l_cell, n_cell,
                                          s_cell, hopping, n_sites, dim_h, mu, u_int)
        
        
    ### storing data of monte carlo variables into a dataframe 
    def datadf(self,file_t,file_dict):
        """
        This function will load the data and put them into a dataframe
        parameters:
            self(object): instance of the class
            file_t(float): key/temperature of the data we want to use
        return:
            datadf(dataframe): dataframe holding data for a particular temperature
        """

        data = pd.read_csv(file_dict[file_t], header=None)
        data.columns = ['data']

        
        ### number of monte carlo configurations for each temperature
        self.nmc = data.shape[0]

        ## list to store the data from the files for sites in the unit cell
        m0, m1, m2 = [], [], []
        theta0, theta1, theta2 = [],[],[]
        phi0, phi1, phi2 = [],[],[]

        for i in range(data.shape[0]):
            data_i = data.loc[i, 'data'].strip().split()
            m0.append(data_i[2])
            m1.append(data_i[3])
            m2.append(data_i[4])
            theta0.append(data_i[5])
            theta1.append(data_i[6])
            theta2.append(data_i[7])
            phi0.append(data_i[8])
            phi1.append(data_i[9])
            phi2.append(data_i[10])

        data_dict = {
            'm1': np.array(m0, dtype='float'),
            'm2': np.array(m1, dtype='float'),
            'm3': np.array(m2, dtype='float'),
            'theta1': np.array(theta0, dtype='float'),
            'theta2': np.array(theta1, dtype='float'),
            'theta3': np.array(theta2, dtype='float'),
            'phi1': np.array(phi0, dtype='float'),
            'phi2': np.array(phi1, dtype='float'),
            'phi3': np.array(phi2, dtype='float')
        }
        return pd.DataFrame(data_dict)


    ### this function will initialize the interacting part of the hamiltonian
    def intpart(self,df,mat,p):
        """
        This function will use the equilibrated monte carlo configuration and generate
        the hamiltonian. The non interacting part of the hamiltonian is already generated
        and the same matrix will be used here
        parameter:
            self(object): instance of the class
            df(dataframe): dataframe that is storing the monte carlo configurations
            mat(matrix): non interacting hamiltonian matrix 
            p(integer): the monte carlo configuration that we want to use
        return:
            None
        """

        ## mx = m0 * cos(phi) * sin(theta)
        ## my = m0 * sin(phi) * sin(theta)
        ## mz = m0 * cos(theta)

        m_x = lambda x,y,z : x*np.cos(y) * np.sin(z)
        m_y = lambda x,y,z : x*np.sin(y) * np.sin(z)
        m_z = lambda x,y : x*np.cos(y) 
        
        for j in range(self.n_cell):
            s0 = j*self.s_cell; s1 = s0 + 1; s2 = s0 +2
            s0n = s0 + (self.s_cell * self.n_cell); s1n = s0n + 1 ; s2n = s0n + 2
            
            ## slice corresponding to a sile
            data_loc = df.iloc[(p*self.n_cell+j),:]
            ### initialize the mx_i,my_i,mz_i for the hamiltonian
            m_0 = data_loc[0] ; m_1 = data_loc[1];m_2 = data_loc[2]
            theta_0 = data_loc[3];theta_1 = data_loc[4];theta_2 = data_loc[5]
            phi_0 = data_loc[6];phi_1 = data_loc[7];phi_2 = data_loc[8]

            m_x0 = m_x(m_0,phi_0,theta_0);m_x1 = m_x(m_1,phi_1,theta_1);m_x2 = m_x(m_2,phi_2,theta_2)
            m_y0 = m_y(m_0, phi_0, theta_0);m_y1 = m_y(m_1, phi_1, theta_1);m_y2 = m_y(m_2,phi_2,theta_2)  
            m_z0 = m_z(m_0,theta_0) ; m_z1 = m_z(m_1,theta_1) ; m_z2 = m_z(m_2,theta_2)
            
            ### diagonal part of the hamiltonian for the up spin
            mat[s0,s0] = -0.5*self.u_int*m_z0 - (self.mu-1.5*self.u_int)
            mat[s1, s1] = -0.5*self.u_int*m_z1 - (self.mu-1.5*self.u_int)
            mat[s2, s2] = -0.5*self.u_int*m_z2 - (self.mu-1.5*self.u_int)

            ### diagonal part of the hamiltonian for the down spin
            mat[s0n, s0n] = 0.5*self.u_int*m_z0  - (self.mu-1.5*self.u_int)
            mat[s1n,s1n] = 0.5*self.u_int*m_z1 - (self.mu-1.5*self.u_int)
            mat[s2n,s2n] = 0.5*self.u_int*m_z2 - (self.mu-1.5*self.u_int)

            ### off diagonal part of the hamiltonian between up,dn spin
            mat[s0,s0n] = -0.5*self.u_int*(m_x0-1j*m_y0)
            mat[s1,s1n] = -0.5*self.u_int*(m_x1-1j*m_y1)
            mat[s2, s2n] = -0.5*self.u_int*(m_x2-1j*m_y2)

            mat[s0n,s0] = -0.5*self.u_int*(m_x0+1j*m_y0)
            mat[s1n,s1] = -0.5*self.u_int*(m_x1+1j*m_y1)
            mat[s2n, s2] = -0.5*self.u_int*(m_x2+1j*m_y2)

        return mat



### this class contains the method and attribute to calculate the observables
### One of the observable we will calculate is the total density
class Observable(Lattice):
    """
    Class initializer.
    """
    def __init__(self,temp,egval,egvec):
        """
        This function will initialize the class or create an object of the type Observable
        parameters:
            temp(float): temperature for which we want to perform the calculation
            egval(array): array containing all the eigenvalues
            egvec(2d array): array containing all the eigenvectors of the system
        """
        self.temp =  temp 
        self.egval = egval 
        self.egvec = egvec


    ## calculate the exponential of the product of beta and E
    def fermiFunc(self):
        """
        This function will calcualate the fermi function for a given temperature
        parameters:
            self(object): Instance of the class
        return:
            ff(array): array storing fermi function for all the energy values
        """
        betaE = np.exp(-(1./self.temp)*self.egval)
        
        self.be = betaE

    
    def feCe(self,j):
        """
        This function calculates the product of fermi function and the wavefunction 
        overlap f_{e}|C|^2
        parameter:
            self(object): instance of the class
            j(int): energy for which we want to calculate the fermi function and wavefn
                    overlap
        return:
            None
        """
        psi  = self.egvec[:,j]
        psi2up = np.conjugate(psi[:self.n_sites])*psi[:self.n_sites]
        psi2dn = np.conjugate(psi[self.n_sites:])*psi[self.n_sites:]
        self.nup = self.be[j]/(1+self.be[j])*psi2up
        self.ndown = self.be[j+self.n_sites]/(1+self.be[j+self.n_sites])*psi2dn
        self.n = self.nup.sum()  + self.ndown.sum()






if __name__=="__main__":

    obs.fermiFunc()
    obs.feCe(1)

    obs.nup.sum()




    obs.n_sites=432



    ### creating an instance of the hamiltonian class
    """
    l_cell = int(12) # number of unit cell in x
    n_cell = int(l_cell**2) # total number of unit cells
    s_cell  = int(3) # number of sites in the unit cell
    n_sites = int(s_cell*n_cell) # total number of site (for one spin type)
    hopping = 1.0   # hopping strength
    dim_h = int(2*n_sites)  # dimensionality of the matrix
    mu = 1.98  ## chemical potential
    u_int = 1.0 ## interaction strength

    ham2 = Hamiltonian(l_cell,n_cell,s_cell,hopping,n_sites,dim_h,mu,u_int)
    ham2.neighbour()

    ham2.haminit()

    evl,evcs = ham2.diag()
    """


    # In[30]:


    l_cell = int(12)  # number of unit cell in x
    n_cell = int(l_cell**2)  # total number of unit cells
    s_cell = int(3)  # number of sites in the unit cell
    n_sites = int(s_cell*n_cell)  # total number of site (for one spin type)
    hopping = 1.0   # hopping strength
    dim_h = int(2*n_sites)  # dimensionality of the matrix
    mu = 13.76  # chemical potential
    u_int = 9.0  # interaction strength




    ## file path were the files are stored
    file_path = './'
    ## pattern we want to search in the file name
    pattern = f'mc*L{l_cell}*Uint{u_int}*.dat'
    #print(pattern)


    ## create instance of the class
    files = Interaction(l_cell, n_cell, s_cell, hopping, n_sites, dim_h, mu, u_int,file_path,pattern)

    ## get names of all the files in the location with the given pattern
    files.getFiles()


    #print(files.data_dict)
    temp = 0.02
    datadf = files.datadf(temp)
    datadf['temp'] = temp

    ### setting up the non interacting hamiltonian
    files.neighbour() ## setting up the neighbour table
    files.haminit() ## initialize the hamiltonian


    ### interacting part of the hamiltonian
    ans = []
    #temp_df = pd.DataFrame()
    #for p in range(1):
    #    mat = files.intpart(datadf, files.ham,p)
    #    evls,evcs = lng.eigh(mat)
        #temp_df = pd.concat([temp_df,pd.Series(evls)],axis=0)
    #    temp_df = pd.Series(evls)
    #    ei_temp = temp_df.describe()
    #    ans.append(ei_temp[ei_temp.index == '50%'].values)



    # In[31]:


    datadf.iloc[:,:-1].plot(kind='kde',subplots=True,layout=(4,3),figsize=(16, 8), sharex=False)


    # In[22]:


    f=plt.figure(figsize=(6,6))
    plt.hist(ans)
    plt.xlabel('E')
    plt.ylabel('H(E)')
    plt.title('L = 12 , T = 0.30 , U = 9')
    #plt.savefig('hist_esymmetric_L12_U9_Temp0p30.eps',dpi = 400,bbox_inches='tight')
    plt.show()


    # In[222]:


    plt.plot(evls,ff,label=0.01)
    #plt.plot(evls2[0::2],ff2[0::2],label='up')
    plt.plot(evls2[1::2],ff2[1::2],label='down')
    plt.axvline(x=1.97,color='')
    plt.legend()

    plt.show()


    # ### Density of states calculations
    # * In this part we will calculate the density of states.
    # * This can be obtained by drawing the histogram of eigenvalues
    # * Steps:    
    #     * Write the hamiltonian using the optimized monte carlo configurations.
    #     * Once hamiltonian is written we can diagonalize the hamiltonian
    #     * One can plot the histogram of all the eigenvalues for different configurations.
    # 
    # 

    # In[30]:


    l_cell = int(12)  # number of unit cell in x
    n_cell = int(l_cell**2)  # total number of unit cells
    s_cell = int(3)  # number of sites in the unit cell
    n_sites = int(s_cell*n_cell)  # total number of site (for one spin type)
    hopping = 1.0   # hopping strength
    dim_h = int(2*n_sites)  # dimensionality of the matrix
    mu = 1.97  # chemical potential
    u_int = 1.0  # interaction strength




    ## file path were the files are stored
    file_path = './'
    ## pattern we want to search in the file name
    pattern = f'mc*L{l_cell}*Uint{u_int}*.dat'
    #print(pattern)


    ## create instance of the class
    files = interaction(l_cell, n_cell, s_cell, hopping, n_sites, dim_h, mu, u_int,file_path,pattern)

    ## get names of all the files in the location with the given pattern
    files.getFiles()


    files.data_dict
    temp = 0.01
    datadf = files.datadf(temp)
    datadf['temp'] = temp

    ### setting up the non interacting hamiltonian
    files.neighbour() ## setting up the neighbour table
    files.haminit() ## initialize the hamiltonian


    ### interacting part of the hamiltonian
    ans = []
    temp_df_U1 = pd.DataFrame()
    for p in range(200):
        mat = files.intpart(datadf, files.ham,p)
        evls,evcs = lng.eigh(mat)
        temp_df_U1 = pd.concat([temp_df_U1,pd.Series(evls)],axis=0)

        
    # In[46]:


    title = f'Density of states L = {l_cell}, T = {temp}, U = {u_int}'
    temp_df_U1.columns = ['ev']
    temp_df_U1.plot(kind='kde',xlabel='E',ylabel='P(E)'     ,title=title)
    #plt.savefig(f'DOS_L{l_cell}_U{u_int}_Temp{temp}.eps',dpi=400,bbox_inches='tight',legend='DOS')
    plt.show()


    # In[38]:


    df =  pd.concat([temp_df_U1,temp_df_U9],axis=1)
    df.columns = ['U1.0','U9.0']


    # In[45]:


    df.plot(kind='kde',legend=['U = 1.0', 'U = 9.0']     ,xlabel='E',ylabel='P(E)',title='Density of States, L = 12, T = 0.01')
    plt.savefig('combinedDos.eps',dpi=400,bbox_inches='tight')
    plt.show()





# %%
