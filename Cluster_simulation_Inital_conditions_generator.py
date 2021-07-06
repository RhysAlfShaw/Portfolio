import numpy as np
import matplotlib.pyplot as plt

"""
Constants and Unit conversions
"""

M_sol = 1.99E30
G = 6.67E-11
AU = 1.49E11



class gen_cluster:

    """
    This class will create a cluster, where options exist for different star distributions
    different IMFs and for other parameters such as N, radius of clusters sphere,
    and the virial ratio q.
    """

    def __init__(self,x):
        self.Rmax = x[0]
        self.N = x[1]
        self.q = x[2]
        self.a = x[0]
        self.b = x[0]
        self.c = x[0]
        self.sigma = x[3]

    

    def Gen_positions(self,distribution,D=1):
        X = np.zeros(self.N)
        Y = np.zeros(self.N)
        Z = np.zeros(self.N)
        if distribution == 'Uniform':
            i = 0
            while i < self.N:
                x = np.random.uniform(-self.Rmax,self.Rmax)
                y = np.random.uniform(-self.Rmax,self.Rmax)
                z = np.random.uniform(-self.Rmax,self.Rmax)                
                if np.sqrt(x**2+y**2+z**2) < self.Rmax :
                    X[i] = x    
                    Y[i] = y
                    Z[i] = z
                    i = i + 1
            self.X = X
            self.Y = Y
            self.Z = Z

        """
        Fractal algorithm goodwin & whitworth (2003)
        creates fractal distribution and also creates associated velocities.
        """

        if distribution == 'Fractal':
    
            gen_number = 0
            n = 0
            while n < self.N:
                if n == 0:
                    parent = np.zeros(1)
                    vel_p = np.array([0,0,0])
                    Mass = np.array([0])
                for i in range(0,len(parent)):
                    vp = vel_p[i]
                    x = self.gen_pos(parent[i],gen_number)
                    vel = self.gen_vel(vp,x)
                    if i == 0 and n == 0:
                        par, vel_p, Mass_new = self.survial(x,vel,D)
                    else:
                        new, vel_n, Mass_new = self.survial(x,vel,D)
                        if len(new) == 0:
                            break
                        else:
                            par = np.vstack((par,new))
                            vel_p = np.vstack((vel_p,vel_n))
                            Mass = np.hstack((Mass,Mass_new))
                parent = par
            
                n = len(parent)
                gen_number += 1
            self.Mass = Mass
            self.X = parent[:,0]
            self.Y = parent[:,1]
            self.Z = parent[:,2]
            print('Number of Generations :',gen_number)
            print('Number of final stars :',n)
            print('Maximum Radius of clsuter:',self.Rmax)
            
            PE_tot = 0
            for k in range(1,len(self.Mass)):
                for j in range(0,k-1):
                    R_res = np.sqrt((parent[k,0]-parent[j,0])**2 + (parent[k,1]-parent[j,1])**2 + (parent[k,2]-parent[j,2])**2)
                    PE_tot += (G * self.Mass[j]*(M_sol**2)*self.Mass[k])/(R_res*AU) 
        
            print('Potential Energy (J) :',PE_tot)
            Ek_tot = 0
            print('Median Mass (Sol Mass)', np.median(self.Mass))
            for i in range(0,len(self.Mass)):
                Ek_tot += 0.5 * self.Mass[i] *M_sol* (np.sqrt(vel_p[i,0]**2 + vel_p[i,1]**2 + vel_p[i,2]**2))**2 
            print('Total kineitc energy (with out added factor): ',Ek_tot)
            a = np.sqrt((self.q*PE_tot)/(Ek_tot))
            print('Velocity factor: ',a)
            self.V = vel_p*a
            print('mean v (m/s):',np.mean(self.V))
            print('median v (m/s):',np.median(self.V))
            print('std v (m/s):',np.std(self.V))
            print('min v (m/s):',np.min(self.V))
            print('max v (m/s):',np.max(self.V))

    

    """
    FRACTAL ALGORITHM FUNCTIONS!
    """


    def Genlen(self,n):
        return self.a*(1/2)**n , self.b*(1/2)**n, self.c*(1/2)**n


    def gen_pos(self,parent,n):
    
        a,b,c = self.Genlen(n)
        one = np.array([[-(1/2)*b , -(1/2)*a , -(1/2)*c],
                        [-(1/2)*b , -(1/2)*a ,  (1/2)*c],
                        [-(1/2)*b ,  (1/2)*a , -(1/2)*c],
                        [ (1/2)*b , -(1/2)*a , -(1/2)*c],
                        [ (1/2)*b ,  (1/2)*a ,  (1/2)*c],
                        [ (1/2)*b ,  (1/2)*a , -(1/2)*c],
                        [ (1/2)*b , -(1/2)*a ,  (1/2)*c],
                        [-(1/2)*b ,  (1/2)*a ,  (1/2)*c]])
        noise = np.random.normal(0,self.sigma, size=(8,3))*self.a  #adds gaussian noise
        new = np.subtract(one, parent) + noise
        return new

    def prob_of_mature(self,D):
        return 2**(D-3)
    

    """
    Gen_vel generates velocities of children around a parent where they are themselves in virial equilibrium.
    """


    def gen_vel(self,vp,x):
        vel = np.random.uniform(-1,1,size=(8,3))
        Mass = self.cust_dis(len(vel),0,20,self.Kroupa_IMF)
        X, Y, Z = x[:,0], x[:,1], x[:,2]
        PE_tot = 0
        for k in range(1,len(X)): 
            for j in range(0,k-1):
                R_res = np.sqrt((X[k]-X[j])**2 + (Y[k]-Y[j])**2 + (Z[k]-Z[j])**2)
                PE_tot += (G*Mass[k]*M_sol**2*Mass[j])/(R_res*AU)     
        Ek_tot = 0
        for i in range(0,len(X)):
            Ek_tot += 0.5 * Mass[i]*M_sol* (np.sqrt(vel[i,0]**2 + vel[i,1]**2 + vel[i,2]**2))**2
        a = np.sqrt((0.5*PE_tot)/(Ek_tot)) 
        vel = a*vel - vp                            
        return vel, Mass


    """
    SURVIVAL
    Determines which of the children will mature and only includes those that survived.
    Returns array of positions, velocity and mass of the mature children.
    """


    def survial(self,children,vel,D):
        parent = []
        vel_sur = []
        v, mass = vel
        Mass = []
        for i in range(0,len(children)):
            rand = np.random.uniform(0,1)
            prob = self.prob_of_mature(D)
            if rand<= prob:
                parent+=[children[i]]
                vel_sur+=[v[i]] 
                Mass += [mass[i]]
        return np.array(parent), np.array(vel_sur), np.array(Mass)


    """
    IMF FUNCTIONS!

    """

    def Gen_mass(self,IMF):
        if IMF == 'constant':
            Mass = np.zeros(self.N)

            for i in range(0,self.N):
                Mass[i] = 0.2
            self.Mass = Mass
        if IMF == 'KROUPA':
            self.Mass = self.cust_dis(self.N,0,20,self.Kroupa_IMF)
            print('Median Mass: ',np.median(self.Mass))
        else:
            print("INPUT ERROR: Please specify 'constant' or 'KROUPA' IMFs")



    def cust_dis(self,N,x0,x1,imf,nControl=10**6):
        sample = []
        nLoop  = 0
        while len(sample)<N and nLoop<nControl:
            x = np.random.uniform(x0,x1)     
            prop = imf(x)
            assert prop>=0
            if np.random.uniform(0,1) <= prop: #26.67 is max probability 
                sample += [x]
            nLoop+= 1
        return np.array(sample)

    
    def Kroupa_IMF(self,m):
        alpha = 0
        if m > 0 and m < 0.08:
            alpha = 0.3
        if m > 0.08 and m < 0.5:
            alpha = 1.3
        if m > 0.5:
            alpha = 2.3
        return m**(-alpha)


    """
    ENERGY BALANCE AND VELOCITY!
    """


    def Gen_velocities(self):
        self.V = np.random.uniform(-1,1,size=(self.N,3))
        PE_tot = 0
        
        # Calculating Graviational Potential!

        for k in range(1,self.N):     #note since 0 is the begining of the index.
             for j in range(0,k-1):
                R_res = np.sqrt((self.X[k]-self.X[j])**2 + (self.Y[k]-self.Y[j])**2 + (self.Z[k]-self.Z[j])**2)
                PE_tot += (G * self.Mass[j]*M_sol**2*self.Mass[k])/(R_res*AU) 
        print('Potential Energy',PE_tot)


        # Calculating Kinetic Energy!
        Ek_tot = 0
        for i in range(0,self.N):
            Ek_tot += 0.5 * self.Mass[i] *M_sol* (np.linalg.norm(self.V))**2
        print('Total kineitc energy (with out added factor): ',Ek_tot)
        a = np.sqrt((self.q*PE_tot)/(Ek_tot))
        print('Velocity factor: ',a)
        self.V = a*self.V



    """
    GRAPHING ;)
    """

    def graph(self,dimentions):
        if dimentions == '3D':
            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            ax.scatter(self.X,self.Y,self.Z,s=1)
            ax.set_box_aspect([1,1,1])
            plt.show()
        if dimentions == '2D':
            plt.scatter(self.X,self.Y,s=1)
            plt.show()




class get_planet:
    """
    This class will create a planetary systems based on input parameters on a clsuters system.
    Takes the inputs of theta (to pick the starting position of the orbit), M_star, M_planet, Phi and eccentricity and semi major axis a.
    Phi should all for the inclination of the orbit to be changed. 
    """
    def __init__(self,x): #x must be a list of each of these parameters! 
        self.theta = x[0] #where on the orbit inradians
        self.M_star = x[1]
        self.M_planet = x[2]
        self.phi = x[3]
        self.a = x[4]
        self.e = x[5]
        self.Np = len(x[2])
        self.theta_1 = x[6]
        self.phi = x[7]
        self.sig = x[8]

    def calculate_params(self):
        self.c = self.a*self.e
        self.b = np.sqrt((self.a**2)- (self.c**2))
        self.R_vec = np.zeros((self.Np,3))
        self.V_vec = np.zeros((self.Np,3))
        for i in range(0,self.Np):
            self.R_vec[i,0] = self.a[i]*np.cos(self.theta[i])-self.c[i]
            self.R_vec[i,1] = self.b[i]*np.sin(self.theta[i])
            self.R_vec[i,2] = 0
            self.R_vec[i,:] = self.rotation(self.R_vec[i,:],self.theta_1,self.phi,self.sig) 
            
            Vmax = np.sqrt(G*(self.M_star+self.M_planet[i])*M_sol * ((2)/(np.linalg.norm(self.R_vec[i])*AU) - (1)/(self.a[i]*AU)))
            self.V_vec[i,0] = Vmax*np.sin(self.theta[i])
            self.V_vec[i,1] = Vmax*np.cos(self.theta[i])
            self.V_vec[i,2] = 0
            self.V_vec[i,:] = self.rotation(self.V_vec[i,:],self.theta_1,self.phi,self.sig)


            
    def ellipse(self,t,a,b,c):
        x = -c + a*np.cos(t)
        y = b*np.sin(t)
        z = np.zeros(len(t))
        X = np.array([x,y,z])
        return  X.T



    def rotation(self,X,theta_1,phi,sig):
        M_rot = np.array([[np.cos(theta_1)*np.cos(phi),np.cos(theta_1)*np.sin(phi), - np.sin(theta_1)],
                [- np.cos(sig)*np.sin(phi) + np.sin(sig)*np.sin(theta_1)*np.cos(phi),np.cos(sig)*np.cos(phi) + np.sin(sig)*np.sin(theta_1)*np.sin(phi),np.sin(sig)*np.cos(theta_1)],
                [np.sin(sig)*np.sin(theta_1) + np.cos(sig)*np.sin(theta_1)*np.cos(phi),-np.sin(sig)*np.cos(theta_1) + np.cos(sig)*np.sin(theta_1)*np.sin(phi),np.sin(sig)*np.cos(theta_1)]])
        return  np.matmul(M_rot,X)



    def graph(self,dimension):

        t = np.linspace(0,2*np.pi, 200)
        if dimension == '2D':
            plt.figure(figsize=[5*np.max(self.a),5*np.max(self.b)]) #figsize changes to appropriate deimensions.
            for i in range(0,self.Np):
                X = self.ellipse(t,self.a[i],self.b[i],self.c[i])
                for j in range(0,200):
                    X[j] = self.rotation(X[j],self.theta_1,self.phi,self.sig)
                plt.plot(X[:,0],X[:,1],linestyle='--')
                plt.scatter(self.R_vec[i,0],self.R_vec[i,1],label='Planet'+str(i+1)+'M ='+str(self.M_planet[i]))
            plt.scatter(0,0,label='Star')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.legend()
            plt.show()



        if dimension == '3D':


            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            ax.scatter(0,0,0,label='star',s=40)
            for i in range(0,self.Np):
                X = self.ellipse(t,self.a[i],self.b[i],self.c[i])
                for j in range(0,200):
                    X[j] = self.rotation(X[j],self.theta_1,self.phi,self.sig)
                ax.plot(X[:,0],X[:,1],X[:,2],linestyle='--')
                ax.scatter(self.R_vec[i,0],self.R_vec[i,1],self.R_vec[i,2],label='Planet'+str(i+1)+'M ='+str(self.M_planet[i]))
            ax.set_box_aspect([1*np.max(self.a),1*np.max(self.b),1]) #figsize changes to appropriate dimentions
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.legend()
            plt.show()


        




class measure:

    """
    allows for measurements of a system. Giving only R_vec and V_vec. This will be important
    when measuring the properties of otbits over time.
    
    """

    def __init__(self,data):
        self.M_star = data[0]
        self.M_planet = data[1]
        self.R_vec = data[2]
        self.V_vec = data[3]
    

    def calculate_a(self):
        mu = self.M_star*M_sol*M_sol*self.M_planet/(self.M_planet*M_sol+self.M_star*M_sol)
        Eb = abs(0.5*mu*np.linalg.norm(self.V_vec)**2  - (G*self.M_star*M_sol*self.M_planet*M_sol)/(np.linalg.norm(self.R_vec)*AU))
        a_calc = (G * self.M_star*M_sol * self.M_planet*M_sol) / (2*Eb*AU)
        return a_calc #units of AU
    

    def calculate_e(self,a):
        e = np.sqrt((1-(np.linalg.norm(self.R_vec))/(abs(a)) )**2 + (np.dot(self.R_vec,self.V_vec)**2)/(abs(a)*G*(self.M_planet+self.M_star)*M_sol))
        return e

    
    def graph(self,dimension,a,e):
        c = a*e #from the definition of eccentricity.
        b = np.sqrt((a**2) - (c**2))
        t = np.linspace(0,2*np.pi, 200)
        x = -c + a*np.cos(t)
        y = b*np.sin(t)

        if dimension == '2D':

            plt.figure(figsize=[5*a,5*b]) #figsize changes to appropriate deimensions.
            plt.plot(x,y,linestyle='--')
            plt.scatter(0,0,label='Star')
            plt.scatter(self.R_vec[0],self.R_vec[1],label='Planet')
            plt.legend()
            plt.show()

        if dimension == '3D':

            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            ax.scatter(0,0,0,label='star',s=40)
            ax.plot(x,y,0,linestyle='--')
            ax.scatter(self.R_vec[0],self.R_vec[1],self.R_vec[2],label='Planet')
            ax.set_box_aspect([1*a,1*b,1]) #figsize changes to appropriate dimentions
            plt.legend()
            plt.show()