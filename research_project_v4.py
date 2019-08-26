import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline
import copy
import time


def grid_nodes(seed_num = (3,3),max_x = 1,max_y = 1):
    node_ls = []

    for i in range(seed_num[0]):
        for j in range(seed_num[1]):
            node_ls.append([float(i),float(j)])
    node_ls = np.array(node_ls)


    node_ls[:,0] = node_ls[:,0]*max_x/(seed_num[0]-1)
    node_ls[:,1] = node_ls[:,1]*max_y/(seed_num[1]-1)
    return node_ls


def grid_ridges(seed_num = (3,3)):
    ridges = []
    sq = seed_num[0]*seed_num[1]
    for i in range(sq):
        if i+1 < sq and (i+1)//seed_num[1] == i//seed_num[1]:
            ridges.append([i,i+1])
        if i+seed_num[1]-1 < sq and (i+seed_num[1]-1)//seed_num[1] == (i//seed_num[1])+1:
            ridges.append([i,i+seed_num[1]-1])
        if i+seed_num[1] < sq:
            ridges.append([i,i+seed_num[1]])
        if i+seed_num[1]+1 < sq and (i+seed_num[1]+1)//seed_num[1] == (i//seed_num[1])+1:
            ridges.append([i,i+seed_num[1]+1])
        # i < another node index
    return np.array(ridges,dtype="int")


#--------------------------------------------------------------------------------------------------------------------------------


class adjacent_list():
    def __init__(self, nodes, ridge_nodes,gird_seed,max_x,max_y,total_mass = 1,k = 10,damp = 5):
        
        self.nodes = nodes
        self.ridge_nodes = ridge_nodes.astype("int")
        self.node_num = nodes.shape[0]
        self.node_num_x = grid_seed[0]
        self.node_num_y = grid_seed[1]
        self.max_x = max_x
        self.max_y = max_y
        self.node_ls = []
        for i in range(self.node_num):
            self.node_ls.append(adjacent_node(i,self.nodes[i],damp))
        self.add_neighbor_global()
        self.distance_update_global("init")
        self.fixed_nodes() # fix all nodes on the bottom
        self.external_force() # apply force to nodes
        self.boundary_sum = 0
        self.boundary_nodes() # add label to boundary nodes
        self.stiffness = k
        self.charged_p = [] # record those charged nodes
        self.charged_n = [] # record those charged nodes
        self.charged = [] # record those charged nodes
        
        self.ridge_ls = []
        for i in self.ridge_nodes:
            n1_num = i[0]
            n2_num = i[1]
            self.ridge_ls.append(adjacent_ridge(self.node_ls[n1_num],self.node_ls[n2_num]))
        self.apply_voltage(pos = "bottom",voltage = -20)
        self.apply_voltage(pos = "top_mid1/3",voltage = 200)
        for i in self.ridge_ls:
            i.update_from_nodes(first = True)
        self.update_charge_density()
        self.update_electric_force()
        

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Geometry
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def add_neighbor_global(self,):
        for rn in self.ridge_nodes:
            n1,n2 = rn[0],rn[1]
            self.node_ls[n1].add_neighbor(n2)
            self.node_ls[n2].add_neighbor(n1)
        

    def distance_update_global(self,mode):
        if mode == "init":
            for ridge in self.ridge_nodes:
                i,j = ridge[0],ridge[1]
                pos_i,pos_j = np.array(self.node_ls[i].position),np.array(self.node_ls[j].position)
                dist = pos_j - pos_i
                self.node_ls[i].init_distance_to_neighbor(j,dist)
                self.node_ls[j].init_distance_to_neighbor(i,-dist)
                self.node_ls[i].update_distance_to_neighbor(j,dist)
                self.node_ls[j].update_distance_to_neighbor(i,-dist)

        elif mode == "update":
            for ridge in self.ridge_nodes:
                i,j = ridge[0],ridge[1]
                pos_i,pos_j = np.array(self.node_ls[i].position),np.array(self.node_ls[j].position)
                dist = pos_j - pos_i
                self.node_ls[i].update_distance_to_neighbor(j,dist)
                self.node_ls[j].update_distance_to_neighbor(i,-dist)
        else:
            print("Error! No such mode for distance calculation!")

    def fixed_nodes(self,): # boundary condition
        for i in range(self.node_num):
            if self.node_ls[i].position[1] == 0:
                self.node_ls[i].fixed_or_not = True
    
    def boundary_nodes(self,):
        self.boundary_sum = 0
        for i in range(self.node_num):
            node = self.node_ls[i].position
            if node[0] == 0 or node[1] == 0 or node[0] == 1 or node[1] == 1:
                self.boundary_sum += 1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Mechanics
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def external_force(self,pos = "top", f = np.array([0,0])):
        if pos == "none":
            pass
        elif pos == "top":
            for i in range(self.node_num):
                if self.node_ls[i].position[1] == self.max_y:
                    self.node_ls[i].fou = f
        elif pos == "topmid1/3":
            for i in range(self.node_num):
                if self.node_ls[i].position[1] == self.max_y and abs(self.node_ls[i].position[0]-0.5) < 1/6+0.01:
                    self.node_ls[i].fou = f
        elif pos == "top_stretch":
            f1 = np.array([-1,1])
            f2 = np.array([1,1])
            for i in range(self.node_num):
                if self.node_ls[i].position[1] == self.max_y:
                    if self.node_ls[i].position[0]-0.5 < -0.01:
                        self.node_ls[i].fou = f1
                    elif self.node_ls[i].position[0]-0.5 > 0.01:
                        self.node_ls[i].fou = f2
        else:
            print("Error! Force mode not supported yet")

    def consitutive_eq(self,length_init,length_now):
        epsilon = (length_now-length_init)/length_init # epsilon>0 => direction of direction.
        # direction: from central node to planet nodes
        Fin = np.tan(epsilon*2/np.pi)*self.stiffness*length_init # when epsilon≈0, Fin≈△x*k
        return Fin
    
    def Fin_update_global(self,):
        for ridge in self.ridge_nodes:
            i,j = ridge[0],ridge[1]
            fin = self.consitutive_eq(self.node_ls[i].length_init[j],self.node_ls[i].length_now[j])
            # if self.node_ls[i].length_now[j] != self.node_ls[j].length_now[i]:
                # print("Error! lengths from two nodes are not mutally equal!")
            self.node_ls[i].update_fin(j,fin)
            self.node_ls[j].update_fin(i,fin)
    
    def one_time_step(self,flag,mode="EC"):
        self.Fin_update_global()
        if mode == "EC":
            for i in range(self.node_num):
                self.node_ls[i].move_up_EC()
        if mode == "LF":
            for i in range(self.node_num):
                self.node_ls[i].move_up_LF()
        self.distance_update_global("update")
        if flag == True:
            for i in self.ridge_ls:
                i.update_from_nodes()
            self.update_charge_density()
            self.update_electric_force()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Electrostatics
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def global_charge_clear(self,):
        for i in range(self.node_num):
            self.node_ls[i].charge_density = 0

    def apply_voltage(self,pos = "top",voltage = 1000,mode = "add"):
        if mode == "add":
            if pos == "bottom":
                for t in range(self.node_num_x):
                    i = t*self.node_num_y
                    self.node_ls[i].voltage = voltage
                for t in range(self.node_num_x-1):
                    i = t*self.node_num_y
                    j = (t+1)*self.node_num_y
                    for index,r in enumerate(self.ridge_nodes):
                        if np.equal(r,np.array([i,j])).all():
                            ridge_num = index
                            break
                    self.charged_n.append(ridge_num)
                    self.charged.append(ridge_num)

            elif pos == "top":
                for t in range(self.node_num_x):
                    i = (t+1)*self.node_num_y-1
                    self.node_ls[i].voltage = voltage
                for t in range(self.node_num_x-1):
                    i = (t+1)*self.node_num_y-1
                    j = i+self.node_num_y
                    for index,r in enumerate(self.ridge_nodes):
                        if np.equal(r,np.array([i,j])).all():
                            ridge_num = index
                            break
                    self.charged_p.append(ridge_num)
                    self.charged.append(ridge_num)

            elif pos == "top_mid1/3":
                i_ls = []
                for t in range(self.node_num_x):
                    i = (t+1)*self.node_num_y-1
                    if self.node_ls[i].position[0] >= self.max_x/3 and self.node_ls[i].position[0] <= self.max_x*2/3:
                        self.node_ls[i].voltage = voltage
                        i_ls.append(i)
                for t,i in enumerate(i_ls[:-1]):
                    j = i_ls[t+1]
                    for index,r in enumerate(self.ridge_nodes):
                        if np.equal(r,np.array([i,j])).all():
                            ridge_num = index
                            break
                    self.charged_p.append(ridge_num)
                    self.charged.append(ridge_num)

    def update_charge_density(self,
                permittivity = 3 * 8.85*10e-12): # relativity permittivity = 3
        # solve this linear algebra problem: L\alpha=g
        # reference https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-635-advanced-electromagnetism-spring-2003/lecture-notes/Mar10.pdf
        # print("\ncharged n p:",self.charged_n,self.charged_p)
        n_sqr_n = len(self.charged_n)
        n_sqr_p = len(self.charged_p)
        # print("\nn num, p num:",n_sqr_n,n_sqr_p)
        num = n_sqr_n**2 + n_sqr_p**2
        L = np.zeros((num,num))
        g = np.zeros(num)
        y_ls_n = []
        y_ls_p = []

        for i in range(n_sqr_n**2):
            r = self.ridge_ls[self.charged_n[i//n_sqr_n]]
            if r.position[0] not in y_ls_n:
                y_ls_n.append(r.position[0]) # prepare for the 3D virtual rebuilding
        for i in range(n_sqr_p**2):
            r = self.ridge_ls[self.charged_p[i//n_sqr_p]]
            if r.position[0] not in y_ls_p:
                y_ls_p.append(r.position[0]) # prepare for the 3D virtual rebuilding
        # print("y_ls_n,y_ls_p:",y_ls_n,y_ls_p)

        # voltage
        for i in range(n_sqr_n**2):
            g[i] = self.ridge_ls[self.charged_n[i//n_sqr_n]].voltage
        for i in range(n_sqr_p**2):
            g[i+n_sqr_n**2] = self.ridge_ls[self.charged_p[i//n_sqr_p]].voltage
        print("g:",g)

        # mutual effect Matrix
        for i in range(n_sqr_n**2):
            for j in range(n_sqr_n**2):
                if i != j:
                    ri = self.ridge_ls[self.charged[i//n_sqr_n]]
                    rj = self.ridge_ls[self.charged[j//n_sqr_n]]
                    xi = ri.position[0]
                    yi = y_ls_n[i%n_sqr_n]
                    zi = ri.position[1]
                    xj = rj.position[0]
                    yj = y_ls_n[j%n_sqr_n]
                    zj = rj.position[1]
                    distance = np.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
                    area = rj.length**2
                    permittivity = 3 * 8.85*10e-12 # relativity permittivity = 3
                    L[i,j] = area/(4*np.pi*permittivity*distance)
                else:
                    ri = self.ridge_ls[self.charged[i//n_sqr_n]]
                    L[i,j] = ri.length * np.log(1+np.sqrt(2)) /(np.pi*permittivity)

        for i in range(n_sqr_n**2):
            for j in range(n_sqr_p**2):
                ri = self.ridge_ls[self.charged[i//n_sqr_n]]
                rj = self.ridge_ls[self.charged[n_sqr_n + j//n_sqr_p]]
                xi = ri.position[0]
                yi = y_ls_n[i%n_sqr_n]
                zi = ri.position[1]
                xj = rj.position[0]
                yj = y_ls_n[j%n_sqr_p]
                zj = rj.position[1]
                distance = np.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
                area = rj.length**2
                permittivity = 3 * 8.85*10e-12 # relativity permittivity = 3
                L[i,j+n_sqr_n**2] = area/(4*np.pi*permittivity*distance)
                L[j+n_sqr_n**2,i] = area/(4*np.pi*permittivity*distance)

        for i in range(n_sqr_p**2):
            for j in range(n_sqr_p**2):
                if i != j:
                    ri = self.ridge_ls[self.charged[n_sqr_n + i//n_sqr_p]]
                    rj = self.ridge_ls[self.charged[n_sqr_n + j//n_sqr_p]]
                    xi = ri.position[0]
                    yi = y_ls_p[i%n_sqr_p]
                    zi = ri.position[1]
                    xj = rj.position[0]
                    yj = y_ls_p[j%n_sqr_p]
                    zj = rj.position[1]
                    distance = np.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
                    area = rj.length**2
                    permittivity = 3 * 8.85*10e-12 # relativity permittivity = 3
                    L[i+n_sqr_n,j+n_sqr_n] = area/(4*np.pi*permittivity*distance)
                else:
                    ri = self.ridge_ls[self.charged[n_sqr_n + i//n_sqr_p]]
                    L[i+n_sqr_n**2,j+n_sqr_n**2] = ri.length * np.log(1+np.sqrt(2)) /(np.pi*permittivity)

        # print("L:",L)
        alpha = np.linalg.solve(L,g)
        # print("\nalpha:",alpha)
        print("electric charge balance?",alpha[:n_sqr_n**2].sum(),alpha[n_sqr_n**2:].sum())
        
        
        for i in range(n_sqr_n):
            self.ridge_ls[self.charged[i]].charge_density = alpha[i*n_sqr_n:(i+1)*n_sqr_n].mean()
        for i in range(n_sqr_p):
            self.ridge_ls[self.charged[n_sqr_n + i]].charge_density = alpha[n_sqr_n**2 + i*n_sqr_p : n_sqr_n**2 + (i+1)*n_sqr_p].mean()

    def update_electric_force(self,Ke = 9*10e9):
        for i in self.charged:
            r = self.ridge_ls[i]
            r.charge = r.length * r.charge_density
        for i1 in self.charged:
            r1 = self.ridge_ls[i1]
            r1.force = 0
            for i2 in self.charged:
                if i2 != i1:
                    r2 = self.ridge_ls[i2]
                    x1 = r1.position[0]
                    x2 = r2.position[0]
                    y1 = r1.position[1]
                    y2 = r2.position[1]
                    direction = r1.position-r2.position
                    direction = direction/np.linalg.norm(direction)
                    if np.isnan(direction).any():
                        print("x1 y1 x2 y2:",x1,y1,x2,y2)
                    r1.force += direction * Ke*r1.charge*r2.charge/((x1-x2)**2 + (y1-y2)**2)
            print("r1.force",r1.force)
        for i in self.charged:
            r = self.ridge_ls[i]
            n1 = r.n1
            n2 = r.n2
            n1.fele = 0
            n2.fele = 0
        for i in self.charged:
            r = self.ridge_ls[i]
            n1 = r.n1
            n2 = r.n2
            n1.fele += 0.5*r.force
            n2.fele += 0.5*r.force
        print("--------------------------------------------------")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Visualization
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

 
    def render_grid(self,mode = "boundary"):
        if mode == "boundary":
            # seed_num = int(np.sqrt(self.node_num))
            plt.figure(figsize=(4,4))
            plt.xlim(-0.5,1.5)
            plt.ylim(-0.5,1.5)
            # left
            for i in range(self.node_num_y-1):
                j = (i+1)
                xi = self.node_ls[i].position[0]
                xj = self.node_ls[j].position[0]
                yi = self.node_ls[i].position[1]
                yj = self.node_ls[j].position[1]
                plt.plot([xi,xj],[yi,yj],c="r")
            # right
            for i in range(self.node_num-self.node_num_y,self.node_num-1):
                j = (i+1)
                xi = self.node_ls[i].position[0]
                xj = self.node_ls[j].position[0]
                yi = self.node_ls[i].position[1]
                yj = self.node_ls[j].position[1]
                plt.plot([xi,xj],[yi,yj],c="r")
            # bottom
            for t in range(self.node_num_x-1):
                i = t*self.node_num_y
                j = (t+1)*self.node_num_y
                xi = self.node_ls[i].position[0]
                xj = self.node_ls[j].position[0]
                yi = self.node_ls[i].position[1]
                yj = self.node_ls[j].position[1]
                plt.plot([xi,xj],[yi,yj],c="r")
            # top
            for t in range(self.node_num_x-1):
                i = (t+1)*self.node_num_y-1
                j = i+self.node_num_y
                xi = self.node_ls[i].position[0]
                xj = self.node_ls[j].position[0]
                yi = self.node_ls[i].position[1]
                yj = self.node_ls[j].position[1]
                plt.plot([xi,xj],[yi,yj],c="r")
            plt.show()
        elif mode == "all":
            plt.figure(figsize=(4,4))
            plt.xlim(-0.5,1.5)
            plt.ylim(-0.5,1.5)
            for i in range(self.node_num):
                plt.scatter(self.node_ls[i].position[0],self.node_ls[i].position[1],c="k")
            plt.show()
        elif mode == "all+boundary" or mode == "boundary+all":
            plt.figure(figsize=(4,4))
            plt.xlim(-0.5,2)
            plt.ylim(-0.5,.5)
            # left
            for i in range(self.node_num_y-1):
                j = (i+1)
                xi = self.node_ls[i].position[0]
                xj = self.node_ls[j].position[0]
                yi = self.node_ls[i].position[1]
                yj = self.node_ls[j].position[1]
                plt.plot([xi,xj],[yi,yj],c="r")
            # right
            for i in range(self.node_num-self.node_num_y,self.node_num-1):
                j = (i+1)
                xi = self.node_ls[i].position[0]
                xj = self.node_ls[j].position[0]
                yi = self.node_ls[i].position[1]
                yj = self.node_ls[j].position[1]
                plt.plot([xi,xj],[yi,yj],c="r")
            # bottom
            for t in range(self.node_num_x-1):
                i = t*self.node_num_y
                j = (t+1)*self.node_num_y
                xi = self.node_ls[i].position[0]
                xj = self.node_ls[j].position[0]
                yi = self.node_ls[i].position[1]
                yj = self.node_ls[j].position[1]
                plt.plot([xi,xj],[yi,yj],c="r")
            # top
            for t in range(self.node_num_x-1):
                i = (t+1)*self.node_num_y-1
                j = i+self.node_num_y
                xi = self.node_ls[i].position[0]
                xj = self.node_ls[j].position[0]
                yi = self.node_ls[i].position[1]
                yj = self.node_ls[j].position[1]
                plt.plot([xi,xj],[yi,yj],c="r")
            for i in range(self.node_num):
                plt.scatter(self.node_ls[i].position[0],self.node_ls[i].position[1],c="k",s=3)
            plt.show()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Define Ridges
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class adjacent_ridge():
    def __init__(self,node1,node2):
        if node1.seq_num > node2.seq_num:
            self.n1 = node2
            self.n2 = node1
        else:
            self.n2 = node2
            self.n1 = node1
        self.voltage = None
        self.charge_density = 0.5*(node1.charge_density + node2.charge_density)
        self.length = None
        self.position = None
        self.update_from_nodes()
        self.charge = None
        self.electric_force = None

    # should update from nodes every iteration since length may get changed
    def update_from_nodes(self,first = False): # geometry and voltage information
        if first == True:
            if self.n1.voltage == None or self.n2.voltage == None:
                self.voltage = None
            elif self.n1.voltage != self.n2.voltage:
                self.voltage = None # prevent when 2 y layers, some ridges has a fake "0" voltage
            else:
                self.voltage = 0.5*(self.n1.voltage + self.n2.voltage)
        # print("n1 n2 voltage:",self.n1.voltage,self.n2.voltage)
        # self.charge_density = 0.5*(self.n1.charge_density + self.n2.charge_density)
        self.length = np.sqrt((self.n1.position[0] - self.n2.position[0])**2 + (self.n1.position[1] - self.n2.position[1])**2)
        self.position = np.array([0.5*(self.n1.position[0] + self.n2.position[0]),0.5*(self.n1.position[1] + self.n2.position[1])])
        # print("position:",self.position)
        # print("ridge voltage:",self.voltage)

    def update_to_nodes(self,): # charge density information
        self.n1.charge_density += self.charge_density*0.5 
        self.n2.charge_density += self.charge_density*0.5

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Define Nodes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class adjacent_node():
    def __init__(self,seq_num,position,damp = 1):
        self.seq_num = seq_num
        self.position = position # should be np array
        self.velocity = np.array([0.0,0.0])
        self.acceleration = np.array([0.0,0.0])
        self.fixed_or_not = False
        self.damp = damp
        self.neighbor = []
        self.fin = {}
        self.length_init = {}
        self.length_now = {}
        self.direction = {}
        self.fou = np.array([0,0])
        self.fele = np.array([0,0])
        self.f_sum = np.array([0,0])
        self.counter = 0
        self.voltage = None
        self.charge_density = 0


    def add_neighbor(self,n2):
        self.neighbor.append(n2)

    def init_distance_to_neighbor(self,n2,distance_n1_to_n2):
        self.length_init[n2] = np.linalg.norm(distance_n1_to_n2)

    def update_distance_to_neighbor(self,n2,distance_n1_to_n2): # only need to traverse every edge in graph
        self.length_now[n2] = np.linalg.norm(distance_n1_to_n2)
        self.direction[n2] = distance_n1_to_n2/self.length_now[n2]

    def update_fin(self,n2,fin):
        self.fin[n2] = fin*self.direction[n2]

    def move_up_EC(self,delta_t = 0.1): # dynamic equation,mode = "Euler-Cromer",
        fin_total = np.array([0.0,0.0])
        if self.fixed_or_not == False:
            # Euler, Euler-Cromer, leap frog, Verlet
            for j in self.fin.keys():
                fin_total += self.fin[j]
            self.f_sum = fin_total+self.fou+self.fele
            self.acceleration = (self.f_sum-self.velocity*self.damp)
            self.velocity += self.acceleration * delta_t
            self.position += self.velocity * delta_t
        if self.fixed_or_not == True:
            self.acceleration = np.array([0.0,0.0])
            self.velocity = np.array([0.0,0.0])

    def move_up_LF(self,delta_t = 0.04):
        fin_total = np.array([0.0,0.0])
        if self.fixed_or_not == False:
            if self.counter == 0:
                for j in self.fin.keys():
                    fin_total += self.fin[j]
                self.f_sum = fin_total+self.fou+self.fele
                self.acceleration = (self.f_sum-self.velocity*self.damp)
                self.velocity += self.acceleration * delta_t * 0.5
                self.position += self.velocity * delta_t
                self.counter += 1
            if self.counter != 0:
                for j in self.fin.keys():
                    fin_total += self.fin[j]
                self.f_sum = fin_total+self.fou+self.fele
                self.acceleration = (self.f_sum-self.velocity*self.damp)
                self.velocity += self.acceleration * delta_t
                self.position += self.velocity * delta_t
        if self.fixed_or_not == True:
            self.acceleration = np.array([0.0,0.0])
            self.velocity = np.array([0.0,0.0])


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Implementation
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


max_x = 1.6
max_y = 0.1

grid_seed = (33,3)
nodes2 = grid_nodes(grid_seed,max_x=max_x,max_y=max_y)
ridges2 = grid_ridges(grid_seed)
# print("grid nodes:",nodes2)
# print("grid ridges:",ridges2)


my_adjacent_list = adjacent_list(nodes2,ridges2,grid_seed,max_x=max_x,max_y=max_y, \
                                 total_mass = max_x*max_y,k=20,damp=10)


y = [] # postition list for ploting

ite = 400
flag = False
for i in range(ite):
    if i % int(ite/5) == 0:
        flag = True
    else:
        flag = False
    my_adjacent_list.one_time_step(flag)
    y.append(my_adjacent_list.node_ls[1].position[1])
my_adjacent_list.render_grid("all+boundary")
# my_adjacent_list.render_grid("boundary")

# for c in my_adjacent_list.charged:
#     r = my_adjacent_list.ridge_ls[c]
#     print("\nridge num:",c,":",r.position)
    # print("nodes position:",r.position)
    # print("voltage:",r.voltage)

# for r in my_adjacent_list.ridge_ls:
#     print("\nridge",r.position,":",r.charge_density)

# for n in my_adjacent_list.node_ls:
#     print("\nnode",n.position,":",n.fele)

plt.plot(y)
plt.show()

