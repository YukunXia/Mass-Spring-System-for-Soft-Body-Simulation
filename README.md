# Mass-Spring-System-for-Soft-Body-Simulation

### 1. Project Background

From continuum assumption, elasticity theory has been perfectly developed for a long time. According to the results from elasticity theory, computational scientists implemented discretization and weak form of PDE to build the Finite Element Method (FEM). FEM has already been successfully and broadly used in structure analysis as well as multiphysical simulation. However, FEM essentially needs to solve a matrix problem ”Ax=b”, while here matrix A is usually pretty large, and thus it’s a computationally intensive task to solve vector x. For a 3D structure, if there are n nodes in each edge, then there are n 3 nodes totally, and solving A −1 b has the complexity O(n 9 ). Another latent problem for FEM is the diﬃculties to keep programs converging when local deformation is large, and nonlinear eﬀect dominates there.

Switching from area or voxel view of material, we can also regard them as a combination of particles. This perspective belongs to mesh-free method. The ﬁrst mesh-free method is Smooth Particle Method. Liquids are divided into particles, with particles moving and interacting. If we know the position of every particle, we can then calculate the surface of liquid by multiplying particle positions with the kernel function. A smooth kernel can always lead to smooth surface.

For solid material, the relationship between particles should have much more limits. One way is to add springs to connect neighborhood particles. Since this method only count local interaction, it’s quite computationally cheap, and can work for real time simulation [1] [2]. In this project, I would ﬁrstly introduce how to build such a spring-mass system, and then demonstrate some simple cases to verify the feasibility, and ﬁnally, electrostatic forces are incorporated for Dielectric Elastomer Actuator study.

### 2. Project Analysis

A heuristic way to solve deformation problem is to regard a block of solid as balls connected by springs, and the linkage only exists between adjacent balls. Discretizing solids into balls, and setting boundary condition are similar to FEM. The next step is to apply force to certain balls, so these balls will have acceleration and start to move, like the spreading of strain wave. For a dynamical system like this, the motion will never end. To stabilize the system, artiﬁcial damping is added to every nodes, so every ball suﬀers a negative force in its direction of moving. Finally, there will only be one stable position of every node, and the deformation has been simulated.

If this system is further compared with neural network, the balls are nodes in network, force or displacement implemented to certain nodes are data input, and when the velocity spreads to the ﬁxed boundary nodes, the wave will stop and turn to propagate backwards.

In the optimal situation, the complexity of this method is only O(n), where n is the number of nodes. The reason is that the time for the dynamical system to stabilize can be limited and constant, and in each time step, there are no more than 10 forces implemented on each node in 1D case. 10 forces include 8 spring forces, 1 external forces, and 1 artiﬁcial viscosity force.

### 3. Forces in Dynamical system

Here’s a node, whose index is (i,j). Node(i,j) connects to 8 neighbors with spring, and 8 spring forces are deﬁned accordingly.

![](https://upload-images.jianshu.io/upload_images/11683600-cc4239b29b8984a8.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/340)

Figure 1. Schematic diagram of nodes and springs in 2D. Black circle, blue circles and green lines are respectively the node being studied now, neighborhood nodes, and springs connected to central node.

### 4. Pseudo Code

![](https://upload-images.jianshu.io/upload_images/11683600-bba3c78d9e39c87d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### 5. Demonstrations

![](https://upload-images.jianshu.io/upload_images/11683600-d79ef38594accdd0.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### References

[1] Peter E Hammer, Michael S Sacks, J Pedro, and Robert D Howe. Mass-spring model for simulation of heart valve tissue mechanical behavior. Annals of biomedical engineering, 39(6):1668–1679, 2011.

[2] Ullrich Meier, Oscar L´opez, Carlos Monserrat, Mari C Juan, and M Alcaniz. Real-time deformable models for surgery simulation: a survey. Computer methods and programs in biomedicine, 77(3):183–197, 2005.
