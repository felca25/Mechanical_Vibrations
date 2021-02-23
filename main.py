import numpy as np 
from math import cos, sin, pi
from sympy import Symbol, lambdify
from scipy.integrate import odeint
from matplotlib import pyplot as plt

class Free_1DoF_Vibrational_System:
    mass = 50
    k = 10
    c = 0.001
    delta = 0.0001
    initial_position = 1
    initial_velocity = 0
    force = 0

    def __init__(self):

        self.time = np.arange(0.0, 180.0, self.delta, float)
        self.initialConditions = [self.initial_position, self.initial_velocity]

    @classmethod
    def setInitialConditions(cls, initial_position, initial_velocity):
        cls.initial_position = initial_position
        cls.initial_velocity = initial_velocity

    @classmethod
    def setTimeInterval(cls, inital_time = 0, end_time = 180, delta = 0.001):
        cls.delta = delta
        cls.time = np.arrange(inital_time, end_time, delta, float)
    
    def setQsi(self, qsi):
        self.c = 2 * qsi * np.sqrt(self.k * self.mass)
        self.qsi = qsi

    def setOmega_n(self, omega):
        self.k = self.mass * omega ** 2
        self.omega_n = omega

class NumericVibrationalSystem(Free_1DoF_Vibrational_System):

    def __init__(self):
        super().__init__()

    def FreeVibrationalSystem(self, array, time):
        # Here the default value for the force is defined as 0 N

        # Creating Symbols to use in operationable functions
        f = Symbol('f')
        x = Symbol('x')
        xp = Symbol('xp')

        # redifining self conditions to make it easier to use 
        mass = self.mass
        k = self.k
        c = self.c

        # Defining matrices to calculate easily
        A = np.array([
            [0, 1],
            [-k/mass, -c/mass]
        ])
        X_matrix = np.array([
            [x],
            [xp]
        ])

        Force_matrix = np.array([
            [0],
            [f/mass]
        ])
        #Defining the function that will define the free vibrational system
        Y_matrix = np.matmul(A, X_matrix) + Force_matrix
        # This will give us a matrix with mathematical symbols, which means that if we want to we can't pass
        # a variable inside the function to get a proper value, so we need to convert the mathsymbols function
        # to a usable function


        # Converting matrix into an usable function
        y = lambdify([x, xp], Y_matrix[0])
        yp = lambdify([x, xp, f], Y_matrix[1])


        u = array[0]
        up = array[-1]

        return [y(u,up)[-1], yp(u, up, self.force)[-1]]
        # [-1] to transform the array given by each element into float value

    def NumericSolution(self):
        #return the result of odeint calculations
        return odeint(self.FreeVibrationalSystem, self.initialConditions, self.time)

    def plotPositionGraphs(self):
        plt.figure()
        plt.subplot(211)
        plt.plot(self.time, self.NumericSolution()[:,0])
        plt.xlabel('Time(s)')
        plt.ylabel('Deslocamento(m)')
        plt.subplot(212)
        plt.plot(self.time, self.NumericSolution()[:,1])
        plt.xlabel('Time(s)')
        plt.ylabel('Velocidade(m/s)')
        # plt.show()

class AnaliticalSolution(Free_1DoF_Vibrational_System):

    def __init__(self):
        super().__init__()

    def FreeVibrationalSystem(self):
        
        if self.qsi == 0:
            A = 1 #amplitude
            phi = 0

        



System = NumericVibrationalSystem()
System.setQsi(0.1)
System.plotPositionGraphs()

systemUndamped = NumericVibrationalSystem()
systemUndamped.setQsi(0)
systemUndamped.plotPositionGraphs()

plt.show()


