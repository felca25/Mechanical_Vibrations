import numpy as np 
from math import cos, sin, pi, e
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
    qsi = 0
    omega_n = np.sqrt(k/mass)

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

    @classmethod
    def setOmega_n(cls, omega = None):
        if omega is not None:
            cls.k = cls.mass * omega ** 2
            cls.omega_n = omega
        else:
            cls.omega_n = np.sqrt(cls.k / cls.mass)
    @classmethod
    def setQsi(cls, qsi):
        cls.c = 2 * qsi * np.sqrt(cls.k * cls.mass)
        cls.qsi = qsi

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

class AnaliticalSystem(Free_1DoF_Vibrational_System):

    def __init__(self):
        super().__init__()
        self.omega_d = self.omega_n * np.sqrt(1 - (self.qsi)**2)

    def VibrationalSystem(self):

        t = Symbol('t')
        #Free Vibrational System
        if self.qsi == 0:

            A = np.sqrt((self.initial_position ** 2) + ((self.initial_velocity / self.omega_n)**2)) #amplitude
            phi = np.arctan(self.initial_velocity/(self.initial_position * self.omega_n))
            
            #analitical solution
            self.u = A * cos(self.omega_n * t - phi)

            u_func = lambdify(t, self.u)

            return u_func
        #Underdamped Vibrational System
        elif 0 < self.qsi < 1 and self.qsi > 0:

            exponential_part = e ** (- self.qsi * self.omega_n * t)
            c1 = self.initial_position
            c2 = (self.initial_velocity + (self.initial_position * self.qsi * self.omega_n)) / self.omega_d
            
            #analitical solution
            self.u = exponential_part * (c1 * cos(self.omega_d * t) + c2 * sin(self.omega_d * t))

            u_func = lambdify(t, self.u)

            return u_func
        #Critically Damped Vibrational System
        elif self.qsi == 1:

            A1 = self.initial_position
            A2 = self.initial_velocity + (A1 * self.omega_n)
            exponential_part = np.exp(- self.omega_n * t)

            self.u = exponential_part * (A1 + A2 * t)

            u_func = lambdify(t, self.u)

            return u_func
        #Overdamped System
        elif self.qsi > 1:

            omega_d = self.omega_n * np.sqrt(self.qsi ** 2 - 1)

            exponential_part = np.exp(- self.qsi * self.omega_n * t)
            C1 = self.initial_position
            C2 = (self.initial_velocity + (C1 * self.qsi * self.omega_n))/(omega_d)

            self.u = exponential_part * (C1 * np.cosh(omega_d * t) + C2 * np.sinh(omega_d * t))

            u_func = lambdify(t, self.u)

            return u_func
        else:
            pass
    def plotPositionGraphs(self):
        plt.figure()
        plt.subplot(211)
        plt.plot(self.time, self.VibrationalSystem())
        plt.xlabel('Time(s)')
        plt.ylabel('Deslocamento(m)')
        plt.subplot(212)
        plt.plot(self.time, self.VibrationalSystem())
        plt.xlabel('Time(s)')
        plt.ylabel('Velocidade(m/s)')


# ----- Portugues -------
# para utilizar esse código faz se com o seguinte exemplo:
# Funcional:
# 
# System_Numeric = NumericVibrationalSystem()
# System_Numeric.setQsi(0.1)
# 
# 
# Primeiro criamos uma variavel que será o nosso sistema, e ele vai ser criado
# baseado no tipo de análise desejada, por exemplo, se quisermos fazer uma análise 
# analítica utilizamos a classe AnaliticalSolution() para iniciar a nossa variável
# dessa forma nós temos um 



System = AnaliticalSystem()
System.setQsi(0.1)
System.plotPositionGraphs()

systemUndamped = NumericVibrationalSystem()
systemUndamped.setQsi(0)
systemUndamped.plotPositionGraphs()

plt.show()