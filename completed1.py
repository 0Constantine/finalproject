from scipy.integrate import solve_ivp
from vpython import *
import numpy as np
import matplotlib.pyplot as plt
import vpython as vp


# Define the equations of motion
def triatomic_oscillator(t, y):
    u1, u2, u3, v1, v2, v3 = y

    k = 1.0  # Spring constant
    m1 = 0.25  # Mass 1
    m2 = 0.4  # Mass 2
    m3 = 0.25  # Mass 3

    du1dt = v1
    du2dt = v2
    du3dt = v3

    dv1dt = -(k / m1) * (u1 - u2)
    dv2dt = -(k / m2) * (u2 - u1) - (k / m2) * (u2 - u3)
    dv3dt = -(k / m3) * (u3 - u2)

    return [du1dt, du2dt, du3dt, dv1dt, dv2dt, dv3dt] #1차원배열만취급해줌

# Set initial conditions
initial_conditions = [0., 2., 4., 0.0, 0.0, 0.0] #??ㅋㅋ;

# Integrate the equations of motion
solution = solve_ivp(triatomic_oscillator, [0, 100], initial_conditions, t_eval=np.linspace(0, 100, 20000))

# Extract solution data
t = solution.t
u1, u2, u3, v1, v2, v3 = solution.y



balls0 =vp.sphere(pos = vp.vector(u1[0],0,0),color = color.green, radius =0.3, make_trail = True, retain = 20)
balls1 = vp.sphere(pos = vp.vector(u2[0],0,0),color = color.blue, radius = 0.3, make_trail = True, retain = 20)
balls2 = vp.sphere(pos = vp.vector(u3[0],0,0),color = color.red, radius = 0.3, make_trail = True, retain = 20)
springs0 = vp.helix(pos=vp.vector(u1[0],0,0),axis = vp.vector(u2[0]-u1[0],0,0), radius = 0.3)
springs1 = vp.helix(pos=vp.vector(u2[0],0,0),axis = vp.vector(u3[0]-u2[0],0,0), radius = 0.3)
base = box(pos=vector(0, 0, 0), axis = vector(0,0.1,0),size=vector(10,0.5,10))

# springs = [vp.helix(pos=atoms[i].pos, axis=atoms[i+1].pos-atoms[i].pos, radius=0.05) for i in range(len(atom))]
print('Start')
i = 1
while True:
    rate(30)
    i = i+1
    i = i % len(u1)
    balls0.pos = vp.vector(u1[i],0,0)
    balls1.pos = vp.vector(u2[i],0,0)
    balls2.pos = vp.vector(u3[i],0,0)
    springs0.pos= vp.vector(u1[i],0,0)
    springs0.axis = vp.vector(u2[i]-u1[i],0,0)
    springs1.pos = vp.vector(u2[i],0,0)
    springs1.axis = vp.vector(u3[i]-u2[i],0,0)
