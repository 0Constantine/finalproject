import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import vpython as vp
from vpython import *
import numpy as np
import pandas as pd

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

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

header = ['u1','v1','u2','v2' ,'u3', 'v3']

data = {'u1': u1, 'v1': v1, 'u2': u2, 'v2': v2, 'u3': u3, 'v3': v3}

# Create a DataFrame from the dictionary
df = pd.DataFrame(data)

# Print the DataFrame




# Animate the motion
fig = plt.figure()




atom = [0]*3




# Set up the VPython scene
scene = vp.canvas(width=800, height=600, background=vp.color.white)
atoms = [vp.sphere(pos=vp.vector(0., 0, 0), radius=0.1, color=vp.color.cyan),
         vp.sphere(pos=vp.vector(3, 0., 0),radius=0.1, color=vp.color.magenta),
         vp.sphere(pos=vp.vector(6., 0, 0), radius=0.1, color=vp.color.blue)]

springs = [vp.helix(pos=atoms[i].pos, axis=atoms[i+1].pos-atoms[i].pos, radius=0.05) for i in range(len(atoms)-1)]


h=0.1
# Animate the motion
#below : I thought only for atom1 designated by an velocity and increase their velocity by increment h,
#But when i saw the motion of this. only blue one(particle 3) can move
# #Why is that?
# for i in range(len(t)):
#     vp.rate(50)
#     A0=atoms[0].pos 
#     print(A0)
#     for j in range(len(atoms)):
#         A0 = A0 + vp.vector(v1[j],0,0)*h
#         atoms[1].pos = atoms[1].pos + vp.vector(v2[j],0,0)*h
#         atoms[2].pos = atoms[2].pos+ vp.vector(v3[j],0,0)*h

#     for j in range(len(springs)):
#         springs[j].pos = atoms[j].pos
#         springs[j].axis = atoms[j+1].pos - atoms[j].pos


###################################################################
#here is another code.
# #Animate version by gpt.
# i = 0
# j = 0
# while j < len(atoms):
#     vp.rate(200)
#     j+=1
#     while i < len(t):
#         i += 1
#         atoms[j].pos = vp.vector(solution.y[2*j][i],0,0)
   
#         springs[j].pos = atoms[j].pos
#         springs[j].axis = atoms[j+1].pos - atoms[j].pos


###################################################
i = 0




while i < len(t):
    vp.rate(200) 
    atoms[0].pos = vp.vector(u1[i],0,0)
    atoms[1].pos = vp.vector(u2[i],0,0)
    atoms[2].pos = vp.vector(u3[i],0,0)

    springs[0].pos = atoms[0].pos
    springs[0].axis = atoms[1].pos - atoms[0].pos
    springs[1].pos = atoms[1].pos
    springs[1].axis = atoms[2].pos-atoms[1].pos

    i += 1