import numpy as np
from scipy.integrate import solve_ivp
import vpython as vp
import matplotlib.pyplot as plt

def N_triatomic_oscillator(t, y):
    [[u], [v]] = y #note : y 가 이제는 

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

    return [du1dt, du2dt, du3dt, dv1dt, dv2dt, dv3dt]

# Set initial conditions
initial_conditions = [-2, 0., 2., 0.0, 0.0, 0.0]

# Integrate the equations of motion
solution = solve_ivp(N_triatomic_oscillator, [0, 10], initial_conditions, t_eval=np.linspace(0, 10, 500))

# Extract solution data
t = solution.t
u1, u2, u3, v1, v2, v3 = solution.y
plt.plot(t, u1, '-')
plt.show()

# Set up the VPython scene
scene = vp.canvas(width=800, height=600, background=vp.color.white)
atoms = [vp.sphere(pos=vp.vector(-3., 0, 0), radius=0.1, color=vp.color.cyan),
         vp.sphere(pos=vp.vector(0, 0, 0), radius=0.1, color=vp.color.magenta),
         vp.sphere(pos=vp.vector(3., 0, 0), radius=0.1, color=vp.color.blue)]

springs = [vp.helix(pos=atoms[i].pos, axis=atoms[i+1].pos-atoms[i].pos, radius=0.05)
           for i in range(len(atoms)-1)]



# while True:
#     rate(1/dt) 
#     # Calculate accelerations of the Lagrangian coordinates=
#     atheta1 = ((E*C/B)*sin(theta1)-F*sin(theta2))/(D-E*A/B)
#     atheta2 = -(A*atheta1+C*sin(theta1))/B
#     # Update velocities of the Lagrangian coordinates=
#     theta1dot = theta1dot+atheta1*dt
#     theta2dot = theta2dot+atheta2*dt
#     # Update Lagrangian coordinates=
#     dtheta1 = theta1dot*dt
#     dtheta2 = theta2dot*dt
#     theta1 = theta1+dtheta1
#     theta2 = theta2+dtheta2
    
#     bar1.rotate( angle=dtheta1, axis=vec(0,0,1), origin=pivot1 )
#     bar1b.rotate( angle=dtheta1, axis=vec(0,0,1), origin=pivot1 )
#     pivot2 = vec(axle2.pos.x, axle2.pos.y, pivot1.z)
#     axle2.rotate( angle=dtheta1, axis=vec(0,0,1), origin=pivot1 )
#     bar2.rotate( angle=dtheta2, axis=vec(0,0,1), origin=pivot2 )
#     pivot2 = vec(axle2.pos.x, axle2.pos.y, pivot1.z)
#     bar2.pos = pivot2 + bar2.axis/2
    
#     t = t+dt
# h= 0.2
# for j in range(len(t)):
#     #j = 1일때,
#     vp.rate(50)
#     A0=atoms[0].pos 
#     print(A0)

#     atoms[1] = A0 + vp.vector(v1[j],0,0)*h
#     atoms[1].pos = atoms[1].pos + vp.vector(v2[j],0,0)*t[j]
#     atoms[2].pos = atoms[2].pos+ vp.vector(v3[j],0,0)*t[j]


#     springs[0].pos = atoms[0].pos
#     springs[1].axis = atoms[2].pos - atoms[1].pos
vp.rate(50)
atoms[1].pos = atoms[0].pos + vp.vector(v1)