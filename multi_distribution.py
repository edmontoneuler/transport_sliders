import numpy as np
import matplotlib.pyplot as plt
from funcs import *
import matplotlib.animation as animation
from matplotlib.widgets import Slider
import time
#plt.ion()

#Constants
kB_eV = 8.617332385e-5
q = 1.609e-19
h_eV = 4.135667e-15
h_J = 6.626068e-34

G0 = 2*q*q/h_J
S0 = -kB_eV
k_e0 = 2*q*kB_eV*kB_eV/h_eV

#Chemical potential/Temperature Grids
dT = 100
T_vec = np.arange(200,800, dT)
dE = 0.01
E = np.arange(-4, 4, dE)
MU, TEMP = np.meshgrid(E, T_vec)

power_list = [0.0, 0.5, 1.0]
moment_list = []

for k in power_list:
    t = base_transport(k, E)
    I0, I1, I2 = moments(t, E, T_vec)
    moment_list.append((t, I0, I1, I2))

gap = 0.38 #Band gap
A=2e8 #Valence band scaling parameter
B=3e8 #Conduction band scaling parameter
kappa_l =2 # Lattice thermal conductivity 
thickness = 1e-9

base_val_t, I0, I1, I2 = moment_list[0]
base_cond_t, J0, J1, J2 = moment_list[2]

V0, V1, V2 = np.flipud(I0), -np.flipud(I1), np.flipud(I2) #Needs minus sign because integral is odd in coordinate inversion
C0, C1, C2 = shift_fermi_int(J0,E, gap), shift_fermi_int(J1, E, gap), shift_fermi_int(J2, E, gap)

K0 = A*V0 + B*C0
K1 = A*V1 + B*C1
K2 = A*V2 + B*C2

t = A*np.flipud(base_val_t) + B*shift_transport(base_cond_t, E, gap)

plt.plot(E, t)
plt.show()

G = G0*K0
S = S0*K1/K0
kappa_e = k_e0*(K2 - K1*K1/K0)*TEMP.T
kappa = kappa_e + kappa_l*np.ones_like(kappa_e)
PF = G*S*S
ZT = TEMP.T*G*S*S/(kappa)

#Initialization
plotlims = [-2, 2, 250,550]
ZT_plot_data = clean_image(ZT, E, T_vec, plotlims=plotlims)

fig, ax_lst = plt.subplots(1,2)

ax_lst[0].imshow(ZT_plot_data, cmap='jet', aspect='auto', extent = plotlims)

#Band Delimiters
ax_lst[0].plot([0,0], [plotlims[2], plotlims[3]], 'r--')
ax_lst[0].plot([gap, gap], [plotlims[2], plotlims[3]], 'r--')

ax_lst[0].set_title('Figure of Merit')
ax_lst[0].set_xlabel('Fermi Level [eV]')
ax_lst[0].set_ylabel('Temperature [K]')

initial_T = 400
T_index = get_index(initial_T, T_vec)
T = T_vec[T_index]
image_line,  = ax_lst[0].plot([plotlims[0], plotlims[1]], [T, T], 'k')

cut, = ax_lst[1].plot(E, ZT[:, T_index])
ax_lst[1].set_ylabel('Figure of Merit')
ax_lst[1].set_xlabel('Fermi Level [eV]')
ax_lst[1].set_xlim(plotlims[0], plotlims[1])
ax_lst[1].plot([0,0], [0, 1.2*np.amax(ZT_plot_data)], 'k--')
ax_lst[1].plot([gap, gap], [0, 1.2*np.amax(ZT_plot_data)], 'k--')
ax_lst[1].set_ylim(0, np.amax(ZT_plot_data))
#  Sliders
T_loc = plt.axes([0.25, 0.00, 0.50, 0.02])
T_amp = Slider(T_loc, 'Temperature [K]', plotlims[2], plotlims[3], valinit = initial_T)
kl_loc = plt.axes([0.25, 0.05, 0.50, 0.02])
kl_amp = Slider(kl_loc, 'Lattice Thermal Conductivity', 0, 1, valinit = 0.5)

def image_update(val):
    time.sleep(0.5)
    T = T_amp.val
    T_index = get_index(T, T_vec)
    T = T_vec[T_index]
    
    kappa_l = kl_amp.val

    kappa = kappa_e + kappa_l*np.ones_like(kappa_e)
    ZT = (TEMP.T)*G*S*S/(kappa)
    ZT_plot_data = clean_image(ZT, E, T_vec, plotlims= plotlims)
    #ax_lst[0].cla()
    ax_lst[0].imshow(ZT_plot_data, cmap='jet', aspect='auto', extent = plotlims)


    image_line.set_ydata([T,T])
    cut.set_ydata(ZT[:, T_index])
    ax_lst[1].set_ylim(0, np.amax(ZT_plot_data))

def temp_update(val):
    T = T_amp.val
    T_index = get_index(T, T_vec)
    T = T_vec[T_index]
    
    image_line.set_ydata([T,T])
    cut.set_ydata(ZT[:, T_index])
    ax_lst[1].set_ylim(0, np.amax(ZT_plot_data))


#fig.canvas.draw_idle()
kl_amp.on_changed(image_update)
T_amp.on_changed(temp_update)
plt.show()

            

