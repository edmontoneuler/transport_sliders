import numpy as np
import matplotlib.pyplot as plt
from funcs import *
from matplotlib.widgets import Slider
import time

start = time.time()

#Constants
kB_eV = 8.617332385e-5
q = 1.609e-19
h_eV = 4.135667e-15
h_J = 6.626068e-34

G0 = 2*q*q/h_J
S0 = -kB_eV
k_e0 = 2*q*kB_eV*kB_eV/h_eV

#Chemical potential/Temperature Grids
dT = 20
T_vec = np.arange(200,1200, dT)
dE = 0.002
E = np.arange(-3.5, 3.5, dE)
plotlims = [-2, 2, 300, 1000] #should leave at least 1.5eV of buffer (integration and gap shifting)
MU, TEMP = np.meshgrid(E, T_vec)
kappa_l = 0.05
m_list = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
moment0_list = []
moment1_list =[]
moment2_list = []

#  Calculation of basic Fermi moments
for k in range(len(m_list)):
    I0, I1, I2 = moments(transport(E, vbands=[], cbands=[np.array([m_list[k], 1, 0])]), E, T_vec)
    moment0_list.append(I0)
    moment1_list.append(I1)
    moment2_list.append(I2)

#  Initial parameters

v1_A = 1
v1_m = 0.0
v1_gap = 0.0
v1 = np.array([v1_m, v1_A, v1_gap])

c1_A = 1
c1_m = 1.0
c1_gap = 0.5
c1 = np.array([c1_m, c1_A, c1_gap])

# TE Calculation from moment lists
v1_m_index = get_index(v1_m, m_list)
c1_m_index = get_index(c1_m, m_list)

V0 = moment0_list[v1_m_index]
V1 = moment1_list[v1_m_index]
V2 = moment2_list[v1_m_index]

C0 = moment0_list[c1_m_index]
C1 = moment1_list[c1_m_index]
C2 = moment2_list[c1_m_index]

V0, V1, V2 = np.flipud(V0), -np.flipud(V1), np.flipud(V2) # Needs minus sign because integral is odd under coordinate inversion
C0, C1, C2 = shift_fermi_int(C0, E, c1_gap), shift_fermi_int(C1, E, c1_gap), shift_fermi_int(C2, E, c1_gap)

K0 = v1_A*V0 + c1_A*C0
K1 = v1_A*V1 + c1_A*C1
K2 = v1_A*V2 + c1_A*C2

G = 1e8*G0*K0
S = S0*K1/K0
kappa_e = 1e8*k_e0*(K2 - K1*K1/K0)*TEMP.T
kappa = kappa_e + kappa_l*np.ones_like(kappa_e)
PF = G*S*S
ZT = TEMP.T*PF/kappa

#  Initialization

ZT_plot_data = clean_image(ZT, E, T_vec, plotlims=plotlims)
fig, ax_lst = plt.subplots(2,2)
ax_lst[1, 0].axis('off')

ax_lst[0,0].imshow(ZT_plot_data, cmap='jet', aspect='auto', extent=plotlims)

#Band Delimiters
v_edge, = ax_lst[0,0].plot([0,0], [plotlims[2], plotlims[3]], 'r--')
c_edge, = ax_lst[0,0].plot([c1_gap, c1_gap], [plotlims[2], plotlims[3]], 'r--')

ax_lst[0,0].set_title('Figure of Merit')
ax_lst[0,0].set_xlabel('Fermi Level [eV]')
ax_lst[0,0].set_ylabel('Temperature [K]')

initial_T = 400
T_index = get_index(initial_T, T_vec)
T = T_vec[T_index]
image_line,  = ax_lst[0,0].plot([plotlims[0], plotlims[1]], [T, T], 'k')

cut, = ax_lst[0,1].plot(E, ZT[:, T_index])
ax_lst[0, 1].set_ylabel('Figure of Merit')
ax_lst[0, 1].set_xlabel('Fermi Level [eV]')
ax_lst[0, 1].set_xlim(plotlims[0], plotlims[1])
ax_lst[0, 1].plot([0,0], [0, 1.2*np.amax(ZT_plot_data)], 'k--')
cut_cband_edge, = ax_lst[0,1].plot([c1_gap, c1_gap], [0, 1.2*np.amax(ZT_plot_data)], 'k--')
cut_vband_edge, = ax_lst[0,1].plot([0, 0], [0, 1.2*np.amax(ZT_plot_data)], 'k--')
ax_lst[0, 1].set_ylim(0, np.amax(ZT_plot_data))

t,  = ax_lst[1,1].plot(E, transport(E, vbands=[v1], cbands=[c1]), lw=2)
ax_lst[1, 1].set_ylabel('Transport Distribution [1e8]')
ax_lst[1, 1].set_xlabel('Fermi Energy (eV)')
ax_lst[1, 1].set_xlim(plotlims[0], plotlims[1])

T_loc = plt.axes([0.15, 0.40, 0.30, 0.02])
T_amp = Slider(T_loc, 'Temperature [K]', plotlims[2], plotlims[3], valinit = initial_T)

c1_gap_loc = plt.axes([0.15, 0.15, 0.30, 0.02])
c1_gap_amp = Slider(c1_gap_loc, 'Band gap [eV]', 0, 1, valinit=c1_gap)

c1_A_loc = plt.axes([0.15, 0.20, 0.30, 0.02])
c1_A_amp = Slider(c1_A_loc, 'Conduction Band Amplitude', 0.1, 10, valinit=c1_A)

v1_A_loc = plt.axes([0.15, 0.25, 0.30, 0.02])
v1_A_amp = Slider(v1_A_loc, 'Valence Band Amplitude', 0.1, 10, valinit=v1_A)

c1_m_loc = plt.axes([0.15, 0.30, 0.30, 0.02])
c1_m_amp = Slider(c1_m_loc, 'Conduction Band Power Law', min(m_list), max(m_list), valinit=c1_m)

v1_m_loc = plt.axes([0.15, 0.35, 0.30, 0.02])
v1_m_amp = Slider(v1_m_loc, 'Valence Band Power Law', min(m_list), max(m_list), valinit=v1_m)

end = time.time()

print('Elapsed time is ', end-start, ' seconds')

def param_update(val):
    global v1
    global c1
    global ZT

    c1_gap = c1_gap_amp.val
    c1[2] = c1_gap
    c1_A = c1_A_amp.val
    c1[1] = c1_A
    v1_A = v1_A_amp.val
    v1[1] = v1_A
    c1_m = c1_m_amp.val
    v1_m = v1_m_amp.val

    v1_m_index = get_index(v1_m, m_list)
    c1_m_index = get_index(c1_m, m_list)
    v1_m = m_list[v1_m_index]
    c1_m = m_list[c1_m_index]
    v1[0] = v1_m
    c1[0] = c1_m

    V0 = moment0_list[v1_m_index]
    V1 = moment1_list[v1_m_index]
    V2 = moment2_list[v1_m_index]

    C0 = moment0_list[c1_m_index]
    C1 = moment1_list[c1_m_index]
    C2 = moment2_list[c1_m_index]

    V0, V1, V2 = np.flipud(V0), -np.flipud(V1), np.flipud(V2)  # Needs minus sign because integral is odd under coordinate inversion
    C0, C1, C2 = shift_fermi_int(C0, E, c1_gap), shift_fermi_int(C1, E, c1_gap), shift_fermi_int(C2, E, c1_gap)
    K0 = v1_A * V0 + c1_A * C0
    K1 = v1_A * V1 + c1_A * C1
    K2 = v1_A * V2 + c1_A * C2

    G = 1e8*G0 * K0
    S = S0 * K1 / K0
    kappa_e = 1e8*k_e0 * (K2 - K1 * K1 / K0) * TEMP.T
    kappa = kappa_e + kappa_l * np.ones_like(kappa_e)
    PF = G * S * S
    ZT = TEMP.T * PF / kappa

    ZT_plot_data = clean_image(ZT, E, T_vec, plotlims=plotlims)
    ax_lst[0, 0].imshow(ZT_plot_data, cmap='jet', aspect='auto', extent=plotlims)
    v_edge.set_xdata([0,0])
    c_edge.set_xdata([c1_gap, c1_gap])
    cut.set_ydata(ZT[:, T_index])

    cut_cband_edge.set_xdata([c1_gap, c1_gap])
    cut_cband_edge.set_ydata([0, 1.2*np.amax(ZT)])
    cut_vband_edge.set_ydata([0, 1.2*np.amax(ZT)])
    ax_lst[0, 1].set_ylim(0, 1.1*np.amax(ZT))

    trans = transport(E, vbands=[v1], cbands=[c1])
    t.set_ydata(trans)
    ax_lst[1,1].set_ylim([0, 1.1*max(trans)])

def temp_update(val):
    global T
    global T_index
    global ZT
    T = T_amp.val
    T_index = get_index(T, T_vec)
    T = T_vec[T_index]

    image_line.set_ydata([T, T])
    cut.set_ydata(ZT[:, T_index])
    ax_lst[0,1].set_ylim(0, np.amax(ZT))

T_amp.on_changed(temp_update)
c1_gap_amp.on_changed(param_update)
c1_A_amp.on_changed(param_update)
v1_A_amp.on_changed(param_update)
c1_m_amp.on_changed(param_update)
v1_m_amp.on_changed(param_update)
plt.show()