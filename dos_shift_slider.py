import numpy as np
import matplotlib.pyplot as plt
from funcs import *
from matplotlib.widgets import Slider

N = 1000
egrid = np.linspace(-3, 3, N)
dE = egrid[2] - egrid[1]

#Product of M(E) and Landauer velocity quantity 
base_quantity = transport(egrid)
band_gap = find_bandgap(base_quantity, egrid)

#Density of states data
raw_dos = np.zeros_like(egrid)

for k in range(N):
    if egrid[k] < 0:
        raw_dos[k] = np.exp(-egrid[k]/0.5) - 1
    elif egrid[k] > band_gap:
        raw_dos[k] = np.exp((egrid[k] - band_gap)/0.5) - 1
    else:
        raw_dos[k] = 0

raw_dos = 1e6*raw_dos

def dos_shift(egrid, dos, numerator, v_shift, c_shift):
    N = len(egrid)
    dE = egrid[2] - egrid[1]
    v_index_shift = int(v_shift/dE)
    c_index_shift = int(c_shift/dE)

    for k in range(N):
        if dos[k] ==0:
            v_edge_index = k
            break
    new_egrid = egrid[v_index_shift:N-c_index_shift-1]
    new_numerator = numerator[v_index_shift:N-c_index_shift-1]
    indices = [k for k in range(v_edge_index-1, v_edge_index+v_index_shift+c_index_shift)]
    new_dos = np.delete(dos, indices)

    return new_egrid, new_dos, new_numerator 
  

fig, ax_lst = plt.subplots(1,2)

dos, = ax_lst[0].plot(egrid, raw_dos, 'k', lw=2)
T, = ax_lst[1].plot(egrid, base_quantity/raw_dos, 'k', lw=2 )

ax_lst[0].set_title('Density of States')
ax_lst[0].set_xlabel('Energy Level [eV]')
ax_lst[0].set_ylabel('DOS (arb. units)')
ax_lst[0].plot(egrid, base_quantity, 'g', lw=2)
ax_lst[0].axvline(x=0, linestyle = '--')
ax_lst[0].axvline(x=band_gap, linestyle='--')

ax_lst[1].set_title('Transport Distribution')
ax_lst[1].set_xlabel('Energy Level [eV]')
ax_lst[1].set_ylabel('T(E) (arb. units)')
ax_lst[1].axvline(x=0, linestyle='--')
ax_lst[1].axvline(x=band_gap, linestyle='--')

v_shift_loc = plt.axes([0.25, 0.00, 0.50, 0.02])
v_shift_amp = Slider(v_shift_loc, 'Valence Band shift (eV)', 0, 0.2, valinit = 0.02)
c_shift_loc = plt.axes([0.25, 0.025, 0.50, 0.02])
c_shift_amp = Slider(c_shift_loc,'Conduction Band shift (eV)', 0, 0.2, valinit = 0.2)

def update(val):
    v_shift = v_shift_amp.val
    c_shift = c_shift_amp.val

    new_egrid, new_dos, new_numerator = dos_shift(egrid, raw_dos, base_quantity, v_shift, c_shift)
    dos.set_data(new_egrid, new_dos)
    T.set_data(new_egrid, new_numerator/new_dos)
    
v_shift_amp.on_changed(update)
c_shift_amp.on_changed(update)
plt.show()
