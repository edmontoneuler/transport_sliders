import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.widgets import TextBox
from funcs import *

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)
plt.title('Transport Distribution')
plt.xlabel('Fermi Level (eV)')
Egrid = np.linspace(-2.0, 2.0, 1000)

#  First Valence Band
initial_v1_m = 0.0
initial_v1_m_text = str(initial_v1_m)
v1_m = initial_v1_m
initial_v1_A = 1  # units of 1e8
initial_v1_A_text = str(initial_v1_A)
v1_A = initial_v1_A
initial_v1_gap = 0.0
initial_v1_gap_text = str(initial_v1_gap)
v1_gap = initial_v1_gap

#  Second Valence Band
initial_v2_m = 0.0
initial_v2_m_text = str(initial_v2_m)
v2_m = initial_v2_m
initial_v2_A = 1 # units of 1e8
initial_v2_A_text = str(initial_v2_A)
v2_A = initial_v2_A
initial_v2_gap = 0.3
initial_v2_gap_text = str(initial_v2_gap)
v2_gap = initial_v2_gap

#  Combining Bands
v1 = np.array([v1_m, 1e8*v1_A, v1_gap])
v2 = np.array([v2_m, 1e8*v2_A, v2_gap])
vbands = [v1, v2]

#  First Conduction Band
initial_c1_m = 1.0
initial_c1_m_text = str(initial_c1_m)
c1_m = initial_c1_m
initial_c1_A = 1 # units of 1e8
initial_c1_A_text = str(initial_c1_A)
c1_A = initial_c1_A
initial_c1_gap = 0.5
initial_c1_gap_text = str(initial_c1_gap)
c1_gap = initial_c1_gap

#  Second Conduction Band
initial_c2_m = 1.0
initial_c2_m_text = str(initial_c2_m)
c2_m = initial_c2_m
initial_c2_A = 2 # units of 1e8
initial_c2_A_text = str(initial_c2_A)
c2_A = initial_c2_A
initial_c2_gap = 0.8
initial_c2_gap_text = str(initial_c2_gap)
c2_gap = initial_c2_gap

#  Combining Bands
c1 = np.array([c1_m, 1e8*c1_A, c1_gap])
c2 = np.array([c2_m, 1e8*c2_A, c2_gap])
cbands = [c1, c2]

t = transport(Egrid, vbands, cbands)
l, = plt.plot(Egrid, t, lw=2)


def v1_m_update(text):
    global vbands
    global v1_m
    v1_m = eval(text)
    v1 = np.array([v1_m, 1e8*v1_A, v1_gap])
    vbands[0] = v1
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))
    plt.draw()


def v1_A_update(text):
    global vbands
    global v1_A
    v1_A = eval(text)
    v1 = np.array([v1_m, 1e8*v1_A, v1_gap])
    vbands[0] = v1
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))
    plt.draw()


def v1_gap_update(text):
    global v1_gap
    global vbands
    v1_gap = eval(text)
    v1 = np.array([v1_m, 1e8*v1_A, v1_gap])
    vbands[0] = v1
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))
    plt.draw()


def v2_m_update(text):
    global v2_m
    global vbands
    v2_m = eval(text)
    v2 = np.array([v2_m, 1e8*v2_A, v2_gap])
    vbands[1] = v2
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))
    plt.draw()


def v2_A_update(text):
    global v2_A
    global vbands
    v2_A = eval(text)
    v2 = np.array([v2_m, 1e8*v2_A, v2_gap])
    vbands[1] = v2
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))
    plt.draw()


def v2_gap_update(text):
    global v2_gap
    global vbands
    v2_gap = eval(text)
    v2 = np.array([v2_m, 1e8*v2_A, v2_gap])
    vbands[1] = v2
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))
    plt.draw()


def c1_m_update(text):
    global c1_m
    global cbands
    c1_m = eval(text)
    c1 = np.array([c1_m, 1e8*c1_A, c1_gap])
    cbands[0] = c1
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))


def c1_A_update(text):
    global c1_A
    global cbands
    c1_A = eval(text)
    c1 = np.array([c1_m, 1e8*c1_A, c1_gap])
    cbands[0] = c1
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))


def c1_gap_update(text):
    global c1_gap
    global cbands
    c1_gap = eval(text)
    c1 = np.array([c1_m, 1e8*c1_A, c1_gap])
    cbands[0] = c1
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))


def c2_m_update(text):
    global c2_m
    global cbands
    c2_m = eval(text)
    c2 = np.array([c2_m, 1e8*c2_A, c2_gap])
    cbands[1] = c2
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))


def c2_A_update(text):
    global c2_A
    global cbands
    c2_A = eval(text)
    c2 = np.array([c2_m, 1e8*c2_A, c2_gap])
    cbands[1] = c2
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))


def c2_gap_update(text):
    global c2_gap
    global cbands
    c2_gap = eval(text)
    c2 = np.array([c2_m, 1e8*c2_A, c2_gap])
    cbands[1] = c2
    t = transport(Egrid, vbands = vbands, cbands=cbands)
    l.set_ydata(t)
    ax.set_ylim(np.min(t), np.max(t))


v1_m_ax = plt.axes([0.3, 0.08, 0.15, 0.02])
v1_m_box = TextBox(v1_m_ax, 'v1_m', initial=initial_v1_m_text)
v1_m_box.on_submit(v1_m_update)

v1_A_ax = plt.axes([0.3, 0.05, 0.15, 0.02])
v1_A_box = TextBox(v1_A_ax, 'v1_A', initial=initial_v1_A_text)
v1_A_box.on_submit(v1_A_update)

v1_gap_ax = plt.axes([0.3, 0.02, 0.15, 0.02])
v1_gap_box = TextBox(v1_gap_ax, 'v1_gap', initial=initial_v1_gap_text)
v1_gap_box.on_submit(v1_gap_update)

v2_m_ax = plt.axes([0.1, 0.08, 0.15, 0.02])
v2_m_box = TextBox(v2_m_ax, 'v2_m', initial=initial_v2_m_text)
v2_m_box.on_submit(v2_m_update)

v2_A_ax = plt.axes([0.1, 0.05, 0.15, 0.02])
v2_A_box = TextBox(v2_A_ax, 'v2_A', initial=initial_v2_A_text)
v2_A_box.on_submit(v2_A_update)

v2_gap_ax = plt.axes([0.1, 0.02, 0.15, 0.02])
v2_gap_box = TextBox(v2_gap_ax, 'v2_gap', initial=initial_v2_gap_text)
v2_gap_box.on_submit(v2_gap_update)

c1_m_ax = plt.axes([0.6, 0.08, 0.15, 0.02])
c1_m_box = TextBox(c1_m_ax, 'c1_m', initial=initial_c1_m_text)
c1_m_box.on_submit(c1_m_update)

c1_A_ax = plt.axes([0.6, 0.05, 0.15, 0.02])
c1_A_box = TextBox(c1_A_ax, 'c1_A', initial=initial_c1_A_text)
c1_A_box.on_submit(c1_A_update)

c1_gap_ax = plt.axes([0.6, 0.02, 0.15, 0.02])
c1_gap_box = TextBox(c1_gap_ax, 'c1_gap', initial=initial_c1_gap_text)
c1_gap_box.on_submit(c1_gap_update)

c2_m_ax = plt.axes([0.8, 0.08, 0.15, 0.02])
c2_m_box = TextBox(c2_m_ax, 'c2_m', initial=initial_c2_m_text)
c2_m_box.on_submit(c2_m_update)

c2_A_ax = plt.axes([0.8, 0.05, 0.15, 0.02])
c2_A_box = TextBox(c2_A_ax, 'c2_A', initial=initial_c2_A_text)
c2_A_box.on_submit(c2_A_update)

c2_gap_ax = plt.axes([0.8, 0.02, 0.15, 0.02])
c2_gap_box = TextBox(c2_gap_ax, 'c2_gap', initial=initial_c2_gap_text)
c2_gap_box.on_submit(c2_gap_update)

plt.show()
