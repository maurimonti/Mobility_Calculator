import numpy as np
# import matplotlib.pyplot as plt
import scipy.constants as cnst
import tkinter as tk
from tkinter import ttk

h = cnst.h
e = cnst.e
Z0 = cnst.physical_constants['characteristic impedance of vacuum'][0]
c = cnst.c
pi = np.pi

N, E, W, S = tk.N, tk.E, tk.W, tk.S


def calculate(*args):
    try:
        P_val = float(P.get())
        R_val = float(R.get())
        d_val = float(d.get())
        wl_val = float(wl.get())

        En.set((P_val * 1e-3 / R_val) * 1e6)
        En_val = float(En.get())

        F.set(En_val / (((d_val * 1e-1)**2) * pi * 2))
        Ep.set(h * c / (wl_val * 1e-9) / e)

        Ep_val = float(Ep.get())
        F_val = float(F.get())

        phi.set((F_val * 1e-6 / (Ep_val * e)) / 1e14)

        phi_val = float(phi.get())
        ni_val = float(ni.get())
        ns_val = float(ns.get())
        dEonE_val = float(dEonE.get())
        # Abs_val = float(Abs.get())
        phi_abs_val = phi_val /  (1 + 0.03**2 / d_val**2)

        mu.set((ni_val + ns_val) * dEonE_val / (Z0 * e * phi_abs_val * 1e14))
    except ValueError:
        pass


root = tk.Tk()
root.title('Example')
mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

P = tk.StringVar()
R = tk.StringVar(value=500)
d = tk.StringVar()
wl = tk.StringVar(value=800)
ni = tk.StringVar(value=1.0)
ns = tk.StringVar(value=2.0)
dEonE = tk.StringVar()
Abs = tk.StringVar(value=3.5)

En = tk.StringVar()
F = tk.StringVar()
Ep = tk.StringVar()
phi = tk.StringVar()
mu = tk.StringVar()


P_entry = ttk.Entry(mainframe, width=7, textvariable=P)
P_entry.grid(column=1, row=1, sticky=(W, E))
R_entry = ttk.Entry(mainframe, width=7, textvariable=R)
R_entry.grid(column=1, row=2, sticky=(W, E))
d_entry = ttk.Entry(mainframe, width=7, textvariable=d)
d_entry.grid(column=1, row=3, sticky=(W, E))
wl_entry = ttk.Entry(mainframe, width=7, textvariable=wl)
wl_entry.grid(column=1, row=4, sticky=(W, E))
ni_entry = ttk.Entry(mainframe, width=7, textvariable=ni)
ni_entry.grid(column=1, row=5, sticky=(W, E))
ns_entry = ttk.Entry(mainframe, width=7, textvariable=ns)
ns_entry.grid(column=1, row=6, sticky=(W, E))
dEonE_entry = ttk.Entry(mainframe, width=7, textvariable=dEonE)
dEonE_entry.grid(column=1, row=7, sticky=(W, E))
# Abs_entry = ttk.Entry(mainframe, width=7, textvariable=Abs)
# Abs_entry.grid(column=1, row=4, sticky=(W, E))

ttk.Label(mainframe,
          textvariable=En).grid(column=3, row=2, sticky=(W, E))
ttk.Label(mainframe,
          textvariable=F).grid(column=3, row=4, sticky=(W, E))
ttk.Label(mainframe,
          textvariable=Ep).grid(column=3, row=6, sticky=(W, E))
ttk.Label(mainframe,
          textvariable=phi).grid(column=3, row=8, sticky=(W, E))
ttk.Label(mainframe,
          textvariable=mu).grid(column=4, row=2, sticky=(W, E))

ttk.Button(mainframe,
           text="Calculate",
           command=calculate).grid(column=4, row=7, sticky=W)


ttk.Label(mainframe, text="Power(mW)").grid(column=2, row=1, sticky=W)
ttk.Label(mainframe, text="Rep rate(Hz)").grid(column=2, row=2, sticky=W)
ttk.Label(mainframe, text="Beam radius(mm)").grid(column=2, row=3, sticky=W)
ttk.Label(mainframe, text="Wavelength(nm)").grid(column=2, row=4, sticky=W)
ttk.Label(mainframe, text="ni").grid(column=2, row=5, sticky=W)
ttk.Label(mainframe, text="ns").grid(column=2, row=6, sticky=W)
ttk.Label(mainframe, text="Delta E /E").grid(column=2, row=7, sticky=W)
# ttk.Label(mainframe, text="Absorbance").grid(column=2, row=8, sticky=W)


# E_entry = ttk.Entry(mainframe, width=7, textvariable=En)
# E_entry.grid(column=3, row=1, sticky=('W', 'E'))
# F_entry = ttk.Entry(mainframe, width=7, textvariable=F)
# F_entry.grid(column=3, row=2, sticky=('W', 'E'))

ttk.Label(mainframe, text="Energy pp(uJ):").grid(column=3, row=1, sticky=W)
ttk.Label(mainframe, text="Fluence(uJ/cm^2):").grid(column=3, row=3, sticky=W)
ttk.Label(mainframe, text="Photon energy(eV):").grid(column=3, row=5, sticky=W)
ttk.Label(mainframe,
          text="Photon flux(1/cm^2) [x 10^14]:").grid(column=3,
                                                      row=7, sticky=W)
ttk.Label(mainframe, text="mu(cm^2/Vs):").grid(column=4, row=1, sticky=W)


for child in mainframe.winfo_children():
    child.grid_configure(padx=5, pady=5)
P_entry.focus()
root.bind('<Return>', calculate)


# d = 1.36  # mm
# R = 500  # Hz
# wl = 800  # nm
# P = 1  # mW
# ni = 1.0
# ns = 2.0
# Re = np.abs((3.4 - 1) / (3.4 + 1))
# A = 3.5
# dEonE = 0.12

# E = P * 1e-3 / R
# print('Energy per pulse:')
# print(str(E * 1e6) + ' uJ\n')
# F = E / ((d * 1e-1)**2 * pi)
# print('Fluence:')
# print(str(F * 1e6) + ' uJ/cm^2\n')
# E_p = h * c / (wl * 1e-9)
# print('Photon Energy:')
# print(str(E_p / e) + ' eV\n')
# phi = F / E_p
# print('Photon flux:')
# print(str(phi * 1e-14) + ' x10$^{14}$ 1/cm^2\n')
# phi_abs = phi * (1 - Re) * (1 - 10**(-A))
# print('Absorbed photons:')
# print(str(phi_abs * 1e-14) + ' x10$^{14}$ 1/cm^2\n')

# # mu = (h * c * R * (d * 1e-1)**2 * pi * (ni + ns)) * dEonE /\
# #     (e * Z0 * wl * 1e-9 * P * 1e-3 * (1 - Re) * (1 - 10**(-A)))

# mu = (ni + ns) * dEonE / (Z0 * e * phi_abs)
# print('Mobility:')
# print(s2tr(mu) + ' cm^2/Vs')

root.mainloop()
