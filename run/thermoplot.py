import numpy as np
import matplotlib.pyplot as plt

t, tot_e, ep, ek = [], [], [], []

with open("thermo.dat", mode='r') as file:
    lines = file.readlines()  # Read all lines into a list
    for line in lines[:-1]:  # Iterate over all but the last line
        str_line = line.split()
        t.append(float(str_line[0]))
        tot_e.append(float(str_line[1]))
        ep.append(float(str_line[2]))
        ek.append(float(str_line[3]))


t = np.array(t)
tot_e = np.array(tot_e)
ep = np.array(ep)
ek = np.array(ek)

plt.figure()
plt.title('Internal energy, potential energy \n and kinetic energy as functions of time')
plt.plot(t, tot_e, 'g-', label='Internal energy')
plt.plot(t, ep, 'b-', label='potential energy', alpha=0.7)
plt.plot(t, ek, 'r-', label='kinetic energy', alpha=0.7)
plt.xlabel('Time (ps*)')
plt.ylabel('Energy (eV*)')
plt.legend()
plt.savefig('Thermodata_plot.png', dpi=100)