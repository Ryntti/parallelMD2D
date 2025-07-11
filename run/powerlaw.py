import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear_fit(x, a, b):
    return a*x + b

nprocs, avgs = [], []
with open("avg_puhtitimes.dat", "r") as file:
    for i, line in enumerate(file):
        avgt = float(line.split()[-1])
        nproc = int(line.split()[7])
        avgs.append(avgt)
        nprocs.append(nproc)


# Prepare data for fitting
log_np = np.log10(nprocs)  # logarithmized data for fitting
log_avgs = np.log10(avgs)

# linear fit on the logarithmized data
popt, pcov = curve_fit(linear_fit, log_np, log_avgs)
a, b = popt

fit_x = np.logspace(np.log10(min(nprocs)), np.log10(max(nprocs)), 100)
fit_y = 10**b * fit_x**a  # switch back to original scale

print()

plt.figure()
plt.title(f't(N) as a loglog plot \n power law a = {a:.5f}')

plt.loglog(nprocs, avgs, 'r-', label='t(N)') 
plt.loglog(fit_x, fit_y, 'b-', label='Line fit', alpha=0.6) 

plt.xlabel('Number of processes')
plt.ylabel('Wall clock time (s)')
plt.legend()
plt.savefig('wc-n.png', dpi=100)
