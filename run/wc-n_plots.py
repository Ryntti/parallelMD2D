import matplotlib.pyplot as plt

nprocs, avgs = [], []
with open("avg_wctimes.dat", "r") as file:
    for i, line in enumerate(file):
        avgt = float(line.split()[-1])
        np = int(line.split()[7])
        avgs.append(avgt)
        nprocs.append(np)

plt.figure()
plt.title('Wall clock time as a function of the number of processes \n using 16000 atoms')
plt.plot(nprocs, avgs, 'r-', label='wc(N)')
plt.xlabel('Number of processes')
plt.ylabel('Wall clock time (s)')
plt.legend()
plt.savefig('wc-n.png', dpi=100)
