import sys

nproc = int(sys.argv[1])
times = []

with open("local_wctimes.temp", "r") as file:
    for line in file:
        t = float(line.split()[3])
        times.append(t)

lt = len(times)
avg = sum(times)/lt
print(f"Average wall clock time using {nproc} processes: {avg} seconds")