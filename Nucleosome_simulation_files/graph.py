import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('rmsd_core.txt', skiprows=2)
time = data[:,0]  # ps
rms  = data[:,1]  # Å

plt.figure(figsize=(10,6))
plt.plot(time, rms, linewidth=2)
plt.xlabel('Time (ps)', fontsize=14)
plt.ylabel('RMSD (Å)', fontsize=14)
plt.title('DNA RMSD Over 10 ns Production', fontsize=16)
plt.axhline(np.mean(rms), linestyle='--', color='gray')
plt.tight_layout()
plt.savefig('rmsd_core.png', dpi=300)
