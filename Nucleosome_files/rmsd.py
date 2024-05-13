import numpy as np
import matplotlib.pyplot as plt
rmsd = np.loadtxt('rms_Nucleosome_no_tails.txt', skiprows=2)
plt.figure(figsize=(12,8))
plt.plot(rmsd[:,0]/1000, rmsd[:,1], label = 'Major', color='#6a137a', linewidth=2)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0,10)
plt.xlabel('Time (ns)', fontsize=30)
plt.ylabel('RMSD nucleosome ($\AA$)', fontsize=30)
plt.axhline(np.average(rmsd[:,1]), color='#2E7D32', linewidth=3)
plt.savefig('rmsd_nucleosome.png', bbox_inches='tight')
print(np.average(rmsd[:,1]))