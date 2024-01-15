import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
def format(fig):
    mpl.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['xtick.labelsize'] = 19
    plt.rcParams['ytick.labelsize'] = 19
    plt.rcParams['font.size'] = 19
    plt.rcParams['figure.figsize'] = [5.6*6, 4*3]
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 6
    plt.rcParams['legend.fontsize'] = 15
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.linewidth'] = 1.5
    # plt.style.use('dark_background')

def ax_format(ax, xmaj, xmin, ymaj, ymin):
    ax.xaxis.set_tick_params(which='major', size=5, width=1,
                            direction='in', top='on')
    ax.xaxis.set_tick_params(which='minor', size=3, width=1,
                            direction='in', top='on')
    ax.yaxis.set_tick_params(which='major', size=5, width=1,
                            direction='in', right='on')
    ax.yaxis.set_tick_params(which='minor', size=3, width=1,
                            direction='in', right='on')
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(xmaj))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(xmin))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(ymaj))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(ymin))

path = 'vmagF.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
vmagF = np.array(list)
print(vmagF.size)

path = 'vmagI.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
vmagI = np.array(list)
print(vmagI.size)

pi = 3.141592654
boltz = 1.3806e-23    # Boltzmann's constant (J/K)
mass = 6.63e-26       # Mass of argon atom (kg)
T = 273               # Temperature (K)

fig = plt.figure(figsize=(16,8))
format(fig)

ax1 = fig.add_subplot(1, 2, 1)
#ax_format(ax1, 0.5e-13, 0.1e-13, 0.5e-13, 0.1e-13)

hist, bins = np.histogram(vmagI, bins = np.linspace(0.99 * np.min(vmagI), 1.01 * np.max(vmagI), 30))
print(hist, bins)
dx = bins[1] - bins[0]
probability = hist / (dx * np.sum(hist))
ax1.plot(bins[:-1], probability, '.', label = 'sampling', color = 'yellowgreen', markersize = 10)
ax1.set_xlabel('Speed (m/s)')
ax1.set_ylabel('Probability')
ax1.set_title('Initial speed distribution')

ax1 = fig.add_subplot(1, 2, 2)
#ax_format(ax1, 0.5e-13, 0.1e-13, 0.5e-13, 0.1e-13)

hist, bins = np.histogram(vmagF, bins = np.linspace(np.min(vmagF), np.max(vmagF), 30))
dx = bins[1] - bins[0]
probability = hist / (dx * np.sum(hist))
ax1.plot(bins[:-1], 4 * pi * (mass / 2 / pi / boltz / T)**(3/2) * bins[:-1] * bins[:-1] * np.exp(- 0.5 * mass * bins[:-1] * bins[:-1] / boltz / T), label = 'theory', color = 'darkslategray')
ax1.plot(bins[:-1], probability, '.', label = 'sampling', color = 'yellowgreen', markersize = 10)
ax1.set_xlabel('Speed (m/s)')
ax1.set_ylabel('Probability')
ax1.set_title('Final speed distribution')
ax1.legend()

plt.tight_layout()
plt.show()
