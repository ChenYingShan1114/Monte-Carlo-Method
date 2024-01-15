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

path = 'rand.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
rand = np.array(list)
print(rand.size)

path = 'N.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
N = np.array(list)
print(N.size)

path = 'P.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
P = np.array(list)
print(P.size)

fig = plt.figure(figsize=(16,8))
format(fig)

ax1 = fig.add_subplot(1, 2, 1)
#ax_format(ax1, 0.5e-13, 0.1e-13, 0.5e-13, 0.1e-13)

ax1.plot(rand, '.', markersize = 1, color = 'dimgray')
ax1.set_xlabel('point number')
ax1.set_ylabel('random number')
ax1.set_title('sampling point')

ax1 = fig.add_subplot(1, 2, 2)
#ax_format(ax1, 0.5e-13, 0.1e-13, 0.5e-13, 0.1e-13)

ax1.plot(N, P, label = 'theory', color = 'darkslategray')
hist, bins = np.histogram(rand, bins = np.linspace(np.min(rand), np.max(rand), 100))
#print(hist, bins)
dx = bins[1] - bins[0]
probability = hist / (dx * np.sum(hist))
#ax1.plot((bins[:-1]+bins[1:])/2, probability, '.', label = 'sampling', color = 'yellowgreen', markersize = 10)
ax1.plot((bins[:-1]), probability, '.', label = 'sampling', color = 'yellowgreen', markersize = 10)
ax1.set_xlim(np.min(rand), np.max(rand))
ax1.set_xlabel('random number')
ax1.set_ylabel('probability')
ax1.set_title('probability distribution')
ax1.legend()

plt.tight_layout()
fig.suptitle(r'Use acceptance-rejection generate distribution between [-1,1] of $P_x(x)dx = \frac{3}{4}(1-x^2)dx$')
#fig.suptitle(r'Binomial distribition $P_k[k] = \frac{N!}{k!(N-k)!} \left( \frac{1}{2} \right)^N$')
#fig.suptitle(r'Poisson distribition $P_i[i] = \frac{e^{-\lambda}\lambda^i}{i!}$, where $\lambda = 10$')
plt.tight_layout()
#plt.savefig('acceptance-rejection.png')
#plt.savefig('binomial-distribution.png')
#plt.savefig('poisson-distribution.png')
plt.show()
