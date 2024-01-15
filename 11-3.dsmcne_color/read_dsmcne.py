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

path = 'xcell.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
xcell = np.array(list)
print(xcell.size)

path = 'ave_n.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
ave_n = np.array(list)
print(ave_n.size)

path = 'ave_ux.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
ave_ux = np.array(list)
print(ave_ux.size)

path = 'ave_uy.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
ave_uy = np.array(list)
print(ave_uy.size)

path = 'ave_uz.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
ave_uz = np.array(list)
print(ave_uz.size)

path = 'ave_T.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
ave_T = np.array(list)
print(ave_T.size)

path = 'white.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
white = np.array(list)
print(white.size)

path = 'black.txt'
list = []
fo = open(path, 'r')
for line in fo:
    list.append(float(line))
fo.close()
black = np.array(list)
print(black.size)

fig = plt.figure(figsize=(16,8))
format(fig)

ax1 = fig.add_subplot(2, 2, 1)
#ax_format(ax1, 0.5e-13, 0.1e-13, 0.5e-13, 0.1e-13)

ax1.plot(xcell, ave_n)
ax1.set_xlabel('position')
ax1.set_ylabel('Number density')

ax1 = fig.add_subplot(2, 2, 2)
#ax_format(ax1, 0.5e-13, 0.1e-13, 0.5e-13, 0.1e-13)

ax1.plot(xcell, ave_ux, label = 'x-component')
ax1.plot(xcell, ave_uy, label = 'y-component')
ax1.plot(xcell, ave_uz, label = 'z-component')
ax1.set_xlabel('position')
ax1.set_ylabel('Velocities')
ax1.legend()

ax1 = fig.add_subplot(2, 2, 3)
#ax_format(ax1, 0.5e-13, 0.1e-13, 0.5e-13, 0.1e-13)

ax1.plot(xcell, ave_T)
ax1.set_xlabel('position')
ax1.set_ylabel('Temperature')

ax1 = fig.add_subplot(2, 2, 4)
#ax_format(ax1, 0.5e-13, 0.1e-13, 0.5e-13, 0.1e-13)

ax1.plot(xcell, black, '.', label = 'black')
ax1.plot(xcell, white, '.', label = 'white')
ax1.set_xlabel('position')
ax1.set_ylabel('number')
ax1.legend()

plt.tight_layout()
plt.show()
