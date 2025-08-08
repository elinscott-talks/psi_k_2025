# Generating a horizontal bar chart of the speed up of AiiDA-enabled koopmans vs the old implementation

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set(context='talk')

grid = [2, 3, 4, 5, 6]
fbz = [g ** 3 for g in grid]
ibz = [4, 6, 13, 19, 32]
factor = [y / x for x, y in zip(fbz, ibz)]
time_with = [
    2.600000000000000000e+01,
    1.320000000000000000e+02,
    6.750000000000000000e+02,
    2.309000000000000000e+03,
    1.020000000000000000e+04,
]

time_without = [
	6.400000000000000000e+01,
	7.660000000000000000e+02,
	4.380000000000000000e+03,
	1.740000000000000000e+04,
	5.124000000000000000e+04,
]
speedup = [x / y for x, y in zip(time_without, time_with)]

min_nk = [
    4.000000000000000000e+00,
    6.000000000000000000e+00,
    1.300000000000000000e+01,
    1.900000000000000000e+01,
    3.200000000000000000e+01,
]
max_nk = [
    1.600000000000000000e+01,
    5.400000000000000000e+01,
    1.280000000000000000e+02,
    2.500000000000000000e+02,
    4.320000000000000000e+02,
]

max_nk = [x / 2 for x in max_nk]

factor2 = [[x / y for x, y in zip(nk, fbz)] for nk in [min_nk, max_nk]]

plt.rcParams['axes.titlepad'] = 20

print(factor2)
def plot(y, title, filename, ylabel='fraction of explicit $k$-points', ylim=[0, 1]):
    fig, ax = plt.subplots(figsize=(5, 5))
    bar_width = 0.6

    colors = sns.color_palette(None, 8)
    # ax.plot(grid, fbz, 'x', label='no symmetry')
    # ax.plot(grid, ibz, '+', label='with symmetry')
    if isinstance(y[0], list):
        ax.fill_between(grid, y[0], y[1], color="#ff8d7d")
    else:
        ax.plot(grid, y, color="#ff2600")

    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xticks(grid, [rf"${x}\times{x}\times{x}$" for x in grid])
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_xlim(min(grid), max(grid))
    ax.set_ylim(ylim)

    #ax.legend(loc='lower right', ncols=len(times), bbox_to_anchor=(1, 1), frameon=False)
    plt.tight_layout()
    plt.savefig(filename)

plot(factor, 'outer loop (symmetry of perturbation)', 'bz-to-ibz-outer.svg', 'fraction of explicit $q$-points')
plot(factor2, 'inner loop (little point of $q$)', 'bz-to-ibz-inner.svg')
plot(speedup, None, 'bz-to-ibz-speedup.svg', 'speed-up', ylim=[2, 8])
