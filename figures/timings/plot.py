# Generating a horizontal bar chart of the speed up of AiiDA-enabled koopmans vs the old implementation

import matplotlib.pyplot as plt
import seaborn as sns

sns.set(context='talk')

timings = {'PBE': 1.6, 'KI': 38.0, 'HSE': 160, 'G₀W₀': 2000}

fig, ax = plt.subplots(figsize=(10, 2.5))
bar_width = 0.6

colors = sns.color_palette(None, 8)
# set marvel red to #ff2600
primary = "#dc005a"
secondary = "#f0f500"

for i, (label, time) in enumerate(timings.items()):
    if label == 'PBE':
        continue
    color = primary # if i % 2 == 1 else secondary
    bottom = i
    ax.barh(label, time / timings['PBE'], bar_width, label=label, left=0, color=color)
    ax.text(time / timings['PBE'] + 10, i - 1.05, f'{time / timings["PBE"]:.0f}×', va='center', ha='left')

ax.set_xlabel('wall time / PBE wall time')
# ax.legend(loc='lower right', ncol=len(times), bbox_to_anchor=(1, 1), frameon=False)
#ax.set_xticklabels([None for _ in ax.get_xticks()])
ax.set_xlim([1, 1.15* max(timings.values()) / timings['PBE']])
plt.tight_layout()
plt.savefig('benchmark.svg')
