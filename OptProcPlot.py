import matplotlib.pyplot as plt
import re
import numpy as np

def parse_log(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    scores = {'initialNRSplit': [], 'upperBound': [], 'lowerBound': [], 'allocStep': []}
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith('Optimization Cycle'):
            for j in range(4):
                i += 1
                line = lines[i].strip()
                match = re.match(r'Optimized (\w+): ([\d.]+) \(Score: ([\d.]+)\)', line)
                if match:
                    param = match.group(1)
                    score = float(match.group(3))
                    if param in scores:
                        scores[param].append(score)
        i += 1
    return scores

def parse_best_params(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    best_params = {}
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith('Best Parameters:'):
            i += 1  # Skip the 'Best Parameters:' line
            for j in range(4):
                if i < len(lines):
                    line = lines[i].strip()
                    match = re.match(r'(\w+): ([\d.]+)', line)
                    if match:
                        param = match.group(1)
                        value = float(match.group(2))
                        best_params[param] = value
                    i += 1
        else:
            i += 1
    return best_params

def parse_best_score(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        match = re.search(r'Best Objective Score: ([\d.]+)', line)
        if match:
            return float(match.group(1))
    return None

def plot_scores(ax, scores, title):
    colors = ['b', 'g', 'r', 'c']
    for param, color in zip(scores.keys(), colors):
        ax.plot(range(1, 21), scores[param], label=param, color=color)
    ax.set_xlabel('Cycle Number')
    ax.set_ylabel('Score')
    ax.set_title(title)
    ax.set_xticks([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    ax.legend(loc='upper right')
    ax.grid(True)

dynamic_scores = parse_log('result.log')
static_scores = parse_log('static.log')

fig1, ax1 = plt.subplots(figsize=(20, 12))
plot_scores(ax1, dynamic_scores, 'Dynamic Optimization Process')
plt.show()

fig2, ax2 = plt.subplots(figsize=(20, 12))
plot_scores(ax2, static_scores, 'Static Optimization Process')
plt.show()

dynamic_best = parse_best_params('result.log')
static_best = parse_best_params('static.log')
dynamic_best_score = parse_best_score('result.log')
static_best_score = parse_best_score('static.log')

params = ['initialNRSplit', 'upperBound', 'lowerBound', 'allocStep', 'bestScore']
x = np.arange(len(params))  # the label locations
width = 0.3  # the width of the bars

fig3, ax3 = plt.subplots(figsize=(24, 12))
ax3.bar(x[:-1] - width/2, [static_best[p] for p in params[:-1]], width, label='Static', color='blue')
ax3.bar(x[:-1] + width/2, [dynamic_best[p] for p in params[:-1]], width, label='Dynamic', color='orange')
ax3.set_ylabel('Parameter Value')
ax3.set_ylim(0, 1)

ax_right = ax3.twinx()
ax_right.bar(x[-1] - width/2, static_best_score, width, color='blue')
ax_right.bar(x[-1] + width/2, dynamic_best_score, width, color='orange')
ax_right.set_ylabel('Best Score')
ax_right.set_ylim(4500, 5500)

ax3.set_xlabel('Parameters')
ax3.set_title('Comparison of Optimized Parameters and Best Scores: Static vs Dynamic')
ax3.set_xticks(x)
ax3.set_xticklabels(params)
ax3.legend()
ax3.grid(True)
plt.show()