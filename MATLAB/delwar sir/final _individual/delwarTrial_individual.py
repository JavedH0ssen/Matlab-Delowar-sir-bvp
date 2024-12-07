import matplotlib.pyplot as plt
import numpy as np

class Solution:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def jobs(param):
    x = np.linspace(0, 10, 100)
    y = np.array([np.sin(x + param), np.cos(x + param), np.tan(x + param), np.exp(x + param)])
    return Solution(x, y)

def addArrowAnnotations(fig):
    # Example implementation for adding arrow annotations
    ax = fig.gca()
    ax.annotate('Annotation', xy=(5, 0.5), xytext=(7, 0.7),
                arrowprops=dict(facecolor='black', shrink=0.05))

def plotting(line, titles, sol, sol1, sol2, aa, bb, cc, name):
    plt.plot(sol.x, sol.y[line, :], '-r', label=f'{name} = {aa}')
    plt.plot(sol1.x, sol1.y[line, :], '--g', label=f'{name} = {bb}')
    plt.plot(sol2.x, sol2.y[line, :], '-.b', label=f'{name} = {cc}')
    plt.legend(fontsize=11.5, loc='upper right')
    data = np.concatenate((sol.y[line, :], sol1.y[line, :], sol2.y[line, :]))
    minValue = np.min(data)
    maxValue = np.max(data)
    plt.ylim([minValue, maxValue + (maxValue - minValue) / 100])
    plt.xlabel(r'$\eta$', fontsize=13, fontweight='bold')
    plt.ylabel(titles, fontsize=13, fontweight='bold')
    plt.gca().tick_params(axis='both', which='major', labelsize=10)
    plt.gca().yaxis.set_label_coords(-0.02, 0.5)
    plt.box(True)

def main():
    aa = 1
    bb = 2
    cc = 3
    sol = jobs(aa)
    sol1 = jobs(bb)
    sol2 = jobs(cc)

    name = '\amma'
    ame = name.replace('\lamma', '')

    fig1 = plt.figure(1)
    plotting(2, "f'", sol, sol1, sol2, aa, bb, cc, name)
    plt.xlim([0, 5])
    addArrowAnnotations(fig1)

    fig2 = plt.figure(2)
    plotting(3, "g", sol, sol1, sol2, aa, bb, cc, name)
    plt.xlim([0, 7.2])
    addArrowAnnotations(fig2)

if __name__ == "__main__":
    main()