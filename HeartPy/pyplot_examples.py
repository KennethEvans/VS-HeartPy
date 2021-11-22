import numpy as np
import matplotlib.pyplot as plt

def show_markers():
    '''show_markers
    Shows marker possibilities for matplotlib.pyplot
    '''
    markers = [0, 1, 2, 3, 4, 5, 6, 7, 'o', 'h', '_', '1', '2', '3', '4',
          '8', 'p', '^', 'v', '<', '>', '|', 'd', ',', '+', 's', '*',
          '|', 'x', 'D', 'H', '.']
    n_markers = len(markers)
    size = 20 * n_markers, 300
    dpi = 72.0
    #figsize= size[0] / float(dpi), size[1] / float(dpi)
    #fig = plt.figure(figsize=figsize, dpi=dpi)
    plt.figure(figsize=(8,6))
    plt.axes([0, 0.01, 1, .9], frameon=False)
    for i, m in enumerate(markers):
        X = i * .5 * np.ones(11)
        Y = np.arange(11)
        plt.plot(X, Y, lw=1, marker=m, ms=10, mfc=(.75, .75, 1, 1),
                mec=(0, 0, 1, 1))
        plt.text(.5 * i, 10.25, repr(m), rotation=0, fontsize=10, va='bottom')
    #plt.xlim(-.2, .2 + .5 * n_markers)
    plt.xticks([])
    plt.yticks([])
    plt.title('Marker Possibilities\n')
    plt.show()

def axes_demo():
    np.random.seed(19680801)  # Fixing random state for reproducibility.

    # create some data to use for the plot
    dt = 0.001
    t = np.arange(0.0, 10.0, dt)
    r = np.exp(-t[:1000] / 0.05)  # impulse response
    x = np.random.randn(len(t))
    s = np.convolve(x, r)[:len(x)] * dt  # colored noise

    fig, main_ax = plt.subplots(figsize=(8, 6))
    main_ax.plot(t, s)
    main_ax.set_xlim(0, 1)
    main_ax.set_ylim(1.1 * np.min(s), 2 * np.max(s))
    main_ax.set_xlabel('time (s)')
    main_ax.set_ylabel('current (nA)')
    main_ax.set_title('Gaussian colored noise')

    # this is an inset axes over the main axes
    right_inset_ax = fig.add_axes([.65, .6, .2, .2], facecolor='k')
    right_inset_ax.hist(s, 400, density=True)
    right_inset_ax.set(title='Probability', xticks=[], yticks=[])

    # this is another inset axes over the main axes
    left_inset_ax = fig.add_axes([.2, .6, .2, .2], facecolor='k')
    left_inset_ax.plot(t[:len(r)], r)
    left_inset_ax.set(title='Impulse response', xlim=(0, .2), xticks=[], yticks=[])

    plt.show()

def main():
    show_markers()
    #axes_demo()

if __name__ == "__main__":
    main()

