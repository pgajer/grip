"""Render the full analytic saddle with two explicitly different axis scalings."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def render(out):
    r = 32
    rho, theta = np.meshgrid(np.linspace(0, 1, 25), np.linspace(0, 2*np.pi, 81), indexing='ij')
    u, v = rho*np.cos(theta), rho*np.sin(theta)
    h = u*u-v*v
    fig = plt.figure(figsize=(12, 7.2))
    cases = [(u/r, v/r, 'Equal spatial units\nAll coordinates divided by r²', ('x/r²', 'y/r²')),
             (u, v, 'Separate horizontal and vertical scaling\nHorizontal ÷ r; vertical ÷ r²', ('x/r', 'y/r'))]
    for k, (x, y, title, labels) in enumerate(cases):
        ax = fig.add_axes([.025 + .50*k, .10, .42, .70], projection='3d')
        ax.plot_surface(x, y, h, cmap='coolwarm', vmin=-1, vmax=1, alpha=.75, linewidth=0)
        ax.plot_wireframe(x, y, h, rstride=4, cstride=8, color='#334155', linewidth=.5, alpha=.5)
        ax.plot(x[-1], y[-1], h[-1], color='#172554', linewidth=1.5)
        ax.set(xlim=(-1.05, 1.05), ylim=(-1.05, 1.05), zlim=(-1.05, 1.05),
               xlabel=labels[0], ylabel=labels[1], zlabel='z/r²', title=title,
               xticks=[-1, -.5, 0, .5, 1], yticks=[-1, -.5, 0, .5, 1], zticks=[-1, -.5, 0, .5, 1])
        ax.set_box_aspect((1, 1, 1))
        ax.set_title(title, fontsize=12, pad=20)
        ax.set_zlabel('z/r²', fontsize=10, labelpad=5)
        ax.view_init(elev=22, azim=-52)
    fig.suptitle('The same full saddle z = x² − y² over the disk of radius 32\nThe dark boundary curve is included in both panels', fontsize=14, y=.97)
    fig.savefig(out)
    plt.close(fig)
