import matplotlib.pyplot as plt
import numpy as np


def plot_tumor_grid(grid, day):
    #grafica estado del tumor.
    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(grid, cmap='magma', vmin=0, vmax=1)
    ax.set_title(f"Concentración Tumoral | Día {int(day)}", fontsize=14)
    ax.axis('off')

    cbar = plt.colorbar(im, ax=ax, shrink=0.8, label='Concentración Celular (0 a 1)')
    cbar.ax.tick_params(labelsize=10)

    plt.tight_layout()
    return fig


def plot_growth_curve(history, dt):
    fig, ax = plt.subplots(figsize=(8, 4))

    total_steps = len(history)
    days_axis = np.linspace(0, total_steps * dt, num=total_steps)

    ax.plot(days_axis, history, linewidth=2, color='#007BFF')
    ax.set_title("Evolución de la Masa Tumoral Total", fontsize=16)
    ax.set_xlabel("Días", fontsize=12)
    ax.set_ylabel("Masa Tumoral Total (u.a.)", fontsize=12)

    ax.grid(True, linestyle='--', alpha=0.6)
    ax.tick_params(labelsize=10)

    plt.tight_layout()
    return fig