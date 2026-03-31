import matplotlib.pyplot as plt
import numpy as np


def graficar_resultados(resultados: list, labels, metodo: str):
    '''Grafica los resultados de un método numérico.'''
    fig, ax = plt.subplots()
    for res, label in zip(resultados, labels):
        if res: 
            t_vals = [tupla[0] for tupla in res]
            p_vals = [tupla[1] for tupla in res]
            ax.plot(t_vals, p_vals, label=label)

    ax.set(xlabel='t = tiempo de crecimiento (años)', ylabel='P(t) = población',
        title=f'Metodo {metodo} - Distintos pasos')
    ax.set_ylim((0, 1200))
    ax.set_xlim((0, 6))
    ax.grid()
    ax.legend()

    fig.savefig(f"metodo_{metodo.lower()}.png")
    plt.show()

def graficar_errores(resultados: list, labels, metodo:str):
    fig, ax = plt.subplots()
    for res, label in zip(resultados, labels):
        if res: 
            t_vals = [tupla[0] for tupla in res]
            p_vals = [tupla[1] for tupla in res]
            ax.plot(t_vals, p_vals, label=label)

    ax.set(title=f'Curva de errores - Distintos pasos')
    ax.set_ylim((0, 700))
    ax.set_xlim((0, 6))
    ax.grid()
    ax.legend()

    fig.savefig(f"curva_errores_{metodo.lower()}.png")
    plt.show()