import matplotlib.pyplot as plt
import numpy as np


def grafico_euler(resultados:list, labels):
    fig, ax = plt.subplots()
    for res, label in zip(resultados, labels):
        if res: 
            t_vals = [tupla[0] for tupla in res]
            p_vals = [tupla[1] for tupla in res]
            ax.plot(t_vals, p_vals, label= label)

    ax.set(xlabel='t = tiempo de crecimiento', ylabel='P(t) = población',
        title='Metodo Euler - Distintos pasos')
    ax.set_ylim((0,1200))
    ax.set_xlim((0,6))
    ax.grid()
    ax.legend()

    fig.savefig("metodoEuler.png")
    plt.show()


def grafico_taylor(resultados:list, labels):
    fig, ax = plt.subplots()
    for res, label in zip(resultados, labels):
        if res: 
            t_vals = [tupla[0] for tupla in res]
            p_vals = [tupla[1] for tupla in res]
            ax.plot(t_vals, p_vals, label= label)

    ax.set(xlabel='t = tiempo de crecimiento', ylabel='P(t) = población',
        title='Metodo Taylor - Distintos pasos')
    ax.set_ylim((0,1200))
    ax.set_xlim((0,6))
    ax.grid()
    ax.legend()

    fig.savefig("metodoTaylor.png")
    plt.show()