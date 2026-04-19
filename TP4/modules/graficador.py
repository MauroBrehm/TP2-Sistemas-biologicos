import matplotlib.pyplot as plt
import numpy as np


def graficar_resultados (simulacion, labels):
    tiempos = [p[0] for p in simulacion] 
    CD4 = [p[1][0] for p in simulacion]
    CD8 = [p[1][1] for p in simulacion]
    V   = [p[1][2]*10000 for p in simulacion] #copias/ml (multiplicamos para que se vea en misma escala que CD4 y CD8)

    plt.plot(tiempos, CD4, label=f'CD4 - {labels}', color='purple')
    plt.plot(tiempos, CD8, label=f'CD8 - {labels}', color='magenta')
    plt.plot(tiempos, V, label=f'Virus - {labels}', color='red')

    plt.xlabel('Tiempo (años)')
    plt.ylabel('Población')
    plt.title('Simulación de la dinámica del VIH')

    plt.ylim(0,1000)

    plt.legend()
    plt.grid(True)
    plt.show()
