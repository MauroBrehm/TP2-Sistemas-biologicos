import matplotlib.pyplot as plt

def graficar_resultados (simulacion):
    '''Grafica modelo'''
    tiempos = [p[0] for p in simulacion]
    G6P = [p[1][0] for p in simulacion]
    F6P = [p[1][1] for p in simulacion]
    FBP = [p[1][2] for p in simulacion]
    ATP = [p[1][3] for p in simulacion]
    ADP = [p[1][4] for p in simulacion]
    AMP = [p[1][5] for p in simulacion]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(tiempos, G6P, label=f'G6P', color='steelblue')
    ax.plot(tiempos, F6P, label=f'F6P', color='deeppink')
    ax.plot(tiempos, FBP, label=f'FBP', color='lightgreen')
    ax.plot(tiempos, ATP, label=f'ATP', color='purple')
    ax.plot(tiempos, ADP, label=f'ADP', color='orange')
    ax.plot(tiempos, AMP, label=f'AMP', color='cyan')

    ax.set_ylim(0, 10)
    ax.set_xlabel('Tiempo (minutos)')
    ax.set_ylabel('Concentración (mM)')   
    ax.set_title(f'Modelo dinamico de glucolisis')

    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

