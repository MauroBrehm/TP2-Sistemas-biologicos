import matplotlib.pyplot as plt

def graficar_resultados(datos, labels):
    tiempo, solucion = datos
    hs = solucion[:, 0]  # Humanos infectados
    ms = solucion[:, 1]  # Mosquitos infectados
    
    plt.plot(tiempo, hs, label=f'Humanos infectados', color='steelblue')
    plt.plot(tiempo, ms, label=f'Mosquitos infectados', color='deeppink')

    plt.xlabel('Tiempo (años)')
    plt.ylabel('Proporción')
    plt.title(f'{labels}')

    plt.legend()
    plt.grid(True)
    plt.show()

def graficar_resultados_dist ( simulacion, labels):
    tiempo, solucion = simulacion
    i = solucion[:, 0]  # Humanos infectados
    a = solucion[:, 1]  # Mosquitos infectados
    r = solucion[:, 2]  

    plt.plot(tiempo, i, label=f'Humanos infectados', color='steelblue')
    plt.plot(tiempo, a, label=f'Mosquitos infectados', color='deeppink')
    plt.plot(tiempo, r, label=f'{labels[1]}', color='lightgreen')

    plt.xlabel('Tiempo (años)')
    plt.ylabel('Proporción')   
    plt.title(f'{labels[0]}')

    plt.legend()
    plt.grid(True)
    plt.show()