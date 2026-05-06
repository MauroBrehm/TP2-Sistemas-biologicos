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
