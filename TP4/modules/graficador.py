import matplotlib.pyplot as plt



def graficar_resultados (simulacion, labels):
    tiempos = [p[0] for p in simulacion] 
    Is= [p[1][0] for p in simulacion]
    Gs = [p[1][1] for p in simulacion]
    
    plt.plot(tiempos, Gs, label=f'Glucosa', color='fuchsia')
    plt.plot(tiempos, Is, label=f'Insulina ', color='turquoise')

    plt.xlabel('Tiempo (minutos)')
    plt.ylabel('Concentración (mg/dL para glucosa, mU/L para insulina)')
    plt.title(f'{labels}')

    plt.ylim(0,300)

    plt.legend()
    plt.grid(True)
    plt.show()

#Grafica para la inyeccion de insulina y para la ingesta de glucosa
def graficar_ingreso(tiempos, Gin, Iext, labels):
    plt.plot(tiempos, Gin, label='Glucosa ingerida', color='magenta')
    plt.plot(tiempos, Iext, label='Insulina inyectada', color='purple')

    plt.xlabel('Tiempo (minutos)')
    plt.ylabel('Cantidad (mg o mU)')
    plt.title(f'Inyección de insulina e Ingreso de Glucosa - {labels}')
    
    plt.legend()
    plt.grid(True)
    plt.show()
