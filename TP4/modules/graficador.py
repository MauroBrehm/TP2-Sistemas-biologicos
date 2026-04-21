import matplotlib.pyplot as plt



def graficar_resultados (simulacion, labels):
    tiempos = [p[0] for p in simulacion] 
    Gs= [p[1][0] for p in simulacion]
    Is = [p[1][1] for p in simulacion]
    
    plt.plot(tiempos, Gs, label=f'Glucosa - {labels}', color='blue')
    plt.plot(tiempos, Is, label=f'Insulina - {labels}', color='orange')

    plt.xlabel('Tiempo (minutos)')
    plt.ylabel('Concentración (mg/dL para glucosa, mU/L para insulina)')
    plt.title('Simulación de Glucosa e Insulina en Sangre')

    plt.ylim(0,1000)

    plt.legend()
    plt.grid(True)
    plt.show()

#Grafica para la inyeccion de insulina y para la ingesta de glucosa
def graficar_ingreso(tiempos, Gin, Iext):
    plt.plot(tiempos, Gin, label='Glucosa ingerida', color='green')
    plt.plot(tiempos, Iext, label='Insulina inyectada', color='red')

    plt.xlabel('Tiempo (minutos)')
    plt.ylabel('Cantidad (mg o mU)')
    plt.title('Ingreso de Glucosa e Insulina')

    plt.legend()
    plt.grid(True)
    plt.show()