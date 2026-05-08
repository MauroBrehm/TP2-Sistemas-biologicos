import matplotlib.pyplot as plt

def graficar_resultados(datos, labels):
    '''Grafica i(t) y a(t)'''
    tiempo, solucion = datos
    i = solucion[:, 0]  # Humanos infectados
    a = solucion[:, 1]  # Mosquitos infectados
    s = 1 - i           # humanos sanos/suceptibles (s = 1 - i)
    v = 1 - a           # mosquitos sanos/ suceptibles (v = 1 - a)
    
    fig, ax = plt.subplots(figsize=(9, 5))

    ax.plot(tiempo, i, label=f'Humanos infectados i(t)', color='steelblue')
    ax.plot(tiempo, a, label=f'Mosquitos infectados a(t)', color='deeppink')
    ax.plot(tiempo, s, label=f'Humanos sanos s(t)', color='lightgreen')
    ax.plot(tiempo, v, label=f'Mosquitos sanos v(t)', color='orange' )

    ax.set_ylim(-0.2, 1.2) #tiene que ir de 0 a 1 porque son proporciones de la población total
    ax.set_xlabel('Tiempo (dias)')
    ax.set_ylabel('Proporción de la población')
    ax.set_title(f'{labels}')

    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

def graficar_resultados_dist ( simulacion, labels):
    '''Grafica modelos con 3 compartimentos (i(t), a(t) y un tercero: e,r o p)'''
    tiempo, solucion = simulacion
    i = solucion[:, 0]  # Humanos infectados
    a = solucion[:, 1]  # Mosquitos infectados
    r = solucion[:, 2]  #Tercer compartimento (expuestos, recuperados o plasmidos)
    s=1-i+r # humanos sanos 
    v=1-a # mosquitos sanos 

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(tiempo, i, label=f'Humanos infectados i(t)', color='steelblue')
    ax.plot(tiempo, a, label=f'Mosquitos infectados a(t)', color='deeppink')
    ax.plot(tiempo, r, label=f'{labels[1]}', color='lightgreen')
    ax.plot(tiempo, s, label=f'Humanos sanos s(t)', color='purple')
    ax.plot(tiempo, v, label=f'Mosquitos sanos v(t)', color='orange')

    ax.set_ylim(-0.2, 1.2) #tiene que ir de 0 a 1 porque son proporciones de la población total
    ax.set_xlabel('Tiempo (dias)')
    ax.set_ylabel('Proporción de la población')   
    ax.set_title(f'{labels[0]}')

    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

def graficar_varias_ci(resultados, equilibrio):
    """Grafica múltiples condiciones iniciales en una sola figura con subplots."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle('Ejercicio 3 — Distintas condiciones iniciales (modelo SIS)', fontsize=13)
 
    for ax, res in zip(axes.flat, resultados):
        t = res['tiempo']
        sol = res['solucion']
        i = sol[:, 0]
        a = sol[:, 1]
        s = 1 - i   # humanos sanos
        v = 1 - a   # mosquitos sanos
 
        ax.plot(t, i, color='steelblue', label='i(t) humanos')
        ax.plot(t, a, color='deeppink', label='a(t) mosquitos')
        ax.plot(t, s, color='lightgreen', label='s(t) humanos')
        ax.plot(t, v, color='orange', label='v(t) mosquitos')

        ax.set_title(res['label'], fontsize=10)
        ax.set_xlabel('Tiempo (días)')
        ax.set_ylabel('Proporción')
        ax.set_ylim(-0.2, 1.2)
        ax.legend()
        ax.grid(True)
 
    plt.tight_layout()
    plt.show()