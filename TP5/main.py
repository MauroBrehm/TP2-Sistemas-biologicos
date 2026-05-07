import numpy as np
from scipy.linalg import eigvals
from modules.graficador import graficar_resultados, graficar_resultados_dist, graficar_varias_ci
from scipy.integrate import odeint #scipy es una biblioteca de Python para cálculos científicos,
#y odeint es una función que resuelve sistemas de EDOs.
#podesmos hacerlo con euler tambien pero es mas lento y menos preciso, es mas facil con esta funcion

#Parametros
beta_h = 1 #tasa a la que los mosquitos infectados contagian personas
gamma_h= 0.4 # tasa a la que las personas se curan
beta_m = 1 #tasa a la que las personas infectadas contagian mosquitos
gamma_m=2 #tasa a la que los mosquitos se curan o mueren

#condiciones iniciales
i0=0.1 #humanos infectados al inicio (entre 0 y 1)
a0=0.1 #mosquitos infectados al inicio (entre 0 y 1)
y0=[i0,a0]


#=====================================================================
#ejercicio 1
#=====================================================================

#sistema de ecuaciones diferenciales --> modelo SIS (Suceptibles-Infectados-Susceptibles)
def modelo_malaria (y,t,beta_h,gamma_h,beta_m,gamma_m):
    i,a=y #proporcion de humanos infectados y mosquitos infectados

    di_dt= beta_h*(1 - i)*a - gamma_h*i
    '(1-i) = s --> proporcion humanos sanos/suceptibles'
    'beta_h*(1 - i)*a --> cuantos humanos se infectan por dia'
    'gamma_h*i --> cuantos humanos se curan por dia'

    da_dt= beta_m*(1 - a)*i - gamma_m*a
    '(1-a) = v --> proporcion mosquitos sanos'
    'beta_m*(1 - a)*i --> cuantos mosquitos se infectan por dia'
    'gamma_m*a --> cuantos mosquitos se curan o mueren por dia'

    return [di_dt,da_dt]

'Interpretacion del modelo:'
'-Los humanos se infectan cuando un mosquito infectado (a) pica a un humano sano (1−i), con tasa βh'
'-Los humanos se curan con tasa γh'
'-Los mosquitos se infectan cuando pican a un humano infectado (i) y el mosquito estaba sano (1−a), con tasa βm'
'-Los mosquitos se curan con tasa γm'

#=====================================================================
#ejercicio 2
#=====================================================================
def equilibrio (beta_h,gamma_h,beta_m,gamma_m):
    'Calcula el punto de equilibrio del sistema'
    
    'El n° reproductivo básico R0 se calcula como:'
    'R0 = sqrt(beta_h * beta_m) / (gamma_h * gamma_m)'
    'R0 > 1 --> epidemia se propaga (equilibrio endémico estable)'
    'R0 <= 1 --> epidemia se extingue (equilibrio trivial estable)'
    R0 = np.sqrt((beta_h * beta_m) / (gamma_h * gamma_m))   
    resultados = {'R0': R0}
    
    #Punto de equilibrio trivial --> Jacobiano evaluado en (0,0)
    J_trivial = np.array([[-gamma_h, beta_h], 
                          [beta_m, -gamma_m]])
    autoval_trivial = eigvals(J_trivial)

    #Estable si todos los autovalores tienen parte real negativa
    estable_trivial = np.all(np.real(autoval_trivial) < 0)

    resultados['Punto_trivial'] = {'i*': 0, 'a*': 0, 'autovalores': autoval_trivial, 'estable': estable_trivial} 

    #Punto de equilibrio endémico --> se calcula resolviendo el sistema en equilibrio
    if R0 > 1:
        i_star =(beta_h * beta_m - gamma_h * gamma_m) / (beta_h * beta_m + gamma_h * beta_m)

        #una vez que tenemos i* podemos calcular a* despejando di/dt=0
        a_star = (gamma_h * i_star) / (beta_h * (1 - i_star))
        
        #calculamos Jac evaluado en (i*, a*)
        J_endemico = np.array([
            [-beta_h * a_star - gamma_h,    beta_h * (1 - i_star)],
            [beta_m * (1 - a_star),        -beta_m * i_star - gamma_m]
        ])
        autoval_endemico = eigvals(J_endemico)
        #estable si todos los autovalores tienen parte real negativa
        estable_endemico = np.all(np.real(autoval_endemico) < 0)    
        resultados['Punto_endemico'] = {'i*': i_star, 'a*': a_star, 'autovalores': autoval_endemico, 'estable': estable_endemico} 
    else:
        #no existe equilibrio endémico estable, el sistema se extingue
        resultados['endemico'] = None

    return resultados

#=====================================================================
#ejercicio 3
#=====================================================================
def simular_varias_ci(condiciones, modelo, tiempo, parametros, equilibrio_endemico):
    """Simula el modelo para varias condiciones iniciales y grafica cada caso."""
    resultados = []
    for cond in condiciones:
        y0_cond = [cond['i0'], cond['a0']]
        solucion = odeint(modelo, y0_cond, tiempo, args=parametros)
        resultados.append({
            'tiempo': tiempo,
            'solucion': solucion,
            'label':f"{cond['label']} \n (i0={cond['i0']}, a0={cond['a0']})"
        })
    graficar_varias_ci(resultados,equilibrio_endemico)
    
#=====================================================================
#ejercicio 4
#=====================================================================

#Modelo de malaria con recuperados que tiene menor tasa de reinfectarse(beta_r) y perdida de inmunidad parcial
# w es la tasa a la que los humanos recuperados pierden su inmunidad y vuelven a ser susceptibles
def modelo_malaria_recuperados (y,t,beta_h,gamma_h,beta_m,gamma_m,beta_r,w):
    i,a,r=y #proporcion de humanos infectados, mosquitos infectados y humanos recuperados

    s=1 - i - r #proporcion de humanos susceptibles
   
    di_dt= beta_h*s*a - gamma_h*i + beta_r*a*r
    
    da_dt= beta_m*(1 - a)*i - gamma_m*a
    
    dr_dt = gamma_h * i - beta_r * a * r - w * r #tasa a la que los humanos se recuperan
    return [di_dt,da_dt,dr_dt]

#=====================================================================
#ejercicio 5 - A
#=====================================================================
def modelo_malaria_plasmido (y, t, beta_h, gamma_h, beta_m, gamma_m, beta_p, gamma_p):
    '''modelo de malaria considerando un vector con plasmido inmaduro (mosquito infectado pero no contagia)'''
    i, a, p = y
    di_dt = beta_h * (1 - i)* a - gamma_h * i #humanos solo se infectan por mosquitos infectados (a), no por los con plasmido inmaduro (p)
    da_dt = beta_p *  p - gamma_m * a 
    'beta_p * p --> los mosquitos con plasmido inmaduro que maduran y se vuelven infecciosos'
    'gamma_m * a --> los mosquitos infectados que se curan o mueren'
    dp_dt = beta_m * (1 - p)* i - gamma_p * p - beta_p * p
    return [di_dt, da_dt, dp_dt]

#=====================================================================
#ejercicio 5 - B
#=====================================================================
def modelo_malaria_latencia (y, t, beta_h, gamma_h, beta_m, gamma_m, sigma):
    '''Modelo de malaria con periodo de latencia en humanos (infectados pero no contagian)'''
    i,a,e = y

    # Expuestos: susceptibles (1-i-e) picados por mosquitos infectados (a)
    # Dejan de serlo al terminar la incubación (sigma*e)
    de_dt = beta_h * (1 - i - e) * a - sigma * e 

    # Humanos infecciosos: vienen de expuestos que terminaron la incubación
    di_dt = sigma * e - gamma_h * i

    #Mosquitos infectados: se infectan al picar humanos infecciosos (i) y dejan de serlo al curarse o morir
    da_dt = beta_m * (1 - a) * i - gamma_m * a
    return [ di_dt, da_dt, de_dt ]

#=====================================================================
#ejercicio 5 - C
#=====================================================================
def mortalidad_estacional(t, gamma_m, u, periodo=1.0, duracion_humeda=0.75):
    '''Mortandad de mosquitos que aumenta después del fin de la estación húmeda
        u: aumento de mortalidad en estaciones secas
    '''
    fase = t % periodo
    if fase >= duracion_humeda:
        return gamma_m + u
    return gamma_m

def modelo_malaria_estacional(y, t, beta_h, gamma_h, beta_m, gamma_m, u, periodo=1.0, duracion_humeda=0.75):
    '''Modelo de malaria con mortalidad estacional de mosquitos'''
    i, a = y
    gamma_m_t = mortalidad_estacional(t, gamma_m, u, periodo, duracion_humeda)

    di_dt = beta_h * (1 - i) * a - gamma_h * i
    da_dt = beta_m * (1 - a) * i - gamma_m_t * a
    return [di_dt, da_dt]

#=====================================================================
#simulaciones
#=====================================================================
tiempo = np.linspace(0, 120, 1000) #tiempo de simulacion de 120 dias
simulacion = odeint(modelo_malaria, y0, tiempo, args=(beta_h, gamma_h, beta_m, gamma_m))
graficar_resultados((tiempo, simulacion), 'Simulación modelo SIS de malaria')

#Equilibrio y estabilidad
eq = equilibrio(beta_h, gamma_h, beta_m, gamma_m)
print(f"\n  R0 = {eq['R0']:.4f}")
print(f"  Equilibrio trivial   → estable: {eq['Punto_trivial']['estable']}")
if eq['Punto_endemico']:
    print(f"  Equilibrio endémico → i* = {eq['Punto_endemico']['i*']:.4f}, "
          f"a* = {eq['Punto_endemico']['a*']:.4f}, "
          f"estable: {eq['Punto_endemico']['estable']}")

#se plantean distintas condiciones iniciales que representan causas que provocan brotes de malaria
condiciones_iniciales = [
    {'label': 'Baja infección inicial', 'i0': 0.01, 'a0': 0.05},
    {'label': 'Bajo acceso sanitario', 'i0': 0.05, 'a0': 0.20},
    {'label': 'Alta densidad de mosquitos', 'i0': 0.02, 'a0': 0.40},
    {'label': 'Brotes recurrentes por movilidad poblacional', 'i0': 0.50, 'a0': 0.15}
]
simular_varias_ci(condiciones_iniciales, modelo_malaria, tiempo, (beta_h, gamma_h, beta_m, gamma_m), eq['Punto_endemico'])

#simulacion del modelo con vector con plasmido inmaduro
beta_p = 0.1 
gamma_p = 0.2 
simulacion_plasmido = odeint(modelo_malaria_plasmido, [i0, 0, a0], #i=0.1, a=0 (ninguno infeccioso), p=0.1 
                             tiempo, args=(beta_h, gamma_h, beta_m, gamma_m, beta_p, gamma_p))
graficar_resultados_dist((tiempo, simulacion_plasmido), ['Modelo con período de latencia del plasmido en el mosquito', 'Mosquitos con plasmido inmaduro p(t)'])

#simulacion del modelo con periodo latencia (humano infectado pero no contagia)
sigma = 0.3
simulacion_latencia = odeint(modelo_malaria_latencia, [i0,a0, 0], #i=0.1, a=0.1, e=0 (nadie incubando al inicio)
                              tiempo, args = (beta_h, gamma_h, beta_m, gamma_m, sigma))
graficar_resultados_dist((tiempo, simulacion_latencia), ['Modelo con período de incubación en humanos', 'Humanos expuestos e(t)'])

#simulacion del modelo considerando la mortalidad estacional de los mosquitos
u = 0.5
simulacion_estacional = odeint(modelo_malaria_estacional, y0, tiempo, args=(beta_h, gamma_h, beta_m, gamma_m, u))
graficar_resultados((tiempo, simulacion_estacional), 'Modelo con mortalidad estacional de mosquitos al final de la estación húmeda')