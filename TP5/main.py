import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
from scipy.integrate import odeint #scipy es una biblioteca de Python para cálculos científicos,
#y odeint es una función que resuelve sistemas de EDOs.
#podesmos hacerlo con euler tambien pero es mas lento y menos preciso, es mas facil con esta funcion

#Parametros
beta_h = 1 #tasa a la que los mosquitos infectados contagian personas
gamma_h=0.4 # tasa a la que las personas se curan
beta_m = 1 #tasa a la que las personas infectadas contagian mosquitos
gamma_m=2 #tasa a la que los mosquitos se curan o mueren

#sistema de ecuaciones diferenciales
def modelo_malaria (y,t,beta_h,gamma_h,beta_m,gamma_m):
    i,a=y #proporcion de humanos infectados y mosquitos infectados

    di_dt= beta_h*(1 - i)*a - gamma_h*i
    '(1-i) = s --> proporcion humanos sanos/suceptibles'

    da_dt= beta_m*(1 - a)*i - gamma_m*a
    '(1-a) = v --> proporcion mosquitos sanos'
    return [di_dt,da_dt]

'Interpretacion del modelo:'
'-Los humanos se infectan cuando un mosquito infectado (a) pica a un humano sano (1−i), con tasa βh'
'-Los humanos se curan con tasa γh'
'-Los mosquitos se infectan cuando pican a un humano infectado (i) y el mosquito estaba sano (1−a), con tasa βm'
'-Los mosquitos se curan con tasa γm'

#condiciones iniciales
i0=0.1 #humanos infectados (entre 0 y 1)
a0=0.1 #mosquitos infectados (entre 0 y 1)
y0=[i0,a0]

def equilibrio (beta_h,gamma_h,beta_m,gamma_m):
    'Calcula el punto de equilibrio del sistema'
    '-Punto 1 (trivial): i*=0, a*=0 --> existe siempre pero solo es estable si R0 <= 1'
    '-Punto 2 (endémico): i*>0, a*>0 --> existe solo si R0 > 1'
    
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
    resultados['Punto_trivial'] = {'i*': 0, 'a*': 0, 'autovalores': autoval_trivial} #creo que deberiamos ver si es estable o no

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
        resultados['Punto_endemico'] = {'i*': i_star, 'a*': a_star, 'autovalores': autoval_endemico} #lo mismo que con el otro creo que deberiamos ver si es estabvle o no pero no se como
    else:
        #no existe equilibrio endémico estable, el sistema se extingue
        resultados['endemico'] = None

    return resultados