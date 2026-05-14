import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


# =============================================================================
# PUNTO 1 – Matriz estequiométrica N y calculo de Kernel
# =============================================================================

# La matriz estequiométrica N (6×8) codifica quien participa en cada reacción:
#   N[i, j] = +1  → el metabolito i es producido en la reacción j
#   N[i, j] = -1  → el metabolito i es consumido en la reacción j
#   N[i, j] =  0  → el metabolito i no participa en la reacción j

# Filas   -->  metabolitos: G6P, F6P, FBP, ATP, ADP, AMP
# Columnas --> reacciones : v1 (HK), v2 (G6Pasa), v3 (PGI), v4 (PFK), v5 (Aldolasa), v6 (ATP), v7 (ADP), v8 (AK)

# Leemos la matriz directamente del esquema del modelo:
#   d[G6P]/dt = v1 - v2 - v3        
#   d[F6P]/dt = v3 - v4            
#   d[FBP]/dt = v4 - v5            
#   d[ATP]/dt = -v1 - v2 - v4 + v6 - v7 - v8
#   d[ADP]/dt = v1 + v2 + v4 - v6 + v7 + 2*v8  (el 2 viene de AK: 2ADP)
#   d[AMP]/dt = -v8


#    v1   v2   v3   v4   v5   v6   v7   v8
N = np.array([
    [ 1,  -1,  -1,   0,   0,   0,   0,   0],   # G6P
    [ 0,   0,   1,  -1,   0,   0,   0,   0],   # F6P
    [ 0,   0,   0,   1,  -1,   0,   0,   0],   # FBP
    [-1,  -1,   0,  -1,   0,   1,  -1,  -1],   # ATP
    [ 1,   1,   0,   1,   0,  -1,   1,   2],   # ADP  (coef 2 en v8)
    [ 0,   0,   0,   0,   0,   0,   0,  -1],   # AMP
    ],dtype=float)

print("Matriz estequiométrica N:")
print(N)

# Si sumamos las filas de ATP, ADP y AMP en cada columna, deberíamos obtener 0
# Eso significa que ATP + ADP + AMP = constante --> coservacion
suma_adenilatos = N[3] + N[4] + N[5]  # fila ATP + fila ADP + fila AMP
print(f"\nSuma filas ATP+ADP+AMP (debe ser todo ceros): {suma_adenilatos}")
print("Conservación del pool de adenilatos verficada" if np.all(suma_adenilatos == 0)
      else "error en conservación")

#Calculo del kernel de N (espacio de soluciones de N*v=0)
#Kernel --> representan estados estacionarios posibles
# Usamos sympy para calcular el kernel de N
N_sym = sp.Matrix(N)
kernel = N_sym.nullspace()
print(f"\nKernel de N (estados estacionarios posibles):")
for i, vec in enumerate(kernel):
    print(f"  Vector {i+1}: {vec}")

# =============================================================================
# PUNTO 2 – Implementacion del modelo
# =============================================================================
#Parametros del modelo 
Glucosa=12.8174 #Concentración externa fija
Vmax1=1398.00 #Hexoquinasa
K_ATP1=0.10 #Hexoquinasa
K_Gluc1=0.37 #Hexoquinasa
k2=2.26 #Consumo de G6P
Vf_max3=140.282 #PGI
Vr_max3=140.282 #PGI
K_G6P3=0.80 #PGI
K_F6P3=0.15 #PGI
Vmax4=44.7287 #PFK
K_F6P4=0.021 #PFK
kappa=0.15 #PFK
k5=6.04662 #Aldolasa
k6=68.48 #Produccion de ATP
k7=3.21 #Consumo de ATP
k8f=432.9 #Adenilato Quinasa
k8r=133.33 #Adenilato Quinasa

#Condiciones iniciales sugeridas 
G6P_0=1.0
F6P_0=0.0
FBP_0=0.0
ATP_0=2.1
ADP_0=1.4
AMP_0=0.1

# Definimos las funciones de velocidad para cada reacción(Creo que deberia pasar tambien los adenilatos)
def Leyes_velocidad(t, y,Glucosa, Vmax1, K_ATP1, K_Gluc1, k2, Vf_max3, Vr_max3, K_G6P3, K_F6P3, Vmax4, K_F6P4, kappa, k5, k6, k7, k8f, k8r):
    G6P, F6P, FBP, ATP, ADP, AMP = y
    
    # v1: Hexoquinasa
    v1 = (Vmax1*ATP*Glucosa)/(1+ATP/K_ATP1+Glucosa/K_Gluc1+
        (ATP*Glucosa)/(K_ATP1*K_Gluc1))
    # v2: Consumo de G6P
    v2 = k2 *ATP* G6P
    
    # v3: PGI
    v3=(((Vf_max3/K_G6P3)*G6P)-((Vr_max3/K_F6P3)*F6P))/(1+G6P/K_G6P3+F6P/K_F6P3)
    # v4: PFK
    v4 = (Vmax4*F6P**2)
    # v5: Aldolasa
    v5 = k5 * FBP
    
    # v6: Producción de ATP
    v6 = k6 * ADP
    
    # v7: Consumo de ATP
    v7 = k7 * ATP
    
    # v8: Adenilato Quinasa
    v8=k8f*ATP*AMP-k8r*ADP**2
    
    return [0,v1, v2, v3, v4, v5, v6, v7, v8]#Le agrego un 0 asi es mas facil para indexar

    def ecuaciones_diferenciales(t, y):
        v = Leyes_velocidad(t, y,Glucosa, Vmax1, K_ATP1, K_Gluc1, k2, Vf_max3, Vr_max3, K_G6P3, K_F6P3, Vmax4, K_F6P4, kappa, k5, k6, k7, k8f, k8r)
        G6P_dt = v[1] - v[2] - v[3]
        F6P_dt = v[3] - v[4]
        FBP_dt = v[4] - v[5]
        ATP_dt = -v[1] - v[2] - v[4] + v[6] - v[7] - v[8]
        ADP_dt = v[1] + v[2] + v[4] - v[6] + v[7] + 2*v[8]
        AMP_dt = -v[8]
        return [G6P_dt, F6P_dt, FBP_dt, ATP_dt, ADP_dt, AMP_dt]