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