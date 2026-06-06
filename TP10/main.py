import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation# Para la animación

 
# ─────────────────────────────────────────────────────────────────────────────
# PARÁMETROS DE LA SIMULACIÓN
# ─────────────────────────────────────────────────────────────────────────────
T  = 0.001   # Periodo de integración [s]: cuamto tiempo avanza cada iteración
 
# Pesos que regulan cuánto influye cada regla en la velocidad final
WA = 0.1     # Peso de alineación  
WS = 0.1     # Peso de separación  
WC = 0.1     # Peso de cohesión    
 
# Radios de vecindad: distancia maxima para considerar a otro boid como vecino
RHO_A = 4.0  # Radio de alineación  [m]
RHO_S = 0.5  # Radio de separación  [m] --> tiene que ser menor que RHO_C para que haya agrupamiento
RHO_C = 4.0  # Radio de cohesión    [m]
 
# Velocidad inicial tipica de un ave ( aprox 5 m/s, valor sugerido en el enunciado)
V0 = 5.0     # [m/s]
 
# Límite de velocidad máxima (evita que los boids se disparen al infinito)
V_MAX = 15.0  # [m/s]
 
# Dominio espacial observable [m] --> region 2D donde se mueven los boids
X_MIN, X_MAX = 0.0, 50.0
Y_MIN, Y_MAX = 0.0, 50.0
historial=[] # Para guardar las posiciones de los boids en cada paso y luego graficar la trayectoria

class Boid:
    def __init__(self, pos, vel):
        self.pos = pos  # Posición del boid (vector 2D)
        self.vel = vel  # Velocidad del boid (vector 2D)

# ─────────────────────────────────────────────────────────────────────────────
# FUNCIONES DE LAS REGLAS DE COMPORTAMIENTO
# ─────────────────────────────────────────────────────────────────────────────
    def calcular_alineacion(self, todos):
        """
    Regla de alineación --> a_i = (1 / #A_i) * Σ_{j ∈ A_i} v_j  −  v_i
 
    Empuja al boid para que iguale su velocidad con el promedio de sus vecinos
    dentro del radio RHO_A
    Retorna el vector de influencia (0 si no hay vecinos)
    """
        vecinos_vel = []
 
        for j, otro in enumerate(todos):
            if otro is self:
                continue  # Un boid no es vecino de si mismo
 
            distancia = np.linalg.norm(self.pos - otro.pos)  # Distancia euclidea [m]
 
            if distancia < RHO_A:  # El boid j esta dentro de la vecindad de alineación
                vecinos_vel.append(otro.vel)
 
        if len(vecinos_vel) == 0:
            return np.zeros(2)  # Sin vecinos --> influencia nula (evita división por 0)
 
        velocidad_promedio = np.mean(vecinos_vel, axis=0)  # Velocidad media de los vecinos
        return velocidad_promedio - self.vel                # Diferencia respecto a la propia
 
    def calcular_separacion(self, todos):
        """
        Regla de separación --> s_i = Σ_{j ∈ S_i} (x_i − x_j)
    
        Aleja al boid de sus vecinos muy cercanos (radio RHO_S) para evitar colisiones
        La dirección del vector apunta DESDE el vecino HACIA el boid --> lo empuja a alejarse
        Retorna el vector de influencia (0 si no hay vecinos)
        """
        influencia = np.zeros(2)
    
        for otro in todos:
            if otro is self:
                continue
    
            distancia = np.linalg.norm(self.pos - otro.pos)
    
            if distancia < RHO_S and distancia > 0:  # Vecino demasiado cerca
                influencia += self.pos - otro.pos     # Vector que apunta hacia afuera del vecino
    
        return influencia
    
    def calcular_cohesion(self, todos):
        """
        Regla de cohesión --> c_i = (1 / #C_i) * Σ_{j ∈ C_i} x_j  −  x_i
    
        Atrae al boid hacia el centro geométrico de su vecindad (radio RHO_C)
        Opuesta a la separación: en lugar de alejarse, se acerca al grupo
        Retorna el vector de influencia (0 si no hay vecinos)
        """
        vecinos_pos = []
    
        for otro in todos:
            if otro is self:
                continue
    
            distancia = np.linalg.norm(self.pos - otro.pos)
    
            if distancia < RHO_C:  # Vecino dentro de la vecindad de cohesión
                vecinos_pos.append(otro.pos)
    
        if len(vecinos_pos) == 0:
            return np.zeros(2)  # Sin vecinos --> influencia nula
    
        centro = np.mean(vecinos_pos, axis=0)  # Centro geométrico del grupo
        return centro - self.pos               # Vector que apunta hacia ese centro

    def calcular_fuerza_borde(self):
        """
        Fuerza de retorno cuando un boid sale del área observable (condición de borde)
    
        Segun el TP --> si el boid cruza el limite izquierdo (x < X_MIN),
        se aplica f = [X_MIN - x_i, 0]; si cruza el derecho (x > X_MAX),
        f = [X_MAX - x_i, 0], si cruza el borde inferior (y < Y_MIN), f = [0, Y_MIN - y_i], 
        y si cruza el borde superior (y > Y_MAX), f = [0, Y_MAX - y_i]
        Esta fuerza se suma a la ec 4 como un termino adicional
        """
        fuerza = np.zeros(2)
    
        # Borde izquierdo
        if self.pos[0] < X_MIN:
            fuerza[0] += X_MIN - self.pos[0]
    
        # Borde derecho
        if self.pos[0] > X_MAX:
            fuerza[0] += X_MAX - self.pos[0]
    
        # Borde inferior
        if self.pos[1] < Y_MIN:
            fuerza[1] += Y_MIN - self.pos[1]
    
        # Borde superior
        if self.pos[1] > Y_MAX:
            fuerza[1] += Y_MAX - self.pos[1]
    
        return fuerza

    def actualizar(self,todos):
        # Calcular las contribuciones de cada regla
        a = self.calcular_alineacion(todos)
        s = self.calcular_separacion(todos)
        c = self.calcular_cohesion(todos)
        f = self.calcular_fuerza_borde()

        # Ecuación (4)
        self.vel = self.vel + WA*a + WS*s + WC*c + f

        # Limitar velocidad máxima
        velocidad = np.linalg.norm(self.vel)

        if velocidad > V_MAX:
            self.vel = (self.vel / velocidad) * V_MAX

        # Actualizar posición
        self.pos = self.pos + T*self.vel

#Ahora ya podemos hacer los casos 
#Creo 2 pajaritos 
boids = [
    Boid(
        np.array([10.0, 25.0]),
        np.array([5.0, 0.0])
    ),

    Boid(
        np.array([40.0, 25.0]),
        np.array([-5.0, 0.0])
    )
]
#La simulacion se ejecuta durante un número determinado de pasos,  
# en cada paso se actualizan las posiciones y velocidades de los boids según las reglas de comportamiento.
#  Para evitar que los boids se actualicen secuencialmente (lo que podría introducir sesgos),
#  se hace una copia temporal del estado actual antes de actualizar cada boid.
PASOS = 5000

for paso in range(PASOS):

    # Guardar estado actual
    posiciones = [b.pos.copy() for b in boids]
    velocidades = [b.vel.copy() for b in boids]

    # Crear una copia temporal
    copia = []

    for pos, vel in zip(posiciones, velocidades):
        copia.append(
            Boid(pos.copy(), vel.copy())
        )

    # Actualizar cada boid usando la copia
    for i, boid in enumerate(boids):
        boid.actualizar(copia)
    # Guardar posiciones para graficar trayectoria
    historial.append(np.array([b.pos.copy() for b in boids]))
