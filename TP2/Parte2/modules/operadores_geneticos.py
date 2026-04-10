from modules.metodo import metod_taylor_segundo_orden, crear_modelo_sir
import random

def decode_chromosome(chromosome, p_min, p_max): #convierte un cromosoma binario en un valor real entre p_min y p_max
    L = len(chromosome)
    decimal_value = sum(bit * (2 ** (L - 1 - i)) for i, bit in enumerate(chromosome))
    p = p_min + (p_max - p_min) * decimal_value / (2**L - 1)
    return p

def decode_beta_gamma(chromosome, beta_min, beta_max, gamma_min, gamma_max):
    mitad = len(chromosome) // 2
    bits_beta = chromosome[:mitad]
    bits_gamma = chromosome[mitad:]
    beta = decode_chromosome(bits_beta, beta_min, beta_max)
    gamma = decode_chromosome(bits_gamma, gamma_min, gamma_max)
    return beta, gamma

def encode_chromosome(p, p_min, p_max, L):
    if p < p_min:
        p = p_min
    if p > p_max:
        p = p_max
    #Convierte valor real a un entero entre 0 y 2^L - 1
    valor_decimal = round((p - p_min) / (p_max - p_min) * (2**L - 1))
    #Convierte a binario y rellena con 0 a la izquierda
    binario_str = bin(valor_decimal)[2:]
    binario = '0'*(L - len(binario_str)) + binario_str
    chromosome = [int(bit) for bit in binario]
    return chromosome

def encode_beta_gamma(beta, gamma, beta_min, beta_max, gamma_min, gamma_max, L_total=16): #lo hago de 16 y despues lo divido en 8 para cada parametro para que no genere confusion
    mitad = L_total // 2
    bits_beta = encode_chromosome(beta, beta_min, beta_max, mitad)
    bits_gamma = encode_chromosome(gamma, gamma_min, gamma_max, mitad)
    return bits_beta + bits_gamma

def calcular_ecm(I_simulado, I_observado):
    """Calcula el Error Cuadrático Medio entre valores simulados y observados"""
    n = len(I_simulado)
    ecm = sum((I_simulado[i] - I_observado[i])**2 for i in range(n)) / n # Ecuacion ECM
    return ecm


def seleccion_ruleta(poblacion, fitnesses, k):
    """Selecciona k individuos de la población usando selección por ruleta
    poblacion: lista de cromosomas (listas de bits)
    fitnesses: lista de valores de fitness correspondientes a cada individuo
    k: número de individuos a seleccionar
    Retorna una lista de k individuos seleccionados
    """
    #Calcular la suma total de los fitnesses
    total_fitness = sum(fitnesses)
    if total_fitness == 0:
        #Si el total es cero, seleccionar aleatoriamente
        return [ind[:] for ind in random.choices(poblacion, k=k)]
 
    #Calcular probabilidades de selección
    probabilidades = [f / total_fitness for f in fitnesses]
    
    #Seleccionar k individuos
    seleccionados = random.choices(poblacion, weights=probabilidades, k=k)
    
    return [ind[:] for ind in seleccionados] #es una lista nueva con copias de los individuos seleccionados

def cruzamiento(padre1, padre2, pc=0.8):
    """Realiza un cruce de un punto con probabilidad pc
    Si no hay cruzamiento, devuelve copias de los padres
    Retorna dos hijos resultantes del cruce
    padre1, padre2: listas de bits (cromosomas)
    """
    if random.random() < pc:
        punto= random.randint(1, len(padre1) - 1) #punto de cruce aleatorio
        hijo1 = padre1[:punto] + padre2[punto:]
        hijo2 = padre2[:punto] + padre1[punto:]
    else:
        hijo1=padre1[:]
        hijo2=padre2[:]
    return hijo1, hijo2

def mutacion(cromosoma, pm=0.02): #tipicamente pm= 1/L (es aprox 0.062 para L=16) pero se usa un valor mas bajo para evitar  aleatoriedad
    """Aplica mutacion bit a bit con probabilidad pm
    Retorna un nuevo cromosoma mutado
    """
    mutado=cromosoma[:]
    for i in range(len(mutado)):
        if random.random() < pm:
            mutado[i] = 1 - mutado[i]  #inivierte el bit
    return mutado

def evaluar_poblacion(poblacion, beta_min, beta_max, gamma_min, gamma_max, a, b, S0, I0, R0, h, indices, I_obs):
    '''Evalua que tan bueno es cada individuo de la poblacion'''
    fitnesses = []
    for cromosoma in poblacion:
        #decodificamos el cromosoma para obtener beta y gamma
        beta_val, gamma_val = decode_beta_gamma(cromosoma, beta_min, beta_max, gamma_min, gamma_max)
        #creamos el modelo SIR para esos valores
        f_individuo, jac_indiviuo = crear_modelo_sir(beta_val, gamma_val)
        #simulamos el modelo con Taylor de orden 2
        resultados_sim = metod_taylor_segundo_orden(f_individuo, jac_indiviuo, a, b, [S0, I0, R0], h)
        #sacamos curva de infectados simulados
        I_sim_totsl = [estado[1] for _, estado in resultados_sim]
        I_sim = [I_sim_totsl[i] for i in indices] #solo los 30 puntos que tenemos de observados
        #calculamos el ECM entre la curva simulada y la observada
        ecm = calcular_ecm(I_sim, I_obs)
        fitness = 1 / (1 + ecm)  #Evitamos división por cero
        fitnesses.append(fitness)
    return fitnesses