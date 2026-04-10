def decode_chromosome(chromosome, p_min, p_max):
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

def encode_beta_gamma(beta, gamma, beta_min, beta_max, gamma_min, gamma_max, L=8):
    bits_beta = encode_chromosome(beta, beta_min, beta_max, L)
    bits_gamma = encode_chromosome(gamma, gamma_min, gamma_max, L)
    return bits_beta + bits_gamma

def calcular_ecm(I_simulado, I_observado):
    """Calcula el Error Cuadrático Medio entre valores simulados y observados"""
    n = len(I_simulado)
    ecm = sum((I_simulado[i] - I_observado[i])**2 for i in range(n)) / n # Ecuacion EMC
    return ecm

import random
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