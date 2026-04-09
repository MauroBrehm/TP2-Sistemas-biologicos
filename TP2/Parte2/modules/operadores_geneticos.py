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
    valor_decimal = round((p - p_min) / (p_max - p_min) * (2**L - 1))
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