import sympy as sp

def derivadas(func,orden):#Esto lo saque de Ecuaciones Diferenciales
    for sym in func.free_symbols:
        t = sym
    derivadas = []
    #derivada[0]=func
    for i in range(orden):
        derivadas.append(func.diff(t,i))
    return derivadas

# r=10.0
# K=1000.0
# p0=10.0
#funcion P(t) = K/(1+(K/p0-1)*exp(-r*t)) ya esta resuelta
t = sp.symbols('t', real=True)
def P(t):
    return 1000.0/(1+(1000.0/10.0-1)*sp.exp(-10.0*t))    
print (P(t))