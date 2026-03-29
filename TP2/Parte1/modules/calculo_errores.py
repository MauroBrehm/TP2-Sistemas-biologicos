
def calcular_error(sol_numerica: list, sol_exacta_func):
    """Calcula el error absoluto solución numérica y exacta"""
    errores = []
    for i in sol_numerica:
        for t, P_num in i:
            P_exact = sol_exacta_func(t)
            error_abs = abs(P_num - P_exact)
            errores.append((t, error_abs))
    
    return errores #no siento q este bien
