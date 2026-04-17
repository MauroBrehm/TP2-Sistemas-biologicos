def obtener_parametros():
    #los valores los saque de la parte donde dice simulacion por si las dudas
    parametros = {"a": 0.25, "b": 50,"c": 0.25, "d": 10, "e":0.01, "f":0.0045,
                  "CD4N": 1000, "CD8N": 550
                  }
    #g = e . CD4N 
    parametros ["g"]= parametros["e"] * parametros["CD4N"] 

    # Validar que la relación CD4N/CD8N esté dentro del rango 1.2 a 2.2
    relacion = parametros["CD4N"] / parametros["CD8N"]
    if not (1.2 <= relacion <= 2.2):
        raise ValueError(f"La relación CD4N/CD8N = {relacion:.2f} no está en el rango permitido [1.2, 2.2]. ")

    return parametros
    