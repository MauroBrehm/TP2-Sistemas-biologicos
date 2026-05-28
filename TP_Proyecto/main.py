import numpy as np
import matplotlib.pyplot as plt
 
# =============================================================================
# DICCIONARIO DE CONFIGURACIÓN — esto lo editamos dependiendo que queremos mostrar
# =============================================================================
CONFIG = {
    # ── Condición experimental ──────────────────────────────────────────────
    # Opciones: "anaerobic" o "aerobic"
    "condicion": "anaerobic",
 
    # ── Qué modelos mostrar en la gráfica ───────────────────────────────────
    "mostrar": {
        "experimental":    True,   # puntos reales del paper
        "clasico_biomasa": True,   # FBA tradicional (maximizar biomasa)
        "paper_original":  True,   # paper con w = [1, 1, 1, 1]
        "adaptativo":      True,   # nuestra propuesta de pesos ajustables
    },
 
    # ── Variable a graficar en el eje Y ─────────────────────────────────────
    # Opciones: "crecimiento" o "etanol"
    "variable_y": "crecimiento",
 
    # ── Pesos por compartimento para el modo ADAPTATIVO ─────────────────────
    # Estos son los w_k del problema (7) del paper: max w^T * C^T * v
    # Cada peso representa cuánto contribuye ese compartimento al objetivo.
    # El script los normaliza automáticamente (suma = 1).
    #
    # Biología detrás de cada compartimento:
    #   citosol      (P1): minimizar prod. NADH → fuerza desvío a etanol
    #   mitocondria  (P2): minimizar cons. NADH/NADPH → eficiencia respiratoria
    #   peroxisoma   (P3): maximizar prod. ácidos grasos → crecimiento lento
    #   virtual      (P4): maximizar biomasa → objetivo clásico
    #
    #   Anaeróbico óptimo  → citosol=0.55, mitocondria=0.05, peroxisoma=0.10, virtual=0.50
    #   Aeróbico óptimo    → citosol=0.48, mitocondria=0.35, peroxisoma=0.08, virtual=0.45
    #   Paper original     → todos = 0.25 (equivalente a w = [1,1,1,1])
    "pesos": {
        "citosol":     0.55,
        "mitocondria": 0.05,
        "peroxisoma":  0.10,
        "virtual":     0.50,
    },
 
    # ── Optimizar pesos automáticamente ─────────────────────────────────────
    # Si True: ignora los pesos de arriba y los busca minimizando el error
    # Si False: usa exactamente los pesos que pusiste arriba
    "optimizar_pesos": False, 
}

# =============================================================================
# DATOS EXPERIMENTALES (Tablas 2 y 3 del paper)
# =============================================================================
DATOS_EXPERIMENTALES = {
    # Nissen et al. (1997)
    # Columnas: glucosa [mmol/gPS·h], etanol [mmol/gPS·h], crecimiento [1/h]
    "anaerobic": np.array([
        [ 5.56,  8.28, 0.10],
        [11.50, 17.12, 0.20],
        [17.37, 25.74, 0.30],
        [23.65, 35.27, 0.40],
    ]),
    # Heyland et al. (2009) 
    "aerobic": np.array([
        [ 7.2,  9.0, 0.16],
        [10.2, 11.2, 0.17],
        [12.2, 15.6, 0.23],
        [12.3, 15.2, 0.21],
        [15.1, 20.1, 0.33],
        [18.4, 28.2, 0.36],
        [19.9, 29.6, 0.40],
        [20.2, 30.0, 0.40],
    ]),
}