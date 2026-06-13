import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, linprog


try:
    datos_modelo = np.load('imm904_procesado.npz', allow_pickle=True)
except FileNotFoundError:
    print("Error: No se encontró 'imm904_procesado.npz'.")
    print("Por favor, ejecuta primero 'python procesar_modelo.py' para generar las matrices.")
    exit()

# Estas son las variables globales que usaremos en las funciones matemáticas
S = datos_modelo['S']
RXNS = list(datos_modelo['rxns'])  # Convertimos a lista para usar .index() fácilmente
METS = list(datos_modelo['mets'])
LB_ORIGINAL = datos_modelo['lb']
UB_ORIGINAL = datos_modelo['ub']
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
# ###########################################################################
# IDENTIFICADORES DE REACCIONES — esto no lo editamos, es para mapear las reacciones del modelo
# ############################################################################

#   _c = citosol   _m = mitocondria   _x = peroxisoma   _e = extracelular
REA = {
    "biomasa":  "biomass_SC5_notrace",   # reacción de crecimiento
    "glucosa":  "EX_glc__D_e",          # consumo de glucosa
    "etanol":   "EX_etoh_e",            # producción de etanol
    "oxigeno":  "EX_o2_e",              # intercambio O2
    "co2":      "EX_co2_e",             # intercambio CO2
    "mantenimiento": "ATPM",            # mantenimiento de ATP
    # Compartimentos para ergosterol/ácidos grasos (anaeróbico)
    "ergosterol":   "EX_ergst_e",
    "zymosterol":   "EX_zymst_e",
    "hdcea":        "EX_hdcea_e",
    "ocdca":        "EX_ocdca_e",
    "ocdcea":       "EX_ocdcea_e",
    "ocdcya":       "EX_ocdcya_e",
}

# Función objetivo por compartimento  (sign: +1 max, -1 min)
OBJ_COMPARTIMENTOS = {
    # P1 – citosol: minimizar producción de NADH  (mejor para anaeróbico)
    "citosol": [
        ("NADH2_u10m",   -1),   # NADH deshidrogenasa (aprox citosol)
        ("GAPD",         -1),   # gliceraldehído-3-P deshidrogenasa → NADH
        ("PDHm",         -1),   # piruvato deshidrogenasa → NADH
        ("CS",           -1),   # citrato sintasa
        ("TPI",          +1),   # triosa fosfato isomerasa
        ("PGK",          +1),   # fosfoglicerato quinasa → ATP
    ],
    # P2 – mitocondria: minimizar consumo de NADH/NADPH  (mejor para aeróbico)
    "mitocondria": [
        ("NADH2_u10m",   -1),
        ("SUCD1m",       -1),   # succinato deshidrogenasa
        ("ICDHyrm",      -1),   # isocitrato deshidrogenasa
        ("MDHm",         -1),   # malato deshidrogenasa
        ("AKGDam",       -1),
    ],
    # P3 – peroxisoma: maximizar producción de ácidos grasos
    "peroxisoma": [
        ("r_0658",       +1),   # beta-oxidación (iMM904 IDs varían)
        ("FACOAL160x",   +1),
        ("ACITL",        +1),
        ("FAO180ACPx",   +1),
    ],
    # P4 – virtual: maximizar biomasa (objetivo clásico)
    "virtual": [
        ("biomass_SC5_notrace", +1),
    ],
}


# ####################################################################################
# FUNCIONES AUXILIARES
# ####################################################################################
def aplicar_condicion_matemática(condicion: str, glucosa_uptake: float):
    """
    Ajusta los vectores lb y ub según la Tabla 5 del paper usando índices de NumPy.
    Retorna copias modificadas de (lb, ub).
    """
    # 1. Clonar los límites originales para no sobreescribir el modelo base
    lb_mod = np.copy(LB_ORIGINAL)
    ub_mod = np.copy(UB_ORIGINAL)
    
    # Helper rápido para obtener el índice de una reacción por su ID en la lista RXNS
    def get_idx(rxn_id):
        return RXNS.index(REA[rxn_id]) if REA[rxn_id] in RXNS else None

    # 2. Fijar consumo de Glucosa (Fila 1 de la Tabla 5)
    idx_glc = get_idx("glucosa")
    if idx_glc is not None:
        # En FBA, el consumo/uptake se modela como flujo negativo
        lb_mod[idx_glc] = -glucosa_uptake
        ub_mod[idx_glc] = -glucosa_uptake  # Forzamos a que sea exactamente el valor experimental

    # 3. Configurar según condiciones Anaeróbica / Aeróbica
    if condicion == "anaerobic":
        # Oxígeno = 0 (Fila 2)
        idx_o2 = get_idx("oxigeno")
        if idx_o2 is not None:
            lb_mod[idx_o2] = 0.0
            ub_mod[idx_o2] = 0.0
            
        # CO2 libre hacia afuera (Fila 3: lb=0, ub=1000)
        idx_co2 = get_idx("co2")
        if idx_co2 is not None:
            lb_mod[idx_co2] = 0.0
            ub_mod[idx_co2] = 1000.0
            
        # Nutrientes anaeróbicos: Esteroles y Ácidos Grasos (Fila 4: lb=-1000, ub=1000)
        # Permite que la célula los absorba del medio simulado ya que no los puede fabricar sin O2
        nutrientes_anaerobicos = ["ergosterol", "zymosterol", "hdcea", "ocdca", "ocdcea", "ocdcya"]
        for nut in nutrientes_anaerobicos:
            idx_nut = get_idx(nut)
            if idx_nut is not None:
                lb_mod[idx_nut] = -1000.0
                ub_mod[idx_nut] = 1000.0

    else:  # aerobic
        # Oxígeno libre ilimitado (Fila 5: lb=-1000, ub=0)
        idx_o2 = get_idx("oxigeno")
        if idx_o2 is not None:
            lb_mod[idx_o2] = -1000.0
            ub_mod[idx_o2] = 0.0  # Solo consumo, la levadura no "exhala" O2 puro

    return lb_mod, ub_mod

# def aplicar_condicion(model, condicion: str, glucosa_uptake: float):
#     """
#     Configura los límites del modelo para la condición experimental.
#     Sigue la Tabla 5 del paper.
#     """
#     # Fijar consumo de glucosa al valor experimental
#     try:
#         model.reactions.get_by_id(REA["glucosa"]).lower_bound = -glucosa_uptake
#     except KeyError:
#         pass
 
#     if condicion == "anaerobic":
#         # Sin oxígeno
#         try:
#             r = model.reactions.get_by_id(REA["oxigeno"])
#             r.lower_bound = 0
#             r.upper_bound = 0
#         except KeyError:
#             pass
#         # CO2 fijo (no libre)
#         try:
#             r = model.reactions.get_by_id(REA["co2"])
#             r.lower_bound = 0
#             r.upper_bound = 1000
#         except KeyError:
#             pass
#         # Permitir uptake de esteroles/ácidos grasos que no se sintetizan sin O2
#         for r_id in ["ergosterol", "zymosterol", "hdcea", "ocdca", "ocdcea", "ocdcya"]:
#             try:
#                 model.reactions.get_by_id(REA[r_id]).lower_bound = -1000
#             except KeyError:
#                 pass
 
#     else:  # aerobic
#         # O2 ilimitado (paper cambia límite de -2 a -1000)
#         try:
#             r = model.reactions.get_by_id(REA["oxigeno"])
#             r.lower_bound = -1000
#         except KeyError:
#             pass
def construir_objetivo_compartimento(model, comp_reacciones: list) -> dict:
    """
    Dado una lista de (rxn_id, signo), devuelve un dict {rxn_id: coeficiente}
    con solo las reacciones que existen en el modelo.
    """
    obj = {}
    rxn_ids = {r.id for r in model.reactions}
    for rxn_id, signo in comp_reacciones:
        if rxn_id in rxn_ids:
            obj[rxn_id] = float(signo)
    return obj

def fba_combinado(pesos: dict, condicion: str, glucosa: float):
    """
    Resuelve el FBA con la función objetivo combinada usando scipy.optimize.linprog.
    """
    # 1. Normalizar pesos
    total = sum(abs(v) for v in pesos.values())
    if total == 0:
        return None, None
    pesos_norm = {k: v / total for k, v in pesos.items()}

    # 2. Construir vector de la función objetivo c
    c_vector = np.zeros(len(RXNS))
    for comp, peso in pesos_norm.items():
        if comp not in OBJ_COMPARTIMENTOS or abs(peso) < 1e-10:
            continue
        for rxn_id, signo in OBJ_COMPARTIMENTOS[comp]:
            if rxn_id in RXNS:
                idx = RXNS.index(rxn_id)
                c_vector[idx] += peso * float(signo)

    # 3. Obtener límites modificados para la simulación
    lb_mod, ub_mod = aplicar_condicion_matemática(condicion, glucosa)
    bounds = list(zip(lb_mod, ub_mod))

    # 4. Resolver el problema lineal: S * v = 0
    # b_eq es un vector de ceros del tamaño de las filas de S (metabolitos)
    b_eq = np.zeros(S.shape[0]) 
    
    res = linprog(c=-c_vector, A_eq=S, b_eq=b_eq, bounds=bounds, method='highs')

    if not res.success:
        return 0.0, 0.0

    # 5. Extraer resultados mediante índices
    idx_biomasa = RXNS.index(REA["biomasa"])
    idx_etanol = RXNS.index(REA["etanol"])
    
    crec = res.x[idx_biomasa]
    etoh = abs(res.x[idx_etanol]) # Producción neta (si sale de la célula es positivo)

    return float(crec), float(etoh)


def fba_clasico(condicion: str, glucosa: float):
    """
    FBA tradicional: Maximiza únicamente la reacción de biomasa.
    """
    c_vector = np.zeros(len(RXNS))
    idx_biomasa = RXNS.index(REA["biomasa"])
    c_vector[idx_biomasa] = 1.0  # Coeficiente 1 solo a la biomasa

    lb_mod, ub_mod = aplicar_condicion_matemática(condicion, glucosa)
    bounds = list(zip(lb_mod, ub_mod))
    b_eq = np.zeros(S.shape[0])

    res = linprog(c=-c_vector, A_eq=S, b_eq=b_eq, bounds=bounds, method='highs')

    if not res.success:
        return 0.0, 0.0

    idx_etanol = RXNS.index(REA["etanol"])
    crec = res.x[idx_biomasa]
    etoh = abs(res.x[idx_etanol])

    return float(crec), float(etoh)

# def fba_combinado(model, pesos: dict, condicion: str, glucosa: float):
#     """
#     Parámetros
    
#         model     : cobra.model
#         pesos     : dict  {compartimento: peso}  se normalizan a suma=1
#         condicion : "anaerobic" o "aerobic"
#         glucosa   : uptake de glucosa en mmol/gPS·h
 
#     Retorna: (crecimiento, etanol)
#     """
 
#     # Normalizar pesos
#     total = sum(abs(v) for v in pesos.values())
#     if total == 0:
#         return None, None
#     pesos_norm = {k: v / total for k, v in pesos.items()}
 
#     with model:
#         aplicar_condicion(model, condicion, glucosa)
 
#         # Construir función objetivo combinada: suma ponderada de compartimentos
#         objetivo_combinado = {}
#         for comp, peso in pesos_norm.items():
#             if comp not in OBJ_COMPARTIMENTOS or abs(peso) < 1e-10:
#                 continue
#             obj_comp = construir_objetivo_compartimento(
#                 model, OBJ_COMPARTIMENTOS[comp]
#             )
#             for rxn_id, coef in obj_comp.items():
#                 objetivo_combinado[rxn_id] = (
#                     objetivo_combinado.get(rxn_id, 0.0) + peso * coef
#                 )
 
#         if not objetivo_combinado:
#             return None, None
 
#         # Asignar objetivo al modelo
#         model.objective = objetivo_combinado
 
#         try:
#             sol = model.optimize()
#             if sol.status != "optimal":
#                 return None, None
#         except Exception:
#             return None, None
 
#         # Extraer crecimiento y etanol
#         try:
#             crec = sol.fluxes[REA["biomasa"]]
#         except KeyError:
#             crec = 0.0
#         try:
#             etoh = abs(sol.fluxes.get(REA["etanol"], 0.0))
#         except Exception:
#             etoh = 0.0
 
#         return float(crec), float(etoh)

def error_promedio(pred, exp):
    """Error relativo promedio |pred - exp| / |exp|. Ignora NaN."""
    pred = np.array(pred, dtype=float)
    exp  = np.array(exp,  dtype=float)
    mask = (np.isfinite(pred)) & (np.isfinite(exp)) & (exp != 0)
    if mask.sum() == 0:
        return np.nan
    return float(np.mean(np.abs(pred[mask] - exp[mask]) / np.abs(exp[mask])))

def optimizar_pesos(condicion: str, datos_exp: np.ndarray):
    '''Busca los pesos w_k que minimizan el error entre predicciones y datos experimentales
    usando scipy.optimize.minimize. Solo optimiza para la condición seleccionada en CONFIG.'''

    col_glc  = 0
    col_crec = 2 if condicion == "anaerobic" else 3

    #contruimos funcion objetivo
    def funcion_objetivo(w):
        # Normalizar pesos
        w = np.array(w)
        w_sum=w.sum()
        w_norm = w / w.sum() if w.sum() > 0 else np.ones_like(w) / len(w)
        pesos_dict = {
            "citosol":     w_norm[0],
            "mitocondria": w_norm[1],
            "peroxisoma":  w_norm[2],
            "virtual":     w_norm[3],
        }
        predicciones = []
        for fila in datos_exp:
            glucosa = fila[col_glc]
            pred_crec, _ = fba_combinado(pesos_dict, condicion, glucosa)
            predicciones.append(pred_crec)
        return error_promedio(predicciones, datos_exp[:, col_crec])
    #Punto de partida(pesos iguales)
    w_i = np.array([0.25, 0.25, 0.25, 0.25])
    result = minimize (funcion_objetivo, w_i, method="Nelder-Mead")

    #Aseguramos valores absolutos no negativos y normalizamos 
    w_opt = np.abs(result.x)
    w_opt_n = w_opt / w_opt.sum()
    pesos_optimos = {
        "citosol":     float(w_opt_n[0]),
        "mitocondria": float(w_opt_n[1]),
        "peroxisoma":  float(w_opt_n[2]),
        "virtual":     float(w_opt_n[3]),
    }
    return pesos_optimos

# def fba_clasico(model, condicion: str, glucosa: float):
#     """FBA tradicional: maximizar biomasa."""

#     with model:
#         aplicar_condicion(model, condicion, glucosa)
#         try:
#             bio_r = model.reactions.get_by_id(REA["biomasa"])
#         except KeyError:
#             return None, None
        
#         model.objective = bio_r
#         try:
#             sol = model.optimize()
#             if sol.status != "optimal":
#                 return None, None
#         except Exception:
#             return None, None
 
#         crec = float(sol.fluxes.get(REA["biomasa"], 0.0))
#         etoh = float(abs(sol.fluxes.get(REA["etanol"], 0.0)))
#         return crec, etoh
    
# ##############################################################################################
# SIMULACION
# ##############################################################################################
def simular(model, condicion: str, datos: np.ndarray, config: dict):
    """
    Corre todos los modos de simulación sobre todos los puntos experimentales.
    """
    col_glc  = 0
    col_etoh = 1
    col_crec = 2 if condicion == "anaerobic" else 3
 
    glc_vals   = datos[:, col_glc]
    crec_exp   = datos[:, col_crec]
    etoh_exp   = datos[:, col_etoh]
 
    resultados = {
        "glucosa_exp": glc_vals,
        "crecimiento_exp": crec_exp,
        "etanol_exp":      etoh_exp,
    }
 
    mostrar = config["mostrar"]
    pesos   = config["pesos"]
 
    # Clásico biomasa 
    if mostrar.get("clasico_biomasa"):
        crec_b, etoh_b = [], []
        # print("  Simulando FBA clásico (max biomasa)...")
        for glc in glc_vals:
            c, e = fba_clasico(model, condicion, glc)
            crec_b.append(c)
            etoh_b.append(e)
        resultados["crecimiento_clasico"] = crec_b
        resultados["etanol_clasico"]      = etoh_b
        eg  = error_promedio(crec_b, crec_exp)
        ee  = error_promedio(etoh_b, etoh_exp)
        # print(f"    e_gw = {eg:.3f}  |  e_etoh = {ee:.3f}")
 
    # w = [1,1,1,1]
    if mostrar.get("paper_original"):
        pesos_paper = {"citosol": 1.0, "mitocondria": 1.0,
                       "peroxisoma": 1.0, "virtual": 1.0}
        crec_p, etoh_p = [], []
        # print("  Simulando función objetivo del paper (w=[1,1,1,1])...")
        for glc in glc_vals:
            c, e = fba_combinado(model, pesos_paper, condicion, glc)
            crec_p.append(c)
            etoh_p.append(e)
        resultados["crecimiento_paper"] = crec_p
        resultados["etanol_paper"]      = etoh_p
        eg  = error_promedio(crec_p, crec_exp)
        ee  = error_promedio(etoh_p, etoh_exp)
        # print(f"    e_gw = {eg:.3f}  |  e_etoh = {ee:.3f}")
 
    # Adaptativo (pesos personalizados)
    if mostrar.get("adaptativo"):
        if config.get("optimizar_pesos"):
            pesos = optimizar_pesos(condicion, datos)
        crec_a, etoh_a = [], []
        # print(f"  Simulando adaptativo con pesos: {pesos}...")
        for glc in glc_vals:
            c, e = fba_combinado(model, pesos, condicion, glc)
            crec_a.append(c)
            etoh_a.append(e)
        resultados["crecimiento_adaptativo"] = crec_a
        resultados["etanol_adaptativo"]      = etoh_a
        eg  = error_promedio(crec_a, crec_exp)
        ee  = error_promedio(etoh_a, etoh_exp)
        # print(f"    e_gw = {eg:.3f}  |  e_etoh = {ee:.3f}")
 
    return resultados
 