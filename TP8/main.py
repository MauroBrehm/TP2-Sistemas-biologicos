import numpy as np
import matplotlib.pyplot as plt

# ─── Parámetros del modelo HH ──────────────────────────────────
# Conductancias máximas de cada tipo de canal (cuando todos están abiertos)
gNaBar = 120.0  # mS/cm² — sodio
gKBar  = 36.0   # mS/cm² — potasio
gLBar  = 0.3    # mS/cm² — fuga (leak): siempre abierto, no tiene compuerta

# Potenciales de equilibrio de Nernst para cada ion
# Son los voltajes a los que la corriente neta de ese ion es cero
ENa = 45.0   # mV
EK  = -82.0  # mV
EL  = -59.0  # mV

# Capacitancia de membrana: cuánta carga necesitás para cambiar el voltaje
Cm = 1.0  # µF/cm²

# ─── Funciones de tasas de transición (α y β) ──────────────────
# α es la tasa de apertura (cerrado → abierto)
# β es la tasa de cierre  (abierto → cerrado)
# Ambas dependen del voltaje V en cada momento
# Las fórmulas vienen de ajustes experimentales de Hodgkin y Huxley

def alpha_m(V):
    'Tasa de apertura de la compuerta m (activación del canal Na)'
    'Tiene una indeterminación 0/0 cuando V = -45 mV → se usa LHopital → resultado = 1.0'
    return np.where(np.abs(V + 45) < 1e-7, 1.0,
                    (V + 45) / 10 / (1 - np.exp(-(V + 45) / 10)))

def beta_m(V):
    'Tasa de cierre de la compuerta m'
    return 4.0 * np.exp(-(V + 70) / 18)

def alpha_h(V):
    'Tasa de apertura de la compuerta h (inactivación del canal Na)'
    'h se llama "compuerta de inactivación": cuando está cerrada, el canal Na no conduce'
    return 0.07 * np.exp(-(V + 70) / 20)

def beta_h(V):
    'Tasa de cierre de la compuerta h'
    return 1.0 / (1 + np.exp(-(V + 40) / 10))

def alpha_n(V):
    'Tasa de apertura de la compuerta n (activación del canal K)'
    'También tiene indeterminación en V = -60 mV → LHopital → resultado = 0.1'
    return np.where(np.abs(V + 60) < 1e-7, 0.1,
                    (V + 60) / 100 / (1 - np.exp(-(V + 60) / 10)))

def beta_n(V):
    'Tasa de cierre de la compuerta n'
    return 0.125 * np.exp(-(V + 70) / 80)

# ─── Condiciones iniciales en estado estable ───────────────────
# En reposo (V = -70 mV), las compuertas están en equilibrio:
# la probabilidad de estar abierta = α / (α + β)
# Esto sale de igualar dm/dt = 0 en la ecuación diferencial
V0 = -70.0
m0 = alpha_m(V0) / (alpha_m(V0) + beta_m(V0))
h0 = alpha_h(V0) / (alpha_h(V0) + beta_h(V0))
n0 = alpha_n(V0) / (alpha_n(V0) + beta_n(V0))

# ─── Configuración de la simulación ────────────────────────────
dt = 0.01   # paso de tiempo en ms (debe ser pequeño para precisión)
T  = 50.0   # duración total de la simulación en ms
t  = np.arange(0, T + dt, dt)  # vector de tiempo

# Reservar memoria para todas las variables (más rápido que append)
V = np.zeros(len(t))
m = np.zeros(len(t))
h = np.zeros(len(t))
n = np.zeros(len(t))

# Asignar condiciones iniciales
V[0], m[0], h[0], n[0] = V0, m0, h0, n0

# ─── Corriente de entrada ───────────────────────────────────────
# Pulso rectangular: 15 µA/cm² entre t=5ms y t=6ms
# Simula un estímulo externo que dispara el potencial de acción
I = np.zeros(len(t))
I[(t >= 5) & (t < 6)] = 15.0

# ─── Integración numérica por método de Euler ──────────────────
for k in range(len(t) - 1):
    Vk, mk, hk, nk = V[k], m[k], h[k], n[k]

    # Conductancias instantáneas según el estado de las compuertas
    # gNa = gNaBar * m³ * h  (3 subunidades m + 1 subunidad h, todas deben estar abiertas)
    # gK  = gKBar  * n⁴      (4 subunidades n, todas deben estar abiertas)
    gNa = gNaBar * mk**3 * hk
    gK  = gKBar  * nk**4

    # Ecuación de membrana (ley de Kirchhoff):
    # I_ext = Cm * dV/dt + I_Na + I_K + I_L
    # Despejando dV/dt:
    dVdt = (I[k] - gNa * (Vk - ENa) - gK * (Vk - EK) - gLBar * (Vk - EL)) / Cm

    # Actualizar voltaje
    V[k+1] = Vk + dVdt * dt

    # Actualizar compuertas (ecuaciones diferenciales de primer orden)
    # dm/dt = αm*(1-m) - βm*m
    #   → el primer término: canales cerrados que se abren
    #   → el segundo término: canales abiertos que se cierran
    m[k+1] = mk + (alpha_m(Vk) * (1 - mk) - beta_m(Vk) * mk) * dt
    h[k+1] = hk + (alpha_h(Vk) * (1 - hk) - beta_h(Vk) * hk) * dt
    n[k+1] = nk + (alpha_n(Vk) * (1 - nk) - beta_n(Vk) * nk) * dt

# ─── Calcular conductancias finales para graficar ──────────────
gNA_vec = gNaBar * m**3 * h
gK_vec  = gKBar  * n**4

# ─── Graficación ───────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
fig.suptitle('Modelo determinístico — Codigo del profe')
axes[0].plot(t, V);       axes[0].set_ylabel('V (mV)')
axes[1].plot(t, gNA_vec); axes[1].set_ylabel('gNa (mS/cm²)')
axes[2].plot(t, gK_vec);  axes[2].set_ylabel('gK (mS/cm²)')
axes[2].set_xlabel('Tiempo (ms)')
plt.tight_layout()
plt.show()

############################################################
#HASTA EL CODIGO DEL PROFE
###################################################

####################################################
#ACT 1 - Implementar modelo estocástico con N=100
######################################################
#se toma el codigo del profe como base y se remplaza gk determinidtico
#por una conductancia estocastica calculada a partir de N canales indicviduales 

N = 100

No = int(np.random.binomial(N, n0))
Nc = N - No

V_arr  = np.zeros(len(t));  V_arr[0]  = V0
m_arr  = np.zeros(len(t));  m_arr[0]  = m0
h_arr  = np.zeros(len(t));  h_arr[0]  = h0
gK_arr = np.zeros(len(t));  gK_arr[0] = gKBar * No / N

for k in range(len(t) - 1):
    Vk, mk, hk = V_arr[k], m_arr[k], h_arr[k]

    p_open  = alpha_n(Vk) * dt
    p_close = beta_n(Vk)  * dt
    n_open  = np.random.binomial(Nc, min(p_open,  1.0))
    n_close = np.random.binomial(No, min(p_close, 1.0))
    No = No + n_open - n_close
    Nc = N - No

    gKest = gKBar * (No / N)**4
    gNa   = gNaBar * mk**3 * hk
    dVdt  = (I[k] - gNa*(Vk-ENa) - gKest*(Vk-EK) - gLBar*(Vk-EL)) / Cm

    V_arr[k+1]  = Vk + dVdt * dt
    m_arr[k+1]  = mk + (alpha_m(Vk)*(1-mk) - beta_m(Vk)*mk) * dt
    h_arr[k+1]  = hk + (alpha_h(Vk)*(1-hk) - beta_h(Vk)*hk) * dt
    gK_arr[k+1] = gKest

########################################################
#ACT 2 - graficar V, gNa y gKest con pulsos de 15ua a t=5ms 
############################################################
#El pulso ya esta definido en el vector I del profe

gNa_arr = gNaBar * m_arr**3 * h_arr

fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
fig.suptitle('Actividad 2 — Modelo estocástico (N=100), pulso 15 µA/cm²')

axes[0].plot(t, V_arr);   axes[0].set_ylabel('V (mV)')
axes[1].plot(t, gNa_arr); axes[1].set_ylabel('gNa (mS/cm²)')
axes[2].plot(t, gK_arr);  axes[2].set_ylabel('gKest (mS/cm²)')
axes[2].set_xlabel('Tiempo (ms)')
plt.tight_layout()
plt.show()


###########################################################
#ACT 3 - Comparar N  entre si y con modelo determinístico
###########################################################
fig,axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
fig.suptitle('Modelo estocastico (fuerza bruta) — Comparacion de N')

colores={10:'pink', 100:'cyan', 1000:'purple'}

for N in [10, 100, 1000]:
        # ─── Parámetro nuevo: número de canales de K ───────────────────
        # Cuanto mayor N, más se parece al modelo determinístico
        # Con N pequeño (ej. 10), el ruido estocástico es muy visible
    
        # Canales de K: al inicio, cada canal está abierto con probabilidad n0
        # np.random.binomial simula cuántos de N ensayos Bernoulli salen "éxito"
        No = int(np.random.binomial(N, n0))  # número inicial de canales abiertos --> la probabilidad de que un canal de K esté abierto es n0^4 (todas las subunidades n deben estar abiertas)
        Nc = N - No                           # número inicial de canales cerrados
    
        # ─── Arrays para guardar resultados ────────────────────────────
        V_arr  = np.zeros(len(t))
        m_arr  = np.zeros(len(t))
        h_arr  = np.zeros(len(t))
        gK_arr = np.zeros(len(t))  # conductancia estocástica de K
    
        V_arr[0] = V0
        m_arr[0] = m0
        h_arr[0] = h0
        gK_arr[0] = gKBar * No / N  # fracción de canales abiertos × conductancia máxima


        # ─── Bucle de integración estocástica ──────────────────────────
        for k in range(len(t) - 1):
            Vk = V_arr[k]
            mk = m_arr[k]
            hk = h_arr[k]

            # --- Paso estocástico: actualizar canales de K ---

            # Probabilidad de que una compuerta CERRADA se ABRA en este dt
            p_open  = alpha_n(Vk) * dt  # debe ser << 1 para que el método sea válido

            # Probabilidad de que una compuerta ABIERTA se CIERRE en este dt
            p_close = beta_n(Vk)  * dt

            # Entre los Nc canales cerrados, cuantos se abren?
            # Cada uno tiene probabilidad p_open --> distribución binomial
            n_open  = np.random.binomial(Nc, min(p_open,  1.0))


            # Entre los No canales abiertos, cuantos se cierran?
            n_close = np.random.binomial(No, min(p_close, 1.0))

            # Actualizar conteos
            No = No + n_open  - n_close
            Nc = N  - No  # siempre se cumple No + Nc = N

            # Conductancia estocástica de K: proporción de canales abiertos
            gKest = gKBar * (No / N)**4

            # --- Paso determinístico: canal de Na sigue siendo continuo ---
            gNa = gNaBar * mk**3 * hk

            # Ecuación de membrana con gK estocástico
            dVdt = (I[k] - gNa*(Vk-ENa) - gKest*(Vk-EK) - gLBar*(Vk-EL)) / Cm

            V_arr[k+1]  = Vk + dVdt * dt
            m_arr[k+1]  = mk + (alpha_m(Vk)*(1-mk) - beta_m(Vk)*mk) * dt
            h_arr[k+1]  = hk + (alpha_h(Vk)*(1-hk) - beta_h(Vk)*hk) * dt
            gK_arr[k+1] = gKest

        # ─── Graficación ───────────────────────────────────────────────
        gNa_arr = gNaBar * m_arr**3 * h_arr

        axes[0].plot(t, V_arr, label=f'N={N}', color=colores[N])
        axes[1].plot(t, gNa_arr, label=f'N={N}', color=colores[N])
        axes[2].plot(t, gK_arr, label=f'N={N}', color=colores[N])

#Curva del modelo determinístico para comparar
axes[0].plot(t, V, label='Determinístico', color='black', linestyle='--')
axes[1].plot(t, gNA_vec, color='black', linestyle='--')
axes[2].plot(t, gK_vec, color='black', linestyle='--')

axes[0].set_ylabel('V (mV)')
axes[1].set_ylabel('gNa (mS/cm²)')
axes[2].set_ylabel('gK (mS/cm²)')
axes[2].set_xlabel('Tiempo (ms)')
plt.tight_layout()
plt.show()

###########################################################
#ACT 4 - V, gNa y gKest aplicando pulsos de corriente con distintas amplitudes y duraciones
###########################################################

def generar_pulso_estocastico(t_inicio, dur, amplitud, N=100, T=50.0, dt=0.01):
    t = np.arange(0, T + dt, dt)

    # Condiciones iniciales en estado estable (V = -70 mV)
    V0 = -70.0
    m0 = alpha_m(V0) / (alpha_m(V0) + beta_m(V0))
    h0 = alpha_h(V0) / (alpha_h(V0) + beta_h(V0))
    n0 = alpha_n(V0) / (alpha_n(V0) + beta_n(V0))

    V = np.zeros(len(t))
    m = np.zeros(len(t))
    h = np.zeros(len(t))
    n = np.zeros(len(t))
    gKest = np.zeros(len(t))
    V[0], m[0], h[0], n[0] = V0, m0, h0, n0

    I = np.zeros(len(t))
    I[(t >= t_inicio) & (t < t_inicio + dur)] = amplitud

    No = int(np.random.binomial(N, n0))
    Nc = N - No
    gKest[0] = gKBar * (No / N)**4

    for k in range(len(t) - 1):
        Vk, mk, hk, nk = V[k], m[k], h[k], n[k]

        p_open = min(alpha_n(Vk) * dt, 1.0)
        p_close = min(beta_n(Vk) * dt, 1.0)

        n_open = np.random.binomial(Nc, p_open)
        n_close = np.random.binomial(No, p_close)

        No = No + n_open - n_close
        Nc = N - No
        gKest[k + 1] = gKBar * (No / N)**4

        gNa = gNaBar * mk**3 * hk
        dVdt = (I[k] - gNa * (Vk - ENa) - gKest[k] * (Vk - EK) - gLBar * (Vk - EL)) / Cm

        V[k+1] = Vk + dVdt * dt
        m[k+1] = mk + (alpha_m(Vk) * (1 - mk) - beta_m(Vk) * mk) * dt
        h[k+1] = hk + (alpha_h(Vk) * (1 - hk) - beta_h(Vk) * hk) * dt
        n[k+1] = nk + (alpha_n(Vk) * (1 - nk) - beta_n(Vk) * nk) * dt

    gNa = gNaBar * m**3 * h
    return t, V, gNa, gKest, I

durac      = [0.5, 1, 3]
amplitudes = [7.5, 15, 30]
colores    = ['deeppink', 'yellowgreen', 'steelblue']  # uno por amplitud

fig_v,   axes_v   = plt.subplots(1, len(durac), figsize=(14, 4), sharex=True, sharey=True)
fig_gNa, axes_gNa = plt.subplots(1, len(durac), figsize=(14, 4), sharex=True, sharey=True)
fig_gK,  axes_gK  = plt.subplots(1, len(durac), figsize=(14, 4), sharex=True, sharey=True)

fig_v.suptitle('V — misma duración, distintas amplitudes')
fig_gNa.suptitle('gNa — misma duración, distintas amplitudes')
fig_gK.suptitle('gKest — misma duración, distintas amplitudes')

for j, dur in enumerate(durac):          # columna = duración
    for amp, color in zip(amplitudes, colores):   # línea = amplitud

        t, V, gNa, gKest, I = generar_pulso_estocastico(t_inicio=5, dur=dur, amplitud=amp)

        axes_v[j].plot(t, V, color=color, lw=1.2, label=f'{amp} µA/cm²')
        axes_gNa[j].plot(t, gNa, color=color, lw=1.2, label=f'{amp} µA/cm²')
        axes_gK[j].plot(t, gKest, color=color, lw=1.2, label=f'{amp} µA/cm²')

    for ax in [axes_v[j], axes_gNa[j], axes_gK[j]]:
        ax.set_title(f'Duración = {dur} ms', fontsize=9)
        ax.set_xlabel('t (ms)')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

axes_v[0].set_ylabel('V (mV)')
axes_gNa[0].set_ylabel('gNa (mS/cm²)')
axes_gK[0].set_ylabel('gKest (mS/cm²)')

for fig in [fig_v, fig_gNa, fig_gK]:
    fig.tight_layout()
plt.show()

