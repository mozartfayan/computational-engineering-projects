import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. ALGORITMO TDMA
# ============================================================

def solve_tdma(a_W, a_P, a_E, b_rhs):
  
    n = len(a_P)
    P = np.zeros(n)
    Q = np.zeros(n)
    
    # 1. Varredura Direta (Forward Sweep)
    P[0] = a_E[0] / a_P[0]
    Q[0] = b_rhs[0] / a_P[0]
    
    for i in range(1, n):
        denom = a_P[i] - a_W[i] * P[i-1]
        if i < n - 1:
            P[i] = a_E[i] / denom
        Q[i] = (b_rhs[i] + a_W[i] * Q[i-1]) / denom
        
    # 2. Substituição Regressiva (Back Substitution)
    T = np.zeros(n)
    T[-1] = Q[-1]
    for i in range(n - 2, -1, -1):
        T[i] = P[i] * T[i+1] + Q[i]
        
    return T, P, Q

# ============================================================
# 2. SOLUÇÃO ANALÍTICA
# ============================================================

def T_analytical(x, u, L=0.1, k=0.595, rho=997.0, T0=150.0, TL=50.0):
    """
    Solução analítica exata da equação de convecção-difusão 1D em regime permanente:
    (T(x) - T0) / (TL - T0) = (exp(rho*u*x / k) - 1) / (exp(rho*u*L / k) - 1)
    """
    exp_term = np.exp((rho * u / k) * x)
    exp_L = np.exp((rho * u / k) * L)
    return T0 + (TL - T0) * (exp_term - 1.0) / (exp_L - 1.0)


# ============================================================
# 3. SIMULAÇÃO MVF EDC
# ============================================================

def solve_edc(u, N, L=0.1, k=0.595, rho=997.0, T0=150.0, TL=50.0):
    """
    Solução por Volumes Finitos via Esquema de Diferenças Centrais (EDC).
    """
    dx = L / N
    x_nodes = np.array([(i + 0.5) * dx for i in range(N)])
    F = rho * u
    D = k / dx
    
    a_W = np.zeros(N)
    a_E = np.zeros(N)
    a_P = np.zeros(N)
    b_rhs = np.zeros(N)
    
    for i in range(N):
        if i == 0:
            a_W[i] = 0.0
            a_E[i] = D - F / 2.0
            a_P[i] = 3.0 * D + F / 2.0
            b_rhs[i] = (2.0 * D + F) * T0
        elif i == N - 1:
            a_W[i] = D + F / 2.0
            a_E[i] = 0.0
            a_P[i] = 3.0 * D - F / 2.0
            b_rhs[i] = (2.0 * D - F) * TL
        else:
            a_W[i] = D + F / 2.0
            a_E[i] = D - F / 2.0
            a_P[i] = 2.0 * D
            b_rhs[i] = 0.0

    T_num, P, Q = solve_tdma(a_W, a_P, a_E, b_rhs)
    return x_nodes, T_num, a_W, a_E, a_P, b_rhs, F, D, P, Q


# ============================================================
# 3. SIMULAÇÃO MVF UPWIND
# ============================================================

def solve_upwind(u, N, L=0.1, k=0.595, rho=997.0, T0=150.0, TL=50.0):
    """
    Solução por Volumes Finitos via Esquema Upwind.
    """
    dx = L / N
    x_nodes = np.array([(i + 0.5) * dx for i in range(N)])
    F = rho * u
    D = k / dx
    
    a_W = np.zeros(N)
    a_E = np.zeros(N)
    a_P = np.zeros(N)
    b_rhs = np.zeros(N)
    
    for i in range(N):
        if i == 0:
            a_W[i] = 0.0
            a_E[i] = D
            a_P[i] = 3.0 * D + F
            b_rhs[i] = (2.0 * D + F) * T0
        elif i == N - 1:
            a_W[i] = D + F
            a_E[i] = 0.0
            a_P[i] = 3.0 * D + F
            b_rhs[i] = 2.0 * D * TL
        else:
            a_W[i] = D + F
            a_E[i] = D
            a_P[i] = 2.0 * D + F
            b_rhs[i] = 0.0

    T_num, P, Q = solve_tdma(a_W, a_P, a_E, b_rhs)
    return x_nodes, T_num, a_W, a_E, a_P, b_rhs, F, D, P, Q


# ============================================================
# 4. EXECUÇÃO PRINCIPAL
# ============================================================

if __name__ == "__main__":
    cases = [
        ("Item a", 0.01, 5),
        ("Item b", 0.10, 5),
        ("Item c", 0.10, 20)
    ]

    print("=" * 80)
    print("RESOLUÇÃO DA LISTA DE EXERCÍCIOS 7 - CONVECÇÃO-DIFUSÃO 1D (MVF + TDMA)")
    print("=" * 80)
    print("Propriedades: k = 0.595 W/mK | rho = 997 kg/m³ | T(0) = 150 °C | T(L) = 50 °C | L = 0.1 m\n")

    for item, u, N in cases:
        x_nodes, T_edc, aW_e, aE_e, aP_e, b_e, F, D, P_e, Q_e = solve_edc(u, N)
        _, T_upw, aW_u, aE_u, aP_u, b_u, _, _, P_u, Q_u = solve_upwind(u, N)
        T_anal = T_analytical(x_nodes, u=u)
        Pe = F / D

        print("-" * 80)
        print(f"--- {item}: u = {u:.2f} m/s | N = {N} volumes | dx = {0.1/N:.4f} m ---")
        print(f"Fluxo Convectivo F = {F:.4f} kg/(m²s) | Condutância Difusiva D = {D:.4f} W/(m²K)")
        print(f"Número de Peclet de Malha Pe = F/D = {Pe:.4f} -> {'STÁVEL (Pe <= 2)' if Pe <= 2 else 'INSTÁVEL / OSCILATÓRIO no EDC (Pe > 2)'}")
        print("-" * 80)
        
        print(f"{'Nó':<5} {'x [m]':<8} {'T_Anal [°C]':<14} {'T_EDC [°C]':<14} {'T_Upwind [°C]':<14} {'Erro EDC [°C]':<14}")
        print("." * 75)
        for i in range(N):
            err_e = abs(T_edc[i] - T_anal[i])
            print(f"{i+1:<5} {x_nodes[i]:<8.4f} {T_anal[i]:<14.4f} {T_edc[i]:<14.4f} {T_upw[i]:<14.4f} {err_e:<14.4f}")
        print()

    # Plotagem dos resultados
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for idx, (item, u, N) in enumerate(cases):
        ax = axes[idx]
        x_nodes, T_edc, _, _, _, _, F, D, _, _ = solve_edc(u, N)
        _, T_upw, _, _, _, _, _, _, _, _ = solve_upwind(u, N)
        
        x_dense = np.linspace(0, 0.1, 200)
        T_an_dense = T_analytical(x_dense, u=u)
        
        ax.plot(x_dense, T_an_dense, 'k--', linewidth=2, label='Analítica Exata')
        
        x_plot = np.concatenate(([0.0], x_nodes, [0.1]))
        T_edc_plot = np.concatenate(([150.0], T_edc, [50.0]))
        T_upw_plot = np.concatenate(([150.0], T_upw, [50.0]))
        
        ax.plot(x_plot, T_edc_plot, 'ro-', linewidth=1.5, markersize=4, label='EDC (Diferenças Centrais)')
        ax.plot(x_plot, T_upw_plot, 'bs-', linewidth=1.5, markersize=4, label='Esquema Upwind')
        
        ax.set_xlabel('Posição x [m]')
        if idx == 0:
            ax.set_ylabel('Temperatura T [°C]')
        ax.set_title(f"{item} (u={u}m/s, N={N})")
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.show()
