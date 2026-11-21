import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# 1. ALGORITMO TDMA
# ============================================================

def solve_tdma(a_W, a_P, a_E, b_rhs):
    """
    Algoritmo TDMA para resolver:
    a_P[i]*T[i] = a_E[i]*T[i+1] + a_W[i]*T[i-1] + b_rhs[i]

    Relação de varredura:
    T[i] = P[i]*T[i+1] + Q[i]
    """

    n = len(a_P)
    P = np.zeros(n)
    Q = np.zeros(n)

    # Forward sweep
    P[0] = a_E[0] / a_P[0]
    Q[0] = b_rhs[0] / a_P[0]

    for i in range(1, n):
        denom = a_P[i] - a_W[i] * P[i-1]

        if i < n - 1:
            P[i] = a_E[i] / denom

        Q[i] = (b_rhs[i] + a_W[i] * Q[i-1]) / denom

    # Back substitution
    T = np.zeros(n)
    T[-1] = Q[-1]

    for i in range(n - 2, -1, -1):
        T[i] = P[i] * T[i+1] + Q[i]

    return T, P, Q


# ============================================================
# 2. SOLUÇÃO ANALÍTICA
# ============================================================

def T_analytical(x, t, L=0.1, alpha=2.5e-6,
                 T_init=150.0, n_terms=200):
    """
    Solução analítica exata em Série de Fourier.
    """

    if t <= 0:
        return np.full_like(x, T_init, dtype=float)

    T = np.zeros_like(x, dtype=float)

    for n in range(1, n_terms + 1):

        lam = (2*n - 1) * np.pi / (2*L)

        C_n = (
            4.0 * T_init * ((-1)**(n - 1))
            / ((2*n - 1) * np.pi)
        )

        T += (
            C_n
            * np.cos(lam * x)
            * np.exp(-alpha * (lam**2) * t)
        )

    return T


# ============================================================
# 3. SIMULAÇÃO MVF TRANSIENTE 1D
# ============================================================

def solve_fvm_transient_1d(
    N=30,
    L=0.1,
    k=25.0,
    rho_cp=10.0e6,
    T_init=150.0,
    T_L=0.0,
    t_final=360.0,
    dt=1.0
):
    """
    Simulação 1D transiente via MVF totalmente implícito.

    Fronteira esquerda (x=0): adiabática.
    Fronteira direita (x=L): temperatura prescrita.
    
    Armazena a solução numérica e analítica a cada segundo.
    """

    dx = L / N

    x_nodes = np.array([
        (i + 0.5) * dx for i in range(N)
    ])

    a_P0 = rho_cp * dx / dt

    # Coeficientes da matriz
    a_W = np.zeros(N)
    a_E = np.zeros(N)
    a_P = np.zeros(N)

    for i in range(N):

        if i == 0:

            a_W[i] = 0.0
            a_E[i] = k / dx
            a_P[i] = a_E[i] + a_P0

        elif i == N - 1:

            a_W[i] = k / dx

            a_E_boundary = k / (dx / 2.0)

            a_E[i] = 0.0
            a_P[i] = a_W[i] + a_E_boundary + a_P0

        else:

            a_W[i] = k / dx
            a_E[i] = k / dx
            a_P[i] = a_W[i] + a_E[i] + a_P0

    # Condição inicial
    T_num = np.full(N, T_init, dtype=float)

    # Armazenamento de todos os instantes
    results_num = {0: T_num.copy()}

    results_anal = {
        0: T_analytical(x_nodes, 0.0)
    }

    max_errors = {0: 0.0}
    mean_errors = {0: 0.0}

    # Condutância da fronteira direita
    a_E_boundary = k / (dx / 2.0)

    num_steps = int(round(t_final / dt))

    for step in range(1, num_steps + 1):

        t_curr = step * dt

        # Termo fonte transiente
        b_rhs = a_P0 * T_num.copy()

        # Condição de contorno T = T_L
        b_rhs[-1] += a_E_boundary * T_L

        # Solução numérica TDMA
        T_num, P, Q = solve_tdma(
            a_W, a_P, a_E, b_rhs
        )

        # Solução analítica no mesmo instante
        T_an = T_analytical(
            x_nodes,
            t_curr,
            L=L,
            alpha=k / rho_cp,
            T_init=T_init
        )

        # Erro absoluto em todos os volumes
        err_abs = np.abs(T_num - T_an)

        # Armazenamento segundo a segundo
        time_key = int(round(t_curr))

        results_num[time_key] = T_num.copy()
        results_anal[time_key] = T_an.copy()

        max_errors[time_key] = np.max(err_abs)
        mean_errors[time_key] = np.mean(err_abs)

    return (
        x_nodes,
        results_num,
        results_anal,
        max_errors,
        mean_errors,
        dx
    )


# ============================================================
# 4. EXECUÇÃO PRINCIPAL
# ============================================================

if __name__ == "__main__":

    # Diretórios locais do projeto no Windows
    BASE_DIR = Path(__file__).resolve().parent

    OUT_DIR = BASE_DIR / "out"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Parâmetros físicos e numéricos
    N = 30
    L = 0.1
    k = 25.0
    rho_cp = 10.0e6
    T_init = 150.0
    T_L = 0.0
    t_final = 360
    dt = 1.0

    # Execução da simulação
    (
        x_nodes,
        res_num,
        res_anal,
        max_errs,
        mean_errs,
        dx
    ) = solve_fvm_transient_1d(
        N=N,
        L=L,
        k=k,
        rho_cp=rho_cp,
        T_init=T_init,
        T_L=T_L,
        t_final=t_final,
        dt=dt
    )

    # ========================================================
    # 5. ESTATÍSTICAS GLOBAIS
    # ========================================================

    seconds_arr = np.arange(1, t_final + 1)

    max_err_array = np.array([
        max_errs[t] for t in seconds_arr
    ])

    mean_err_array = np.array([
        mean_errs[t] for t in seconds_arr
    ])

    print("=" * 100)
    print("LISTA DE EXERCÍCIOS 6 - CONDUÇÃO TRANSIENTE 1D")
    print("MÉTODO DOS VOLUMES FINITOS + TDMA")
    print("=" * 100)

    print(f"Malha: N = {N} volumes")
    print(f"Comprimento: L = {L} m")
    print(f"Espaçamento: dx = {dx:.6f} m")
    print(f"Passo temporal: dt = {dt:.2f} s")
    print(f"Instantes avaliados: {t_final + 1}")

    t_peak_err = seconds_arr[np.argmax(max_err_array)]

    print("\nRESUMO DOS ERROS:")
    print(f"Maior erro máximo: {np.max(max_err_array):.6f} °C")
    print(f"Instante do maior erro: {t_peak_err} s")
    print(f"Erro máximo em t = 360 s: {max_errs[t_final]:.6f} °C")
    print(f"Erro médio em t = 360 s: {mean_errs[t_final]:.6f} °C")
    print(f"Erro médio global: {np.mean(mean_err_array):.6f} °C")

    # ========================================================
    # 6. TABELA DE COMPARAÇÃO SEGUNDO A SEGUNDO
    # ========================================================

    idx_mid = N // 2

    print("\n")
    print("=" * 110)
    print("COMPARAÇÃO NUMÉRICA E ANALÍTICA - TODOS OS SEGUNDOS")
    print("=" * 110)

    header = (
        f"{'Tempo [s]':>10}"
        f"{'T_num(0)':>14}"
        f"{'T_anal(0)':>14}"
        f"{'T_num(L/2)':>14}"
        f"{'T_anal(L/2)':>14}"
        f"{'Erro Máx.':>14}"
        f"{'Erro Médio':>14}"
    )

    print(header)
    print("-" * 110)

    for t in range(0, t_final + 1):

        Tn_0 = res_num[t][0]
        Ta_0 = res_anal[t][0]

        Tn_mid = res_num[t][idx_mid]
        Ta_mid = res_anal[t][idx_mid]

        err_max = max_errs[t]
        err_mean = mean_errs[t]

        print(
            f"{t:>10d}"
            f"{Tn_0:>14.4f}"
            f"{Ta_0:>14.4f}"
            f"{Tn_mid:>14.4f}"
            f"{Ta_mid:>14.4f}"
            f"{err_max:>14.4f}"
            f"{err_mean:>14.4f}"
        )

    print("=" * 110)


    # ========================================================
    # 7. GRÁFICOS
    # ========================================================

    plt.style.use("default")

    times_plot = [0, 10, 30, 60, 120, 240, 360]

    # --------------------------------------------------------
    # GRÁFICO 1 - EVOLUÇÃO TEMPORAL DA TEMPERATURA NA BARRA
    # --------------------------------------------------------

    fig1, ax1 = plt.subplots(figsize=(9, 6))

    colors = plt.cm.tab10.colors

    for idx, t in enumerate(times_plot):

        # Solução numérica - MVF
        ax1.plot(
            x_nodes,
            res_num[t],
            "o-",
            color=colors[idx],
            linewidth=1.5,
            markersize=4,
            label=f"MVF - {t} s"
        )

        # Solução analítica
        ax1.plot(
            x_nodes,
            res_anal[t],
            "--",
            color=colors[idx],
            linewidth=1.5,
            label=f"Analítica - {t} s"
        )

    ax1.set_xlabel("Posição na barra, x [m]", fontsize=11)
    ax1.set_ylabel("Temperatura, T [°C]", fontsize=11)

    ax1.set_title(
        "Evolução Temporal da Temperatura na Barra\n"
        "Comparação MVF e Solução Analítica",
        fontsize=13
    )

    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8, ncol=2)

    fig1.tight_layout()

    fig1.savefig(
        OUT_DIR / "evolucao_temperatura_barra.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig1)


    # --------------------------------------------------------
    # GRÁFICO 2 - EVOLUÇÃO TEMPORAL DO ERRO ABSOLUTO
    # --------------------------------------------------------

    fig2, ax2 = plt.subplots(figsize=(9, 6))

    ax2.plot(
        seconds_arr,
        max_err_array,
        color="red",
        linewidth=1.8,
        label="Erro máximo"
    )

    ax2.plot(
        seconds_arr,
        mean_err_array,
        color="blue",
        linestyle="--",
        linewidth=1.8,
        label="Erro médio"
    )

    ax2.set_xlabel("Tempo, t [s]", fontsize=11)
    ax2.set_ylabel("Erro absoluto [°C]", fontsize=11)

    ax2.set_title(
        "Evolução Temporal do Erro Absoluto\n"
        "Comparação entre MVF e Solução Analítica",
        fontsize=13
    )

    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)

    fig2.tight_layout()

    fig2.savefig(
        OUT_DIR / "evolucao_erro_absoluto.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close(fig2)

    print("\nGráficos salvos com sucesso!")
    print("1. evolucao_temperatura_barra.png")
    print("2. evolucao_erro_absoluto.png")
