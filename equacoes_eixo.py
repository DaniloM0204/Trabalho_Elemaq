import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def calcula_c_superf(S_ut):
    """
    Analisa TODOS os acabamentos e retorna o MAIOR fator (Melhor cenário).
    """
    A_b = [
        (1.58, -0.085),  # Retificado
        (4.51, -0.265),  # Usinado
        (57.7, -0.718),  # Laminado
        (272, -0.995),
    ]  # Forjado

    nomes = ["retificado", "usinado", "laminado", "forjado"]
    c_lista = []

    for i in range(len(A_b)):
        A = A_b[i][0]
        b = A_b[i][1]
        fator = A * (S_ut**b)
        if fator > 1.0:
            fator = 1.0
        c_lista.append(fator)

    maior_valor = max(c_lista)
    indice_maior = c_lista.index(maior_valor)
    melhor_acabamento = nomes[indice_maior]

    return maior_valor, melhor_acabamento


def calcula_fator_tamanho(d):
    """Calcula C_tamanho segundo Norton para Flexão/Torção"""
    if d <= 8:
        return 1.0
    elif d <= 250:
        return 1.189 * (d**-0.097)
    else:
        return 0.6


def dimensiona_eixo_por_fadiga(
    M_max,
    Tm,
    Sut,
    Sy,
    Se_linha,
    vida_util_horas=4000,
    Nf=2.0,
    tipo_eixo="simples",
):
    """
    Dimensiona eixo utilizando o critério de Goodman Modificado (Fadiga)
    e o critério de Langer (Escoamento Estático).
    """
    # 1. CONVERSÃO DE UNIDADES
    M_max_Nmm = M_max  # Já em N.mm
    Tm_Nmm = Tm * 1000 # N.m -> N.mm

    # 2. DEFINIÇÃO DOS FATORES DE CONCENTRAÇÃO DE TENSÃO (Kt)
    if tipo_eixo == "simples":
        Kt_flexao = 2.0
        Kt_torcao = 1.7
    else: # "duplo" ou complexo
        Kt_flexao = 2.2
        Kt_torcao = 1.8

    # 3. SENSIBILIDADE AO ENTALHE (q)
    # Para Aço AISI 1020 (Sut ~470 MPa), q é aprox 0.75
    q = 0.75

    # Fatores Dinâmicos (Kf)
    Kf_flexao = 1 + q * (Kt_flexao - 1)
    Kf_torcao = 1 + q * (Kt_torcao - 1)

    # 4. FATORES DE MARIN (CONSTANTES)
    C_carreg = 1.0  # Flexão Rotativa
    C_sup = 0.88    # Usinado (Valor médio para 1020)
    C_temp = 1.0    # Temp. Ambiente
    C_conf = 0.814  # 99% Confiabilidade

    # 5. LOOP DE ITERAÇÃO
    d_guess = 20.0 # Chute inicial (mm)
    tolerance = 0.01
    max_iter = 100
    iter_count = 0

    print(f"Dimensionando {tipo_eixo} (M={M_max:.0f}, T={Tm_Nmm:.0f})...")

    while iter_count < max_iter:
        # A. Fator de Tamanho (Kb) - Depende do diâmetro!
        if d_guess <= 8:
            C_tam = 1.0
        elif d_guess <= 250:
            C_tam = 1.189 * (d_guess ** -0.097)
        else:
            C_tam = 0.6

        # B. Limite de Resistência à Fadiga Corrigido (Se)
        Se = Se_linha * C_carreg * C_tam * C_sup * C_temp * C_conf

        # C. Tensões Atuantes (Baseadas no diâmetro atual)
        # Tensão Alternada (Flexão)
        sigma_a = Kf_flexao * (32 * M_max_Nmm) / (np.pi * d_guess**3)

        # Tensão Média (Torção)
        # Kfm = 1.0 para materiais dúcteis (não multiplica por Kf_torcao)
        tau_m = (16 * Tm_Nmm) / (np.pi * d_guess**3)

        # D. Critério 1: Goodman Modificado (Fadiga)
        # (sigma_a / Se) + (tau_m / Sut) = 1 / Nf
        utilizacao_fadiga = (sigma_a / Se) + (tau_m / Sut)

        # E. Critério 2: Langer (Escoamento Estático)
        # Von Mises das tensões nominais máximas
        sigma_max = (32 * M_max_Nmm) / (np.pi * d_guess**3)
        tau_max = (16 * Tm_Nmm) / (np.pi * d_guess**3)
        sigma_von_mises = np.sqrt(sigma_max**2 + 3 * tau_max**2)

        utilizacao_escoamento = sigma_von_mises / Sy

        # F. Convergência
        utilizacao_maxima = max(utilizacao_fadiga, utilizacao_escoamento)
        objetivo = 1.0 / Nf

        if abs(utilizacao_maxima - objetivo) < (tolerance * objetivo):
            break # Achou o diâmetro ideal!

        # Ajuste do diâmetro para a próxima tentativa
        if utilizacao_maxima > objetivo:
            d_guess *= 1.01 # Aumenta devagar (1%)
        else:
            d_guess *= 0.99 # Diminui

        iter_count += 1

    # 6. RESULTADOS FINAIS PARA O RELATÓRIO
    FS_fadiga = 1 / ((sigma_a / Se) + (tau_m / Sut))
    FS_escoamento = Sy / sigma_von_mises

    return {
        "diametro_minimo": d_guess,
        "FS_fadiga": FS_fadiga,
        "FS_escoamento": FS_escoamento,
        "Se_corrigido": Se,
        "sigma_a": sigma_a,
        "tau_m": tau_m,
        "Kf_flexao": Kf_flexao,
        "C_tam": C_tam,
        "iteracoes": iter_count,
    }


def analisa_deflexao_eixo(M_max, comprimento, diametro, material="aco"):
    E = 200000
    Inercia = (np.pi * diametro**4) / 64
    deflexao_max = (M_max * comprimento**2) / (10 * E * Inercia)
    deflexao_admissivel = 0.0005 * comprimento
    status_deflexao = (
        "APROVADO" if deflexao_max <= deflexao_admissivel else "REPROVADO"
    )

    return {
        "deflexao_max_mm": deflexao_max,
        "deflexao_admissivel_mm": deflexao_admissivel,
        "status_deflexao": status_deflexao,
    }


def calcula_velocidade_critica(diametro, comprimento):
    E = 200000
    rho = 7850
    diametro_m = diametro / 1000
    Inercia = (np.pi * diametro_m**4) / 64
    A = (np.pi * diametro_m**2) / 4
    L = comprimento / 1000
    m = rho * A
    omega_n = (np.pi**2 / L**2) * np.sqrt(E * 1e6 * Inercia / m)
    N_critica = omega_n * 60 / (2 * np.pi)
    return {"velocidade_critica_rpm": N_critica}


# --- FUNÇÕES DE CÁLCULO DE ESFORÇOS ---
def calcula_esforcos_analitico(L_total, posicoes, forcas, torque, nome_eixo):
    """
    Calcula Reações e Momento Máximo para qualquer configuração de eixo.
    Entradas:
      posicoes: dict {'A': 0, 'B': 100, ...}
      forcas: list [(pos, Fv, Fh), ...]
    Retorna: Dicionário com resultados escalares para log.
    """
    xA = posicoes["A"]
    xB = posicoes["B"]
    dist_mancais = xB - xA

    # 1. CÁLCULO DAS REAÇÕES (Estática)
    # Somatório de Momentos em A para achar Rb
    soma_M_A_v = 0
    soma_M_A_h = 0
    soma_F_v = 0
    soma_F_h = 0

    for pos, Fv, Fh in forcas:
        soma_M_A_v += Fv * (pos - xA)
        soma_M_A_h += Fh * (pos - xA)
        soma_F_v += Fv
        soma_F_h += Fh

    # Reações em B (Sinal oposto ao momento gerado pelas cargas)
    R_Bv = -soma_M_A_v / dist_mancais
    R_Bh = -soma_M_A_h / dist_mancais

    # Reações em A (Soma das forças = 0)
    R_Av = -soma_F_v - R_Bv
    R_Ah = -soma_F_h - R_Bh

    R_A_res = np.sqrt(R_Av**2 + R_Ah**2)
    R_B_res = np.sqrt(R_Bv**2 + R_Bh**2)

    # 2. ENCONTRAR MOMENTO MÁXIMO (Varredura)
    # Criamos a lista completa de cargas (Externas + Reações)
    todas_cargas = forcas.copy()
    todas_cargas.append((xA, R_Av, R_Ah))
    todas_cargas.append((xB, R_Bv, R_Bh))
    todas_cargas.sort(key=lambda x: x[0])  # Ordenar por posição

    # Vamos calcular o momento nos pontos de aplicação de força (onde ocorrem os picos)
    # e em alguns pontos intermediários para garantir
    pontos_interesse = sorted(
        list(set([p[0] for p in todas_cargas] + [xA, xB]))
    )
    max_M_res = 0
    pos_max = 0

    for x in pontos_interesse:
        Mv, Mh = 0, 0
        # Método das seções (olhando para esquerda de x)
        for pos, Fv, Fh in todas_cargas:
            if pos <= x:  # Se a força está à esquerda
                braço = x - pos
                Mv += Fv * braço
                Mh += Fh * braço

        M_res = np.sqrt(Mv**2 + Mh**2)
        if M_res > max_M_res:
            max_M_res = M_res
            pos_max = x

    return {
        "M_max": max_M_res,
        "posicao_M_max": pos_max,
        "R_A_resultante": R_A_res,
        "R_B_resultante": R_B_res,
        "R_Av": R_Av,
        "R_Ah": R_Ah,
        "R_Bv": R_Bv,
        "R_Bh": R_Bh,
    }


# --- FUNÇÃO GRÁFICA NOVA (ILUSTRAÇÃO + DIAGRAMAS) ---
def plotar_diagramas_completo(
    L_total, posicoes, forcas, torque, nome_eixo, d_estimado=30
):
    """
    posicoes: dict {'A': 0, 'B': 100, 'Eng1': 40...}
    forcas: list of tuples (pos, F_vertical, F_horizontal)
    d_estimado: diâmetro visual para o desenho
    """
    # 1. Resolver Reações (Genérico para qualquer viga bi-apoiada)
    xA, xB = posicoes["A"], posicoes["B"]

    # Soma de momentos em A (considerando sentido positivo horario)
    # Momento = F * braco. Se F for positiva (pra cima), momento anti-horario (-)
    # Vamos usar convenção: Força pra cima (+), Momento anti-horario (+)
    soma_M_A_v = 0
    soma_M_A_h = 0
    soma_F_v = 0
    soma_F_h = 0

    for pos, Fv, Fh in forcas:
        soma_M_A_v += Fv * (pos - xA)
        soma_M_A_h += Fh * (pos - xA)
        soma_F_v += Fv
        soma_F_h += Fh

    # Reações em B
    R_Bv = -soma_M_A_v / (xB - xA)  # Sinal negativo para contrabalancear
    R_Bh = -soma_M_A_h / (xB - xA)

    # Reações em A
    R_Av = -soma_F_v - R_Bv
    R_Ah = -soma_F_h - R_Bh

    # Lista para plotagem (Cargas + Reações)
    cargas_plot = forcas.copy()
    cargas_plot.append((xA, R_Av, R_Ah))
    cargas_plot.append((xB, R_Bv, R_Bh))
    cargas_plot.sort(key=lambda x: x[0])

    # 2. Calcular vetores V e M
    x = np.linspace(0, L_total, 1000)
    V_res = np.zeros_like(x)
    M_res = np.zeros_like(x)
    T_res = np.zeros_like(x)

    # Torque
    elementos_torque = sorted(
        [v for k, v in posicoes.items() if k not in ["A", "B"]]
    )
    if len(elementos_torque) == 1:  # Entrada ou Saida
        if (
            "Entrada" in nome_eixo or "1" in nome_eixo
        ):  # Torque entra em A e sai na engrenagem
            mask = (x >= xA) & (x <= elementos_torque[0])
        else:  # Torque entra na engrenagem e sai em B (ou vice versa, simplificação: trecho carregado)
            mask = (x >= min(xA, elementos_torque[0])) & (
                x <= max(xB, elementos_torque[0])
            )  # Ajustar conforme lógica real
            # Melhor lógica: Torque flui da Engrenagem para a Saida do eixo.
            # Vou assumir fluxo padrão: Eng -> Ponta ou Ponta -> Eng
            mask = x <= elementos_torque[0]  # Simplificação visual
    elif len(elementos_torque) == 2:  # Intermediário
        mask = (x >= elementos_torque[0]) & (x <= elementos_torque[1])
    else:
        mask = np.zeros_like(x, dtype=bool)

    T_res[mask] = torque

    for i, xi in enumerate(x):
        Vv_local, Vh_local = 0, 0
        Mv_local, Mh_local = 0, 0

        for pos, Fv, Fh in cargas_plot:
            if pos < xi:  # Método das seções olhando para a esquerda
                Vv_local += Fv
                Vh_local += Fh
                Mv_local += Fv * (xi - pos)
                Mh_local += Fh * (xi - pos)

        V_res[i] = np.sqrt(Vv_local**2 + Vh_local**2)
        M_res[i] = np.sqrt(Mv_local**2 + Mh_local**2)

    M_max = np.max(M_res)

    # 3. PLOTAGEM (4 Subplots: Esquema, V, M, T)
    fig, (ax0, ax1, ax2, ax3) = plt.subplots(
        4,
        1,
        figsize=(10, 14),
        gridspec_kw={"height_ratios": [1, 1, 1, 1]},
        sharex=True,
    )

    # --- AX0: ESQUEMA DO EIXO (Ilustração) ---
    ax0.set_title(
        f"{nome_eixo} - Esquema Geométrico", fontsize=12, fontweight="bold"
    )
    ax0.set_ylim(-d_estimado * 1.5, d_estimado * 1.5)
    ax0.axis("off")  # Desliga eixos numéricos

    # Desenha Eixo
    rect_eixo = patches.Rectangle(
        (0, -d_estimado / 2),
        L_total,
        d_estimado,
        linewidth=1,
        edgecolor="black",
        facecolor="lightgray",
    )
    ax0.add_patch(rect_eixo)
    # Linha de centro
    ax0.axhline(0, color="black", linestyle="-.", linewidth=0.5)

    # Desenha Mancais (A e B)
    w_mancal = 15
    h_mancal = d_estimado + 10
    for nome, pos in posicoes.items():
        if nome in ["A", "B"]:
            rect_mancal = patches.Rectangle(
                (pos - w_mancal / 2, -h_mancal / 2),
                w_mancal,
                h_mancal,
                linewidth=1,
                edgecolor="black",
                facecolor="lightblue",
                hatch="///",
            )
            ax0.add_patch(rect_mancal)
            ax0.text(pos, h_mancal / 2 + 2, nome, ha="center")

    # Desenha Engrenagens
    for nome, pos in posicoes.items():
        if nome not in ["A", "B"]:
            # Pinhão é menor, Coroa é maior (Visual)
            if "Pinhao" in nome:
                h_eng = d_estimado * 2.5
            else:
                h_eng = d_estimado * 4.0
            w_eng = 20  # Largura visual

            rect_eng = patches.Rectangle(
                (pos - w_eng / 2, -h_eng / 2),
                w_eng,
                h_eng,
                linewidth=1,
                edgecolor="black",
                facecolor="orange",
                alpha=0.7,
            )
            ax0.add_patch(rect_eng)
            ax0.text(pos, -h_eng / 2 - 5, nome, ha="center")

    # --- AX1: CORTANTE ---
    ax1.plot(x, V_res, color="blue", linewidth=2)
    ax1.fill_between(x, V_res, color="blue", alpha=0.1)
    ax1.set_ylabel("Cortante Resultante (N)")
    ax1.grid(True, alpha=0.3)

    # --- AX2: MOMENTO ---
    ax2.plot(x, M_res, color="red", linewidth=2)
    ax2.fill_between(x, M_res, color="red", alpha=0.1)
    ax2.set_ylabel("Momento Resultante (N.mm)")
    ax2.grid(True, alpha=0.3)
    ax2.annotate(
        f"Mmax: {M_max:.0f}",
        xy=(x[np.argmax(M_res)], M_max),
        xytext=(10, 5),
        textcoords="offset points",
    )

    # --- AX3: TORQUE ---
    ax3.plot(x, T_res / 1000, color="purple", linewidth=2)
    ax3.fill_between(x, T_res / 1000, color="purple", alpha=0.1)
    ax3.set_ylabel("Torque (N.m)")
    ax3.set_xlabel("Comprimento do Eixo (mm)")
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    import os

    if not os.path.exists("Plots"):
        os.makedirs("Plots")
    plt.savefig(f"Plots/Diagrama_{nome_eixo.replace(' ', '_')}.png", dpi=100)
    plt.close()

    return {
        "M_max": M_max,
        "R_A": np.sqrt(R_Av**2 + R_Ah**2),
        "R_B": np.sqrt(R_Bv**2 + R_Bh**2),
    }


def calcula_velocidade_critica_rayleigh(d_mm, L_mm):
    """
    Calcula a primeira velocidade crítica (Rayleigh) para eixo uniforme em aço.
    Entrada: d (mm), L (mm)
    Saída: RPM crítica
    """
    # Constantes do Aço
    E = 200000 * 1e6  # Pa (200 GPa)
    rho = 7850        # kg/m^3

    # Conversão para SI (metros)
    d = d_mm / 1000.0
    L = L_mm / 1000.0

    # Propriedades Geométricas
    area = (np.pi * d**2) / 4
    inercia = (np.pi * d**4) / 64
    massa_linear = rho * area  # kg/m

    # Frequência natural (rad/s) - Modelo Viga Bi-Apoiada
    # wn = (pi^2 / L^2) * sqrt(EI / m_linear)
    wn = (np.pi**2 / L**2) * np.sqrt((E * inercia) / massa_linear)

    # Conversão para RPM
    n_critica = (wn * 60) / (2 * np.pi)

    return n_critica
