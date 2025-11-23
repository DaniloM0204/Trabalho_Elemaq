import numpy as np


# --- FUNÇÃO AUXILIAR PARA O FATOR J (AGMA) ---
def obtem_fator_J_AGMA_20graus(N_pinhao, N_coroa):
    """
    Retorna o Fator Geométrico J (AGMA 908-B89) para 20 graus.
    Carga na ponta (Conservador).
    """
    if N_pinhao <= 17:
        return 0.24
    elif N_pinhao <= 19:
        return 0.25
    elif N_pinhao <= 21:
        return 0.26
    elif N_pinhao <= 25:
        return 0.27
    elif N_pinhao <= 30:
        return 0.28
    else:
        return 0.29


# --- FUNÇÕES ORIGINAIS AJUSTADAS ---


def calcula_grupo_engrenagens(i_alvo, m, N_p, b_face_fator):
    phi_graus = 20.0
    phi_rad = np.radians(phi_graus)

    N_c = round(N_p * i_alvo)
    i_efetiva = N_c / N_p

    d_p = m * N_p
    d_c = m * N_c

    r_p = d_p / 2
    r_c = d_c / 2

    b_face = m * b_face_fator
    C = (d_p + d_c) / 2

    d_ext_p = m * (N_p + 2)
    d_ext_c = m * (N_c + 2)

    r_ext_p = d_ext_p / 2
    r_ext_c = d_ext_c / 2

    p_circ = m * np.pi
    p_base = p_circ * np.cos(phi_rad)

    r_b_p = r_p * np.cos(phi_rad)
    r_b_c = r_c * np.cos(phi_rad)

    Z = (
        np.sqrt(r_ext_p**2 - r_b_p**2)
        + np.sqrt(r_ext_c**2 - r_b_c**2)
        - C * np.sin(phi_rad)
    )

    mp = Z / p_base

    return {
        "m": m,
        "N_p": N_p,
        "N_c": N_c,
        "i_efetiva": i_efetiva,
        "d_p": d_p,
        "d_c": d_c,
        "d_ext_p": d_ext_p,
        "d_ext_c": d_ext_c,
        "C": C,
        "b_face": b_face,
        "mp": mp,
        "phi_rad": phi_rad,
    }


def calcula_forcas(T_p_in, n_p_in, grupoengre, eta_estagio):
    d_p = grupoengre["d_p"]
    phi_rad = grupoengre["phi_rad"]
    i_efetiva = grupoengre["i_efetiva"]

    # Se o torque vier em N.m (< 1000), converte para N.mm para o cálculo da força
    if T_p_in < 1000:
        T_mm = T_p_in * 1000
    else:
        T_mm = T_p_in

    # Wt = 2 * T / d
    W_t = (2 * T_mm) / d_p
    W_r = W_t * np.tan(phi_rad)

    omega_p_in = n_p_in * (2 * np.pi / 60)
    V_t = (d_p / 1000) * omega_p_in / 2

    # Mantém unidade de entrada para o torque de saída
    T_p_out = T_p_in * i_efetiva * eta_estagio
    n_p_out = n_p_in / i_efetiva

    return {
        "W_t": W_t,
        "W_r": W_r,
        "V_t": V_t,
        "T_p_out": T_p_out,
        "n_p_out": n_p_out,
        "i_efetiva": i_efetiva,
    }


def calcula_Km_AGMA(
    b_face, d_p, montagem_tipo="comercial", precisao_engrenagem="Q6"
):
    b_face_pol = b_face / 25.4
    d_p_pol = d_p / 25.4

    if b_face_pol <= 1.0:
        C_pf = (b_face_pol / (10 * d_p_pol)) - 0.025
    else:
        C_pf = (b_face_pol / (10 * d_p_pol)) - 0.0375 + 0.0125 * b_face_pol
    C_pf = max(C_pf, 0)

    if montagem_tipo == "precisao":
        C_pm = 1.0
    elif montagem_tipo == "comercial":
        C_pm = 1.1
    else:
        C_pm = 1.2

    if montagem_tipo == "precisao":
        A, B, C = 0.0675, 0.0128, -0.0000765
    elif montagem_tipo == "comercial":
        A, B, C = 0.127, 0.0158, -0.0000765
    else:
        A, B, C = 0.247, 0.0167, -0.0000765

    C_ma = A + (B * b_face_pol) + (C * b_face_pol**2)

    if montagem_tipo == "precisao":
        C_e = 0.8
    elif montagem_tipo == "comercial":
        C_e = 1.0
    else:
        C_e = 1.0

    K_m = 1.0 + C_pf * C_pm + C_ma * C_e
    return max(1.0, min(K_m, 2.5))


def calcula_Kv_AGMA(V_t, Qv, phi_rad):
    V_t_ftmin = V_t * 196.85
    B_coef = 0.25 * (12 - Qv) ** (2 / 3)
    A_coef = 50 + 56 * (1 - B_coef)

    if V_t_ftmin == 0:
        K_v = 1.0
    else:
        K_v = ((A_coef + np.sqrt(V_t_ftmin)) / A_coef) ** B_coef

    return max(1.0, min(K_v, 2.0))


def calcula_largura_face_automatica(
    m, d_p, W_t, S_at, K_m, K_v, J, FS_desejado=1.5
):
    """
    Calcula a largura de face.
    CORREÇÃO: Kv multiplicando no numerador.
    """
    sigma_b_adm = S_at / FS_desejado

    # Fórmula Lewis/AGMA corrigida: b = (Wt * Kv * Km) / (m * J * sigma_adm)
    b_min_flexao = (W_t * K_m * K_v) / (m * J * sigma_b_adm)

    b_min_aspecto = m * 8
    b_max_aspecto = m * 16
    b_max_diametro = 1.2 * d_p

    b_calc = max(b_min_flexao, b_min_aspecto)
    b_calc = min(b_calc, b_max_aspecto, b_max_diametro)

    b_face = round(b_calc / 5) * 5
    return max(b_face, 10)


def calcula_grupo_engrenagens_auto(
    i_alvo,
    m,
    N_p,
    W_t_estimado,
    S_at,
    K_m_estimado=1.5,
    K_v_estimado=1.2,
    J_estimado=0.24,
):
    phi_graus = 20.0
    phi_rad = np.radians(phi_graus)

    N_c = round(N_p * i_alvo)
    i_efetiva = N_c / N_p

    d_p = m * N_p
    d_c = m * N_c

    # J correto baseado na tabela
    J_real = obtem_fator_J_AGMA_20graus(N_p, N_c)

    b_face = calcula_largura_face_automatica(
        m=m,
        d_p=d_p,
        W_t=W_t_estimado,
        S_at=S_at,
        K_m=K_m_estimado,
        K_v=K_v_estimado,
        J=J_real,
    )

    C = (d_p + d_c) / 2
    d_ext_p = m * (N_p + 2)
    d_ext_c = m * (N_c + 2)

    p_circ = m * np.pi
    p_base = p_circ * np.cos(phi_rad)
    r_p, r_c = d_p / 2, d_c / 2
    r_b_p, r_b_c = r_p * np.cos(phi_rad), r_c * np.cos(phi_rad)
    r_ext_p, r_ext_c = d_ext_p / 2, d_ext_c / 2

    Z = (
        np.sqrt(r_ext_p**2 - r_b_p**2)
        + np.sqrt(r_ext_c**2 - r_b_c**2)
        - C * np.sin(phi_rad)
    )
    mp = Z / p_base

    return {
        "m": m,
        "N_p": N_p,
        "N_c": N_c,
        "i_efetiva": i_efetiva,
        "d_p": d_p,
        "d_c": d_c,
        "d_ext_p": d_ext_p,
        "d_ext_c": d_ext_c,
        "C": C,
        "b_face": b_face,
        "mp": mp,
        "phi_rad": phi_rad,
        "J_utilizado": J_real,
    }


def estima_W_t_inicial(T_in, d_p_estimado):
    # Converte para N.mm para estimar Wt corretamente
    if T_in < 1000:
        T_mm = T_in * 1000
    else:
        T_mm = T_in
    return (2 * T_mm) / d_p_estimado


def calcula_vida_util_engrenagens(
    grupoengre,
    forcas,
    sigma_b,
    sigma_c,
    S_at,
    S_ac,
    vida_util_desejada_horas=4000,
):
    n_p = forcas.get("n_p_in", 1000)
    ciclos_por_hora = n_p * 60

    N_base = 10**7
    b_flexao = -0.09
    c_pitting = -0.06

    # Vida Flexão
    if sigma_b > S_at:
        ratio = S_at / sigma_b
        if ratio < 0.1:
            N_flexao = 100
        else:
            N_flexao = N_base * (ratio ** (1 / abs(b_flexao)))
    else:
        N_flexao = 10**10

    vida_flexao_horas = N_flexao / ciclos_por_hora

    # Vida Pitting
    if sigma_c > S_ac:
        ratio = S_ac / sigma_c
        if ratio < 0.1:
            N_pitting = 100
        else:
            N_pitting = N_base * (ratio ** (1 / abs(c_pitting)))
    else:
        N_pitting = 10**10

    vida_pitting_horas = N_pitting / ciclos_por_hora
    fator_total = 0.5

    vida_min = min(vida_flexao_horas, vida_pitting_horas) * fator_total

    if vida_min >= vida_util_desejada_horas:
        status = "APROVADO"
    elif vida_min >= vida_util_desejada_horas * 0.7:
        status = "MARGINAL"
    else:
        status = "REPROVADO"

    # --- CORREÇÃO AQUI: Retornando as chaves que o utilities.py exige ---
    return {
        "vida_flexao_horas": vida_flexao_horas
        * fator_total,  # Adicionado de volta
        "vida_pitting_horas": vida_pitting_horas
        * fator_total,  # Adicionado de volta
        "vida_util_minima": vida_min,
        "status_vida": status,
        "vida_util_conservadora": vida_min,
    }


def analisa_fatores_com_vida(
    grupoengre,
    forcas,
    S_at,
    S_ac,
    C_p,
    vida_util_horas=4000,
    montagem_tipo="comercial",
    precisao_engrenagem="Q6",
):
    W_t = forcas["W_t"]
    V_t = forcas["V_t"]
    d_p = grupoengre["d_p"]
    b_face = grupoengre["b_face"]
    m = grupoengre["m"]
    N_p = grupoengre["N_p"]
    phi_rad = grupoengre["phi_rad"]
    J = grupoengre.get("J_utilizado", 0.24)

    K_v = calcula_Kv_AGMA(V_t, 6, phi_rad)
    K_m = calcula_Km_AGMA(b_face, d_p, montagem_tipo, precisao_engrenagem)

    K_o = 1.25
    K_s = 1.0
    K_B = 1.0
    K_I = 1.0

    sigma_b = (W_t * K_o * K_v * K_s * K_m * K_B * K_I) / (b_face * m * J)
    # Assumindo KR e KT
    K_R = 0.85  # Exemplo para 90% confiabilidade (Aula 08 Slide 46) ou 1.0 para 99%
    K_T = 1.0  # Se T < 120 graus C

    # Ajuste no cálculo do FS
    FS_flexao = abs((S_at) / (sigma_b * K_T * K_R))

    mg = grupoengre["N_c"] / grupoengre["N_p"]
    I_geo = (np.sin(phi_rad) * np.cos(phi_rad) / 2) * (mg / (mg + 1))

    term_pitting = (W_t * K_o * K_v * K_s * K_m) / (b_face * d_p * I_geo)
    sigma_c = C_p * np.sqrt(term_pitting)
    FS_pitting = (S_ac / sigma_c) ** 2

    vida_util = calcula_vida_util_engrenagens(
        grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_horas
    )

    status_completo = "REPROVADO"
    # Ajuste de critérios de aprovação
    if (
        FS_flexao >= 2.0
        and FS_pitting >= 1.2
        and vida_util["status_vida"] != "REPROVADO"
    ):
        status_completo = "APROVADO"
    elif FS_flexao >= 1.2 and FS_pitting >= 1.0:
        status_completo = "MARGINAL"

    return {
        "FS_flexao": FS_flexao,
        "FS_pitting": FS_pitting,
        "sigma_b": sigma_b,
        "sigma_c": sigma_c,
        "K_m": K_m,
        "K_v": K_v,
        "status_completo": status_completo,
        **vida_util,  # Isso desempacota todas as chaves, incluindo as que o utilities.py precisa
    }


def dimensionamento_engrenagens_com_vida(
    T_in, n_in, i_estagio, eta_estagio, S_at, S_ac, C_p, vida_util_horas=4000
):
    guardaengrenagem1 = []
    guardaengrenagem2 = []

    # Módulos estendidos para garantir que aguente a carga
    modulos_padronizados = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]
    N_pinhao_padronizado = [17, 18, 19, 20, 21, 22, 23, 24, 25]

    resultado_estagio_1 = None
    resultado_estagio_2 = None

    # ESTÁGIO 1
    for Np in N_pinhao_padronizado:
        for m in modulos_padronizados:
            d_p_estimado = m * Np
            W_t_estimado = estima_W_t_inicial(T_in, d_p_estimado)

            grupo_engrenagens = calcula_grupo_engrenagens_auto(
                i_estagio, m, Np, W_t_estimado, S_at
            )

            forcas = calcula_forcas(T_in, n_in, grupo_engrenagens, eta_estagio)
            forcas["n_p_in"] = n_in

            fatores = analisa_fatores_com_vida(
                grupo_engrenagens, forcas, S_at, S_ac, C_p, vida_util_horas
            )

            dados = {**grupo_engrenagens, **forcas, **fatores}
            dados["erro_estagio"] = abs(
                (dados["i_efetiva"] - i_estagio) / i_estagio
            )

            if dados["erro_estagio"] > 0.05:
                continue

            if dados["status_completo"] == "APROVADO":
                dados["status_calc"] = "Aprovado"
                guardaengrenagem1.append(dados)
            elif dados["status_completo"] == "MARGINAL":
                dados["status_calc"] = "Marginal"
                guardaengrenagem1.append(dados)

    aprovadas = [p for p in guardaengrenagem1 if p["status_calc"] == "Aprovado"]
    marginais = [p for p in guardaengrenagem1 if p["status_calc"] == "Marginal"]

    if aprovadas:
        aprovadas.sort(key=lambda x: (x["m"], x["erro_estagio"]))
        resultado_estagio_1 = aprovadas[0]
    elif marginais:
        marginais.sort(key=lambda x: (x["m"], x["erro_estagio"]))
        resultado_estagio_1 = marginais[0]
        print("AVISO: Estágio 1 usando configuração MARGINAL")
    else:
        print("FALHA CRÍTICA: Nenhuma engrenagem aprovada para Estágio 1.")
        return None, None, [], []

    # ESTÁGIO 2
    T_in_2 = resultado_estagio_1["T_p_out"]
    n_in_2 = resultado_estagio_1["n_p_out"]

    for Np in N_pinhao_padronizado:
        for m in modulos_padronizados:
            d_p_estimado = m * Np
            W_t_estimado = estima_W_t_inicial(T_in_2, d_p_estimado)

            grupo_engrenagens = calcula_grupo_engrenagens_auto(
                i_estagio, m, Np, W_t_estimado, S_at
            )

            forcas = calcula_forcas(
                T_in_2, n_in_2, grupo_engrenagens, eta_estagio
            )
            forcas["n_p_in"] = n_in_2

            fatores = analisa_fatores_com_vida(
                grupo_engrenagens, forcas, S_at, S_ac, C_p, vida_util_horas
            )

            dados = {**grupo_engrenagens, **forcas, **fatores}
            dados["erro_estagio"] = abs(
                (dados["i_efetiva"] - i_estagio) / i_estagio
            )

            if dados["erro_estagio"] > 0.05:
                continue

            if dados["status_completo"] == "APROVADO":
                dados["status_calc"] = "Aprovado"
                guardaengrenagem2.append(dados)
            elif dados["status_completo"] == "MARGINAL":
                dados["status_calc"] = "Marginal"
                guardaengrenagem2.append(dados)

    aprovadas2 = [
        p for p in guardaengrenagem2 if p["status_calc"] == "Aprovado"
    ]
    marginais2 = [
        p for p in guardaengrenagem2 if p["status_calc"] == "Marginal"
    ]

    if aprovadas2:
        aprovadas2.sort(key=lambda x: (x["m"], x["erro_estagio"]))
        resultado_estagio_2 = aprovadas2[0]
    elif marginais2:
        marginais2.sort(key=lambda x: (x["m"], x["erro_estagio"]))
        resultado_estagio_2 = marginais2[0]
        print("AVISO: Estágio 2 usando configuração MARGINAL")
    else:
        print("FALHA CRÍTICA: Nenhuma engrenagem aprovada para Estágio 2.")
        resultado_estagio_2 = None

    return (
        resultado_estagio_1,
        resultado_estagio_2,
        guardaengrenagem1,
        guardaengrenagem2,
    )
