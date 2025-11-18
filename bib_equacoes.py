import numpy as np


# Dimensionar engrenagem
def dimensionar_estagio(T_p_in, n_p_in, i_alvo, m, N_p, b_face_fator, eta_estagio):
    phi_graus = 20.0 # Angulo de pressao
    phi_rad = np.radians(phi_graus)

    N_c = round(N_p * i_alvo)
    i_efetiva = N_c / N_p

    d_p = m * N_p
    d_c = m * N_c
    r_p = d_p / 2
    r_c = d_c / 2
    b_face = m * b_face_fator
    C = (d_p + d_c) / 2 # Distancia entre Centros

    d_ext_p = m * (N_p + 2) # Diametro Externo Pinhao
    d_ext_c = m * (N_c + 2) # Diametro Externo Coroa
    r_ext_p = d_ext_p / 2
    r_ext_c = d_ext_c / 2

    # Razao de Contato
    p_circ = m * np.pi
    p_base = p_circ * np.cos(phi_rad)
    r_b_p = r_p * np.cos(phi_rad) # Raio de base pinhao
    r_b_c = r_c * np.cos(phi_rad) # Raio de base coroa

    # Comprimento de ação (Z)
    Z = (np.sqrt(r_ext_p**2 - r_b_p**2) +
        np.sqrt(r_ext_c**2 - r_b_c**2) -
         C * np.sin(phi_rad))

    mp = Z / p_base # Razao de contato

    # Forca Tangencial no pinhao
    W_t = (T_p_in * 1000) / (d_p / 2)
    # Forca Radial no pinhao
    W_r = W_t * np.tan(phi_rad)

    # Velocidade tangencial
    omega_p_in = n_p_in * (2 * np.pi / 60)
    V_t = (d_p / 1000) * omega_p_in / 2

    # Torque de saida e rotacao
    T_p_out = T_p_in * i_efetiva * eta_estagio
    n_p_out = n_p_in / i_efetiva

    # Aço AISI 4140 Normalizado
    S_at = 240  # MPa
    S_ac = 860  # MPa
    C_p = 191   # MPa^0.5

    # --- Fatores AGMA ---
    K_s = 1.0   # Fator de Tamanho
    K_R = 1.25  # Fator de Confiabilidade (Flexão, off-road)
    C_R = 1.25  # Fator de Confiabilidade (Superfície)
    K_T = 1.0   # Fator de Temperatura
    J = 0.32    # Fator de Geometria

    # Fator Dinâmico (K_v)
    Qv = 6
    B = np.pow(12 - Qv, 2/3) / 4
    A = 50 + 56 * (1 - B)
    K_v = (A / (A + np.sqrt(V_t * 200)))**B

    # Fator Distribuicao Carga (K_m)
    K_m = 0
    if b_face <= 50:
        K_m = 1.6
    elif 50 < b_face <= 150:
        K_m = 1.7
    elif 250 < b_face <= 500:
        K_m = 2.0
    else:
        K_m = 2.2 # Suposição para valores > 500

    # Analise de Flexao
    sigma_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
    sigma_b_adm = S_at / (K_T * K_R)

    FS_flexao = np.inf if sigma_b == 0 else sigma_b_adm / sigma_b

    # Analise de Contato/Pitting
    mg = N_c / N_p

    # Fator de Geometria (Superfície)
    I_geo = (np.sin(phi_rad) * np.cos(phi_rad) / 2) * (mg / (mg + 1))

    pitting = (W_t * K_s * K_m) / (b_face * d_p * I_geo * K_v)

    sigma_c = np.inf
    if pitting >= 0 and (b_face * d_p * I_geo * K_v) != 0:
        sigma_c = C_p * np.sqrt(pitting)

    sigma_c_adm = S_ac / C_R

    FS_pitting = 0.0
    if sigma_c == np.inf:
        FS_pitting = 0.0
    elif sigma_c == 0:
        FS_pitting = np.inf
    else:
        FS_pitting = (sigma_c_adm / sigma_c)**2

    return {
        "m": m,
        "N_p": N_p,
        "N_c": N_c,
        "b_face": b_face,
        "d_p": d_p,
        "d_c": d_c,
        "d_ext_p": d_ext_p,
        "d_ext_c": d_ext_c,
        "C": C,
        "Z": Z,
        "p_base": p_base,
        "mp": mp,
        "i_efetiva": i_efetiva,
        "FS_flexao": FS_flexao,
        "FS_pitting": FS_pitting,
        "W_t": W_t,
        "W_r": W_r,
        "T_p_out": T_p_out,
        "n_p_out": n_p_out
    }
