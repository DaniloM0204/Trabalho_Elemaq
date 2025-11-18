import numpy as np


def calcula_grupo_engrenagens(i_alvo, m, N_p, b_face_fator):
    phi_graus = 20.0
    # Passando para radianos pois os calculos sao feitos em radianos
    phi_rad = np.radians(phi_graus)

    N_c = round(N_p * i_alvo) # Calculo do numero de dentes da coroa, sendo o numero arredondado pro inteiro mais próximo
    i_efetiva = N_c / N_p

    # Calculo dos diametros primitivos
    d_p = m * N_p
    d_c = m * N_c

    # Calculo dos raios primitivos
    r_p = d_p / 2
    r_c = d_c / 2

    # Utilizacao das larguras das faces para o cálculo dos fatores
    b_face = m * b_face_fator
    C = (d_p + d_c) / 2 # Distancia entre Centros

    # Diametros externos do pinhão e da coroa
    d_ext_p = m * (N_p + 2)
    d_ext_c = m * (N_c + 2)

    # Raios externos do pinhão e da coroa
    r_ext_p = d_ext_p / 2
    r_ext_c = d_ext_c / 2

    """Aqui é feito o calculo da Razao de Contato das engrenagens"""
    p_circ = m * np.pi

    # Diametro de base
    p_base = p_circ * np.cos(phi_rad)

    # Raios de base do par
    r_b_p = r_p * np.cos(phi_rad)
    r_b_c = r_c * np.cos(phi_rad)

    # Comprimento de ação que é utilizado para calcular a razão de contato
    Z = (np.sqrt(r_ext_p**2 - r_b_p**2) +
        np.sqrt(r_ext_c**2 - r_b_c**2) -
         C * np.sin(phi_rad))

    mp = Z / p_base # Razao de contato

    return {
        "m": m, "N_p": N_p, "N_c": N_c, "i_efetiva": i_efetiva,
        "d_p": d_p, "d_c": d_c, "d_ext_p": d_ext_p, "d_ext_c": d_ext_c,
        "C": C, "b_face": b_face, "mp": mp, "phi_rad": phi_rad
    }



# Dimensionar engrenagem
def calcula_forcas(T_p_in,n_p_in,grupoengre,eta_estagio):
    # Dados dos pares
    d_p = grupoengre['d_p']
    phi_rad = grupoengre['phi_rad']
    i_efetiva = grupoengre['i_efetiva']


    """A partir dos Torques de entrada, velocidades e etc, sao calculados todos os esforcos nas engrenagens que serao passados para as analises de flexao, e o desgates dos dentes (Pitting)"""

    # Forca Tangencial no pinhao
    W_t = (T_p_in * 1000) / (d_p / 2)
    # Forca Radial no pinhao
    W_r = W_t * np.tan(phi_rad)

    # Velocidade tangencial
    omega_p_in = n_p_in * (2 * np.pi / 60) # Em rad/s
    V_t = (d_p / 1000) * omega_p_in / 2 # Em m/s

    # Torque de saida e rotacao
    T_p_out = T_p_in * i_efetiva * eta_estagio
    n_p_out = n_p_in / i_efetiva

    return {
        "W_t": W_t, "W_r": W_r, "V_t": V_t,
        "T_p_out": T_p_out, "n_p_out": n_p_out, "i_efetiva": i_efetiva
            }


def analisa_fatores(grupoengre,forcas,S_at,S_ac,C_p):
    """Aqui eu peguei das tabelas 12-(14,16,17,18) do Norton os fatores pra simplificar as contas e tentar condizer com o jeito que as engrenagens sao fabricadas e pra que, entao tem Kr que escolhi pra carro offroad, o fator de geometria e oque representa o quando o formato do dente influencia nas tensoes e etc"""

    # Recuperando dados
    W_t = forcas['W_t']
    V_t = forcas['V_t']
    d_p = grupoengre['d_p']
    b_face = grupoengre['b_face']
    m = grupoengre['m']
    N_p,N_c = grupoengre['N_p'],grupoengre['N_c']
    phi_rad = grupoengre['phi_rad']

    # --- Fatores AGMA ---
    K_s = 1.0   # Fator de Tamanho
    K_R = 1.25  # Fator de Confiabilidade (Flexão, off-road)
    C_R = 1.25  # Fator de Confiabilidade (Superfície)
    K_T = 1.0   # Fator de Temperatura
    J = 0.32    # Fator de Geometria

    # O fator dinamico e calculado por essa equação que no livro ta na figura 12.22
    # Fator Dinâmico (K_v)
    Qv = 6
    B = np.pow(12 - Qv, 2/3) / 4
    A = 50 + 56 * (1 - B)
    K_v = (A / (A + np.sqrt(V_t * 200)))**B

    # Fator Distribuicao Carga (K_m) representa a distribuicao da carga ao longo da largura da face
    K_m = 0
    if b_face <= 50:
        K_m = 1.6
    elif 50 < b_face <= 150:
        K_m = 1.7
    elif 250 < b_face <= 500:
        K_m = 2.0
    else:
        K_m = 2.2 # Para largura de face maior que 500mm

    """A partir de todas as contas acima, sai feitas as analises de flexao e pitting, fazendo aqui o calculo do fator de seguranca pra iteracao das engrenagens"""


    # Analise de Flexao
    sigma_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
    sigma_b_adm = S_at / (K_T * K_R)

    FS_flexao = np.inf if sigma_b == 0 else sigma_b_adm / sigma_b

    # Analise de Contato/Pitting
    mg = N_c / N_p

    # Fator de Geometria (Superfície)
    I_geo = (np.sin(phi_rad) * np.cos(phi_rad) / 2) * (mg / (mg + 1))

    # Desgaste dos dentes das engrenagens
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

    return {"FS_flexao": FS_flexao, "FS_pitting": FS_pitting}
