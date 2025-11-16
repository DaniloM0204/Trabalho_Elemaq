import numpy as np

#Dados de Entrada
P_in = 1580      # Potência de entrada
n_in = 1450   # Rotacao de entrada
F_cabo = 1988    # Forca de tracao
D_tambor = 0.12  # Diametro do tambor
η_total = 0.85   # Eficiencia total do redutor

# Parametros iniciais
ω_in = n_in * (2 * np.pi / 60)
T_in = P_in / ω_in
T_out = F_cabo * (D_tambor / 2)
i_total = T_out / (T_in * η_total)
i_alvo_1 = np.sqrt(i_total)
i_alvo_2 = np.sqrt(i_total)
η_1 = np.sqrt(η_total)
η_2 = np.sqrt(η_total)


def dimensionar_estagio(T_p_in, n_p_in, i_alvo, m, N_p, b_face_fator, eta_estagio):
    φ_graus = 20.0 # Angulo de pressao
    φ_rad = np.radians(φ_graus)

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
    p_base = p_circ * np.cos(φ_rad)
    r_b_p = r_p * np.cos(φ_rad) # Raio de base pinhao
    r_b_c = r_c * np.cos(φ_rad) # Raio de base coroa

    # Comprimento de ação (Z)
    Z = (np.sqrt(r_ext_p**2 - r_b_p**2) +
        np.sqrt(r_ext_c**2 - r_b_c**2) -
         C * np.sin(φ_rad))

    mp = Z / p_base # Razão de contato

    # Forca Tangencial no pinhao
    W_t = (T_p_in * 1000) / (d_p / 2)
    # Forca Radial no pinhao
    W_r = W_t * np.tan(φ_rad)

    # Velocidade tangencial
    ω_p_in = n_p_in * (2 * np.pi / 60)
    V_t = (d_p / 1000) * ω_p_in / 2

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
    σ_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
    σ_b_adm = S_at / (K_T * K_R)

    FS_flexao = np.inf if σ_b == 0 else σ_b_adm / σ_b

    # Analise de Contato/Pitting
    mg = N_c / N_p

    # Fator de Geometria (Superfície)
    I_geo = (np.sin(φ_rad) * np.cos(φ_rad) / 2) * (mg / (mg + 1))

    pitting = (W_t * K_s * K_m) / (b_face * d_p * I_geo * K_v)

    σ_c = np.inf
    if pitting >= 0 and (b_face * d_p * I_geo * K_v) != 0:
        σ_c = C_p * np.sqrt(pitting)

    σ_c_adm = S_ac / C_R

    FS_pitting = 0.0
    if σ_c == np.inf:
        FS_pitting = 0.0
    elif σ_c == 0:
        FS_pitting = np.inf
    else:
        FS_pitting = (σ_c_adm / σ_c)**2

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



estagio1_iteracoes = []
estagio2_iteracoes = []

modulos_padronizados = [1.5, 2, 2.5, 3, 4, 5, 6, 8]
Np_1 = 18
b_Face_1 = 10
Np_2 = 18
b_Face_2 = 10

resultado_estagio_1 = None
resultado_estagio_2 = None

#
def escreve_estagio(file, nome_estagio, res):
    if res is None:
        file.write("Nenhum módulo padrão foi suficiente.\n")
        return

    file.write(f"\n--- {nome_estagio}: aprovado ---\n")
    file.write(f"  Modulo (m): {res['m']:.1f} mm\n")
    file.write(f"  Largura da Face (b): {res['b_face']:.1f} mm\n")
    file.write(f"  Angulo de Pressao: {20.0} graus\n")
    file.write(f"  Razao de Contato (mp): {res['mp']:.3f}\n")
    file.write("\n  Pinhao:\n")
    file.write(f"    N de Dentes (Np): {res['N_p']}\n")
    file.write(f"    Diametro Primitivo (dp): {res['d_p']:.2f} mm\n")
    file.write(f"    Diametro Externo (dep): {res['d_ext_p']:.2f} mm\n")
    file.write("  Coroa:\n")
    file.write(f"    N de Dentes (Nc): {res['N_c']}\n")
    file.write(f"    Diametro Primitivo (dc): {res['d_c']:.2f} mm\n")
    file.write(f"    Diametro Externo (dec): {res['d_ext_c']:.2f} mm\n")
    file.write("  Conjunto:\n")
    file.write(f"    Razao de Transmissao (i): {res['i_efetiva']:.3f}\n")
    file.write(f"    Distancia entre Centros (C): {res['C']:.2f} mm\n")
    file.write("  Fatores de Seguranca:\n")
    file.write(f"    FS (Flexao): {res['FS_flexao']:.2f}\n")
    file.write(f"    FS (Superfície): {res['FS_pitting']:.2f}\n")


# Estágio 1
for modulo in modulos_padronizados:
    resultado = dimensionar_estagio(T_in, n_in, i_alvo_1, modulo, Np_1, b_Face_1, η_1)

    if resultado['FS_flexao'] > 1.5 and resultado['FS_pitting'] > 1.5:
        resultado_estagio_1 = resultado
        break

# Estágio 2
if resultado_estagio_1 is not None:
    # Saída do estágio 1 é a entrada do estágio 2
    T_p2_in = resultado_estagio_1['T_p_out']
    n_p2_in = resultado_estagio_1['n_p_out']

    for modulo in modulos_padronizados:
        resultado = dimensionar_estagio(T_p2_in, n_p2_in, i_alvo_2, modulo, Np_2, b_Face_2, η_2)

        if resultado['FS_flexao'] > 1.5 and resultado['FS_pitting'] > 1.5:
            resultado_estagio_2 = resultado
            break

with open("resultados_engrenagem.txt", "w") as f:
    # Estágio 1
    escreve_estagio(f, "Estágio 1 (Entrada)", resultado_estagio_1)

    # Estágio 2
    escreve_estagio(f, "Estágio 2 (Saída)", resultado_estagio_2)

    if resultado_estagio_1 and resultado_estagio_2:
        i_total_efetiva_final = resultado_estagio_1['i_efetiva'] * resultado_estagio_2['i_efetiva']
        erro_final = (i_total_efetiva_final - i_total) / i_total

        f.write(f"Razao de Transmissao Total Efetiva: {i_total_efetiva_final:.3f}\n")
        f.write(f"Erro Percentual: {erro_final * 100:.2f}%\n")

        if abs(erro_final) <= 0.05:
            f.write("Dentro do limite de+/- 5%. Projeto aprovado.\n")
        else:
            f.write("Fora do limite de +/- 5%. Projeto com falha.\n")
    else:
        f.write("Dimensionamento Falhou\n")


    # ---  Parametros para outros setores ---
    if resultado_estagio_1 and resultado_estagio_2:
        # Eixo 1 Entrada
        f.write(f"  Rotacao (n_eixo1): {n_in:.1f} rpm\n")
        f.write(f"  Torque (T_eixo1): {T_in:.2f} N.m\n")
        f.write("  Forcas no Pinhao 1 (aplicadas no Eixo 1):\n")
        f.write(f"    Forca Tangencial (W_t1): {resultado_estagio_1['W_t']:.2f} N\n")
        f.write(f"    Forca Radial (W_r1): {resultado_estagio_1['W_r']:.2f} N\n")

        # Eixo 2
        f.write(f"  Rotacao (n_eixo2): {resultado_estagio_1['n_p_out']:.1f} rpm\n")
        f.write(f"  Torque (T_eixo2): {resultado_estagio_1['T_p_out']:.2f} N.m\n")
        f.write("  Forças da Coroa 1 (aplicadas no Eixo 2):\n")
        f.write(f"    Força Tangencial (W_t_c1): {resultado_estagio_1['W_t']:.2f} N\n")
        f.write(f"    Força Radial (W_r_c1): {resultado_estagio_1['W_r']:.2f} N\n")
        f.write("  Forças no Pinhão 2 (aplicadas no Eixo 2):\n")
        f.write(f"    Força Tangencial (W_t_p2): {resultado_estagio_2['W_t']:.2f} N\n")
        f.write(f"    Força Radial (W_r_p2): {resultado_estagio_2['W_r']:.2f} N\n")

        # Eixo 3 Saida
        f.write(f"  Rotacao (n_eixo3): {resultado_estagio_2['n_p_out']:.1f} rpm\n")
        f.write(f"  Torque (T_eixo3): {resultado_estagio_2['T_p_out']:.2f} N.m\n")
        f.write("  Forças da Coroa 2 (aplicadas no Eixo 3):\n")
        f.write(f"    Força Tangencial (W_t_c2): {resultado_estagio_2['W_t']:.2f} N\n")
        f.write(f"    Força Radial (W_r_c2): {resultado_estagio_2['W_r']:.2f} N\n")

    else:
        f.write("\nSem parametros adicionais para outros setores devido a falha no dimensionamento.\n")
