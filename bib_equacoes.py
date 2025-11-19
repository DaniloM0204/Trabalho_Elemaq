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

def calcula_c_superf(S_ut):
    """
    Analisa TODOS os acabamentos e retorna o MAIOR fator (Melhor cenário).
    """
    A_b = [(1.58, -0.085),   # Retificado
           (4.51, -0.265),   # Usinado
           (57.7, -0.718),   # Laminado
           (272, -0.995)]    # Forjado
    
    nomes = ['retificado', 'usinado', 'laminado', 'forjado']
    c_lista = []
    
    for i in range(len(A_b)):
        A = A_b[i][0]
        b = A_b[i][1]
        
        # Cálculo da fórmula: C_superf = A * (Sut)^b
        fator = A * (S_ut ** b)
        
        # O fator não deve ser maior que 1.0
        if fator > 1.0: fator = 1.0
            
        c_lista.append(fator)

    # Encontrar o maior valor
    maior_valor = max(c_lista)
    indice_maior = c_lista.index(maior_valor)
    melhor_acabamento = nomes[indice_maior]
        
    return maior_valor, melhor_acabamento

def calcula_fator_tamanho(d):
    """Calcula C_tamanho segundo Norton para Flexão/Torção"""
    if d <= 8:
        return 1.0
    elif d <= 250:
        return 1.189 * (d ** -0.097)
    else:
        return 0.6

def calcula_diametro_goodman(Nf, Kf, Ma, Kfs, Ta, Sf, Kfm, Mm, Kfsm, Tm, Sut):
    """Equação de Goodman Modificado (Eq 43 do relatório)"""
    
    # Tensão Alternada
    termo_alternado_num = np.sqrt((Kf * Ma)**2 + 0.75 * (Kfs * Ta)**2)
    termo_alternado = termo_alternado_num / Sf
    
    # Tensão Média
    termo_medio_num = np.sqrt((Kfm * Mm)**2 + 0.75 * (Kfsm * Tm)**2)
    termo_medio = termo_medio_num / Sut
    
    multiplicador = (32 * Nf) / np.pi
    
    d_min = (multiplicador * (termo_alternado + termo_medio)) ** (1/3)
    return d_min

def dimensionar_eixo_completo(Ma, Tm, Sut, Se_linha, Nf=1.5, tolerancia=0.05):
    """
    Realiza a iteração do diâmetro considerando o erro aceitável de 5%.
    Automaticamente seleciona o melhor acabamento superficial.
    """
    
    # 1. Fatores constantes
    C_carreg = 1.0 # Flexão rotativa padrão
    C_temp = 1.0
    C_conf = 0.814 # 99% confiabilidade
    
    # Fator de superfície (Otimizado)
    C_sup, acabamento = calcula_c_superf(S_ut=Sut)
    
    # 2. Loop de Iteração
    d_atual = 30.0 # Chute inicial em mm
    erro = 1.0
    iteracao = 0
    
    print(f"--- Iniciando Iteração do Eixo (Acabamento: {acabamento}) ---")
    
    while erro > tolerancia and iteracao < 20:
        # Calcula fator de tamanho com o diâmetro atual
        C_tam = calcula_fator_tamanho(d_atual)
        
        # Atualiza Limite de Fadiga (Sf)
        Sf_corrigido = Se_linha * C_carreg * C_tam * C_sup * C_temp * C_conf
        
        # Calcula novo diâmetro via Goodman
        # Assumindo Kf e Kfs genéricos (ex: chaveta) se não passados
        d_novo = calcula_diametro_goodman(
            Nf=Nf, Kf=2.0, Ma=Ma, Kfs=1.6, Ta=0, 
            Sf=Sf_corrigido, Kfm=1.0, Mm=0, Kfsm=1.0, Tm=Tm, Sut=Sut
        )
        
        # Cálculo do erro relativo
        erro = abs(d_novo - d_atual) / d_atual
        
        print(f"Iter {iteracao+1}: d_calc={d_novo:.2f}mm (C_tam={C_tam:.3f}, Sf={Sf_corrigido:.1f}) - Erro: {erro*100:.2f}%")
        
        d_atual = d_novo
        iteracao += 1
        
    return d_atual, Sf_corrigido
