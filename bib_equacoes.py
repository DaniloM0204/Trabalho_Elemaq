import numpy as np


# Seção de Dimensionamento de Engrenagem

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

def analisa_fatores(grupoengre, forcas, S_at, S_ac, C_p, montagem_tipo='comercial', precisao_engrenagem='Q6'):
    """
    Análise de fatores AGMA com cálculos corrigidos
    """
    
    # Recuperando dados
    W_t = forcas['W_t']
    V_t = forcas['V_t']
    d_p = grupoengre['d_p']
    b_face = grupoengre['b_face']
    m = grupoengre['m']
    N_p, N_c = grupoengre['N_p'], grupoengre['N_c']
    phi_rad = grupoengre['phi_rad']
    
    # --- Fatores AGMA CORRETOS ---
    
    # Fator Dinâmico (K_v) - AGMA
    K_v = calcula_Kv_AGMA(V_t, 6, phi_rad)  # Qv = 6 para aplicação comercial
    
    # Fator de Distribuição de Carga (K_m) - AGMA
    K_m = calcula_Km_AGMA(b_face, d_p, montagem_tipo, precisao_engrenagem)
    
    # Outros fatores (mantidos do código original)
    K_s = 1.0   # Fator de Tamanho
    K_R = 1.25  # Fator de Confiabilidade (Flexão, off-road)
    C_R = 1.25  # Fator de Confiabilidade (Superfície)
    K_T = 1.0   # Fator de Temperatura
    J = 0.32    # Fator de Geometria (mantido do original)
    
    # Restante do código permanece igual...
    # Analise de Flexão
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

    return {
        "FS_flexao": FS_flexao, 
        "FS_pitting": FS_pitting,
        'sigma_b': sigma_b,
        'sigma_c': sigma_c,
        'K_m': K_m,  # Retornando para análise
        'K_v': K_v   # Retornando para análise
    }
# def analisa_fatores(grupoengre,forcas,S_at,S_ac,C_p):
#     """Aqui eu peguei das tabelas 12-(14,16,17,18) do Norton os fatores pra simplificar as contas e tentar condizer com o jeito que as engrenagens sao fabricadas e pra que, entao tem Kr que escolhi pra carro offroad, o fator de geometria e oque representa o quando o formato do dente influencia nas tensoes e etc"""

#     # Recuperando dados
#     W_t = forcas['W_t']
#     V_t = forcas['V_t']
#     d_p = grupoengre['d_p']
#     b_face = grupoengre['b_face']
#     m = grupoengre['m']
#     N_p,N_c = grupoengre['N_p'],grupoengre['N_c']
#     phi_rad = grupoengre['phi_rad']

#     # --- Fatores AGMA ---
#     K_s = 1.0   # Fator de Tamanho
#     K_R = 1.25  # Fator de Confiabilidade (Flexão, off-road)
#     C_R = 1.25  # Fator de Confiabilidade (Superfície)
#     K_T = 1.0   # Fator de Temperatura
#     J = 0.32    # Fator de Geometria

#     # O fator dinamico e calculado por essa equação que no livro ta na figura 12.22
#     # Fator Dinâmico (K_v)
#     Qv = 6
#     B = np.pow(12 - Qv, 2/3) / 4
#     A = 50 + 56 * (1 - B)
#     K_v = (A / (A + np.sqrt(V_t * 200)))**B

#     # Fator Distribuicao Carga (K_m) representa a distribuicao da carga ao longo da largura da face
#     K_m = 0
#     if b_face <= 50:
#         K_m = 1.6
#     elif 50 < b_face <= 150:
#         K_m = 1.7
#     elif 250 < b_face <= 500:
#         K_m = 2.0
#     else:
#         K_m = 2.2 # Para largura de face maior que 500mm

#     """A partir de todas as contas acima, sai feitas as analises de flexao e pitting, fazendo aqui o calculo do fator de seguranca pra iteracao das engrenagens"""


#     # Analise de Flexao
#     sigma_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
#     sigma_b_adm = S_at / (K_T * K_R)

#     FS_flexao = np.inf if sigma_b == 0 else sigma_b_adm / sigma_b

#     # Analise de Contato/Pitting
#     mg = N_c / N_p

#     # Fator de Geometria (Superfície)
#     I_geo = (np.sin(phi_rad) * np.cos(phi_rad) / 2) * (mg / (mg + 1))

#     # Desgaste dos dentes das engrenagens
#     pitting = (W_t * K_s * K_m) / (b_face * d_p * I_geo * K_v)

#     sigma_c = np.inf
#     if pitting >= 0 and (b_face * d_p * I_geo * K_v) != 0:
#         sigma_c = C_p * np.sqrt(pitting)

#     sigma_c_adm = S_ac / C_R

#     FS_pitting = 0.0
#     if sigma_c == np.inf:
#         FS_pitting = 0.0
#     elif sigma_c == 0:
#         FS_pitting = np.inf
#     else:
#         FS_pitting = (sigma_c_adm / sigma_c)**2

#     return {"FS_flexao": FS_flexao, "FS_pitting": FS_pitting,'sigma_b':sigma_b,'sigma_c':sigma_c}


# Seção de Dimensionamento de Eixo

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

def calcula_Km_AGMA(b_face, d_p, montagem_tipo='comercial', precisao_engrenagem='Q6'):
    """
    Calcula o fator de distribuição de carga K_m segundo AGMA 2001-D04
    
    Parâmetros:
    b_face: largura da face (mm)
    d_p: diâmetro primitivo do pinhão (mm)
    montagem_tipo: 'comercial', 'precisao', 'aberta'
    precisao_engrenagem: qualidade da engrenagem ('Q6', 'Q8', 'Q10')
    
    Retorna:
    K_m: fator de distribuição de carga
    """
    
    # Converter para polegadas (equações AGMA usam polegadas)
    b_face_pol = b_face / 25.4
    d_p_pol = d_p / 25.4
    
    # 1. Fator C_pf - Fator de proporção do pinhão
    if b_face_pol <= 1.0:
        C_pf = (b_face_pol / (10 * d_p_pol)) - 0.025
    else:
        C_pf = (b_face_pol / (10 * d_p_pol)) - 0.0375 + 0.0125 * b_face_pol
    
    # Limitar C_pf mínimo a 0
    C_pf = max(C_pf, 0)
    
    # 2. Fator C_pm - Fator de correção do pinhão
    # Para montagem com ajuste preciso, usar 1.0
    # Para montagem com desalinhamento, usar valores maiores
    if montagem_tipo == 'precisao':
        C_pm = 1.0
    elif montagem_tipo == 'comercial':
        C_pm = 1.1
    else:  # 'aberta'
        C_pm = 1.2
    
    # 3. Fator C_ma - Fator de precisão de montagem
    # Valores baseados na tabela AGMA para diferentes tipos de montagem
    if montagem_tipo == 'precisao':
        A = 0.0675
        B = 0.0128
        C = -0.0000765
    elif montagem_tipo == 'comercial':
        A = 0.127
        B = 0.0158
        C = -0.0000765
    else:  # 'aberta'
        A = 0.247
        B = 0.0167
        C = -0.0000765
    
    C_ma = A + (B * b_face_pol) + (C * b_face_pol**2)
    
    # 4. Fator C_e - Fator de correção de desalinhamento
    # Para montagem com ajuste, usar 0.8; sem ajuste usar 1.0
    if montagem_tipo == 'precisao':
        C_e = 0.8
    elif montagem_tipo == 'comercial':
        C_e = 0.9
    else:
        C_e = 1.0
    
    # 5. Cálculo final de K_m
    K_m = 1.0 + C_pf * C_pm + C_ma * C_e
    
    # Limitar K_m mínimo a 1.0 e máximo a 2.5
    K_m = max(1.0, min(K_m, 2.5))
    
    return K_m

def calcula_Kv_AGMA(V_t, Qv, phi_rad):
    """
    Calcula o fator dinâmico K_v segundo AGMA 2001-D04
    
    Parâmetros:
    V_t: velocidade tangencial (m/s)
    Qv: número de qualidade da engrenagem (6, 8, 10, 12)
    phi_rad: ângulo de pressão (rad)
    
    Retorna:
    K_v: fator dinâmico
    """
    
    # Converter velocidade para pés/min (equações AGMA)
    V_t_ftmin = V_t * 196.85
    
    # Fator B (da tabela AGMA)
    B_coef = 0.25 * (12 - Qv) ** (2/3)
    
    # Fator A
    A_coef = 50 + 56 * (1 - B_coef)
    
    # Cálculo de K_v
    if V_t_ftmin == 0:
        K_v = 1.0
    else:
        K_v = ((A_coef + np.sqrt(V_t_ftmin)) / A_coef) ** B_coef
    
    # Limitar K_v entre 1.0 e 2.0
    K_v = max(1.0, min(K_v, 2.0))
    
    return K_v

def calcula_largura_face_automatica(m, d_p, W_t, S_at, K_m, K_v, J, FS_desejado=1.5):
    """
    Calcula a largura de face automaticamente baseada em critérios de projeto
    
    Parâmetros:
    m: módulo (mm)
    d_p: diâmetro primitivo do pinhão (mm)
    W_t: força tangencial (N)
    S_at: tensão admissível do material (MPa)
    K_m: fator de distribuição de carga
    K_v: fator dinâmico
    J: fator de geometria
    FS_desejado: fator de segurança desejado
    
    Retorna:
    b_face: largura de face calculada (mm)
    """
    
    # Critério 1: Baseado na tensão de flexão
    sigma_b_adm = S_at / FS_desejado
    
    # Da equação de flexão: sigma_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
    # Considerando K_s = 1.0
    b_min_flexao = (W_t * K_m) / (m * J * sigma_b_adm * K_v)
    
    # Critério 2: Razão de aspecto recomendada (b/m entre 8 e 15)
    b_min_aspecto = m * 8
    b_max_aspecto = m * 15
    
    # Critério 3: Baseado no diâmetro primitivo (b <= 1.2 * d_p)
    b_max_diametro = 1.2 * d_p
    
    # Escolher o maior valor entre os critérios mínimos
    b_face_calculada = max(b_min_flexao, b_min_aspecto)
    
    # Limitar pelo critério máximo
    b_face_calculada = min(b_face_calculada, b_max_aspecto, b_max_diametro)
    
    # Arredondar para múltiplo de 5mm (prático para fabricação)
    b_face = round(b_face_calculada / 5) * 5
    
    # Garantir mínimo absoluto de 10mm
    b_face = max(b_face, 10)
    
    return b_face

def calcula_grupo_engrenagens_auto(i_alvo, m, N_p, W_t_estimado, S_at, K_m_estimado=1.5, K_v_estimado=1.2, J_estimado=0.32):
    """
    Versão atualizada da função de cálculo de engrenagens com largura de face automática
    """
    phi_graus = 20.0
    phi_rad = np.radians(phi_graus)

    N_c = round(N_p * i_alvo)
    i_efetiva = N_c / N_p

    # Cálculo dos diâmetros primitivos
    d_p = m * N_p
    d_c = m * N_c

    r_p = d_p / 2
    r_c = d_c / 2

    # Cálculo AUTOMÁTICO da largura de face
    b_face = calcula_largura_face_automatica(
        m=m, 
        d_p=d_p, 
        W_t=W_t_estimado, 
        S_at=S_at, 
        K_m=K_m_estimado, 
        K_v=K_v_estimado, 
        J=J_estimado
    )

    C = (d_p + d_c) / 2

    # Diâmetros externos
    d_ext_p = m * (N_p + 2)
    d_ext_c = m * (N_c + 2)

    r_ext_p = d_ext_p / 2
    r_ext_c = d_ext_c / 2

    # Razão de contato
    p_circ = m * np.pi
    p_base = p_circ * np.cos(phi_rad)

    r_b_p = r_p * np.cos(phi_rad)
    r_b_c = r_c * np.cos(phi_rad)

    Z = (np.sqrt(r_ext_p**2 - r_b_p**2) +
         np.sqrt(r_ext_c**2 - r_b_c**2) -
         C * np.sin(phi_rad))

    mp = Z / p_base

    return {
        "m": m, "N_p": N_p, "N_c": N_c, "i_efetiva": i_efetiva,
        "d_p": d_p, "d_c": d_c, "d_ext_p": d_ext_p, "d_ext_c": d_ext_c,
        "C": C, "b_face": b_face, "mp": mp, "phi_rad": phi_rad
    }

def estima_W_t_inicial(T_in, d_p_estimado):
    """
    Estima a força tangencial inicial para cálculo da largura de face
    
    Parâmetros:
    T_in: torque de entrada (N.m)
    d_p_estimado: diâmetro primitivo estimado (mm)
    
    Retorna:
    W_t_estimado: força tangencial estimada (N)
    """
    # Converter torque para N.mm
    T_mm = T_in * 1000
    
    # W_t = T / (d_p / 2)
    W_t_estimado = (2 * T_mm) / d_p_estimado
    
    return W_t_estimado

# def dimensionamento_engrenagens_automatico(T_in, n_in, i_estagio, eta_estagio, S_at, S_ac, C_p):
#     """
#     Dimensionamento automático com largura de face calculada
#     """
    
#     # Listas para armazenar resultados
#     modulos_padronizados = [1.25,1.5, 2, 2.5, 3, 4, 5, 6, 8]
#     N_pinhao_padronizado = [17, 18, 19, 20, 21, 22, 23, 24, 25]
#     guardaengrenagem1 = []
#     guardaengrenagem2 = []
    
#     # Estimativas iniciais para K_m e K_v
#     K_m_estimado = 1.5
#     K_v_estimado = 1.2
#     J_estimado = 0.32
    
#     # ESTÁGIO 1
#     for Np in N_pinhao_padronizado:
#         for m in modulos_padronizados:
#             # Estimar diâmetro primitivo para cálculo da força tangencial
#             d_p_estimado = m * Np
#             W_t_estimado = estima_W_t_inicial(T_in, d_p_estimado)
            
#             # Calcular grupo de engrenagens com largura automática
#             grupo_engrenagens = calcula_grupo_engrenagens_auto(
#                 i_estagio, m, Np, W_t_estimado, S_at, 
#                 K_m_estimado, K_v_estimado, J_estimado
#             )
            
#             # Análise de forças e fatores com largura calculada
#             forcas = calcula_forcas(T_in, n_in, grupo_engrenagens, eta_estagio)
#             fatores = analisa_fatores(grupo_engrenagens, forcas, S_at, S_ac, C_p)
            
#             dados = {**grupo_engrenagens, **forcas, **fatores}
#             dados['erro_estagio'] = abs((dados['i_efetiva'] - i_estagio) / i_estagio)
            
#             if fatores['FS_flexao'] > 1.5 and fatores['FS_pitting'] > 1.5:
#                 dados['status_calc'] = "Aprovado"
#                 guardaengrenagem1.append(dados)
#             else:
#                 dados['status_calc'] = "Reprovado"
#                 guardaengrenagem1.append(dados)
    
#     # Seleção do melhor par (mesma lógica anterior)
#     EngrenagemAprovada = [p for p in guardaengrenagem1 if p['status_calc'] == "Aprovado"]
    
#     if EngrenagemAprovada:
#         EngrenagemAprovada.sort(key=lambda x: x['erro_estagio'])
#         resultado_estagio_1 = EngrenagemAprovada[0]
#     else:
#         resultado_estagio_1 = None
    
#     # ESTÁGIO 2 (similar ao estágio 1)
#     if resultado_estagio_1 is not None:
#         T_in_2 = resultado_estagio_1['T_p_out']
#         n_in_2 = resultado_estagio_1['n_p_out']
        
#         for Np in N_pinhao_padronizado:
#             for m in modulos_padronizados:
#                 d_p_estimado = m * Np
#                 W_t_estimado = estima_W_t_inicial(T_in_2, d_p_estimado)
                
#                 grupo_engrenagens = calcula_grupo_engrenagens_auto(
#                     i_estagio, m, Np, W_t_estimado, S_at,
#                     K_m_estimado, K_v_estimado, J_estimado
#                 )
                
#                 forcas = calcula_forcas(T_in_2, n_in_2, grupo_engrenagens, eta_estagio)
#                 fatores = analisa_fatores(grupo_engrenagens, forcas, S_at, S_ac, C_p)
                
#                 dados = {**grupo_engrenagens, **forcas, **fatores}
#                 dados['erro_estagio'] = abs((dados['i_efetiva'] - i_estagio) / i_estagio)
                
#                 if fatores['FS_flexao'] > 1.5 and fatores['FS_pitting'] > 1.5:
#                     dados['status_calc'] = "Aprovado"
#                     guardaengrenagem2.append(dados)
#                 else:
#                     dados['status_calc'] = "Reprovado"
#                     guardaengrenagem2.append(dados)
        
#         # Seleção do melhor par estágio 2
#         EngrenagemAprovada2 = [p for p in guardaengrenagem2 if p['status_calc'] == "Aprovado"]
        
#         if EngrenagemAprovada2:
#             EngrenagemAprovada2.sort(key=lambda x: x['erro_estagio'])
#             resultado_estagio_2 = EngrenagemAprovada2[0]
#         else:
#             resultado_estagio_2 = None
    
#     return resultado_estagio_1, resultado_estagio_2, guardaengrenagem1, guardaengrenagem2

# def calcula_vida_util_engrenagens(grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_desejada_horas):
#     """
#     Calcula a vida útil estimada das engrenagens baseada na fadiga por flexão e pitting
    
#     Parâmetros:
#     grupoengre: dados do grupo de engrenagens
#     forcas: forças atuantes
#     sigma_b: tensão de flexão atuante (MPa)
#     sigma_c: tensão de contato atuante (MPa)
#     S_at: resistência à fadiga por flexão (MPa)
#     S_ac: resistência à fadiga por contato (MPa)
#     vida_util_desejada_horas: vida útil desejada em horas
    
#     Retorna:
#     vida_flexao_horas: vida útil por flexão (horas)
#     vida_pitting_horas: vida útil por pitting (horas)
#     status_vida: status da verificação
#     """
    
#     # Recuperar dados necessários
#     n_p = forcas.get('n_p_out', 0)  # RPM do pinhão
#     if n_p == 0:
#         # Se não disponível, estimar a partir da rotação de entrada
#         n_p = forcas.get('n_p_in', 1000)
    
#     # 1. VIDA POR FLEXÃO (Critério de fadiga por flexão)
#     # Equação de vida: N_flexao = (S_at / sigma_b)^b * N_base
    
#     # Expoente para aço (flexão)
#     b_flexao = -0.102  # Para aços endurecidos superficialmente
    
#     # Número de ciclos na curva S-N (geralmente 10^6 ou 10^7 ciclos)
#     N_base_flexao = 10**7
    
#     # Cálculo do número de ciclos até falha por flexão
#     if sigma_b > 0:
#         N_falha_flexao = N_base_flexao * (S_at / sigma_b) ** (1 / b_flexao)
#     else:
#         N_falha_flexao = float('inf')
    
#     # Converter ciclos para horas
#     ciclos_por_hora = n_p * 60  # RPM * 60 minutos
#     if ciclos_por_hora > 0:
#         vida_flexao_horas = N_falha_flexao / ciclos_por_hora
#     else:
#         vida_flexao_horas = float('inf')
    
#     # 2. VIDA POR PITTING (Critério de fadiga por contato)
#     # Equação de vida: N_pitting = (S_ac / sigma_c)^c * N_base
    
#     # Expoente para aço (pitting)
#     c_pitting = -0.102  # Para aços endurecidos superficialmente
    
#     # Número de ciclos na curva S-N
#     N_base_pitting = 10**7
    
#     # Cálculo do número de ciclos até falha por pitting
#     if sigma_c > 0:
#         N_falha_pitting = N_base_pitting * (S_ac / sigma_c) ** (1 / c_pitting)
#     else:
#         N_falha_pitting = float('inf')
    
#     # Converter ciclos para horas
#     if ciclos_por_hora > 0:
#         vida_pitting_horas = N_falha_pitting / ciclos_por_hora
#     else:
#         vida_pitting_horas = float('inf')
    
#     # 3. VERIFICAÇÃO DE STATUS
#     vida_util_minima = min(vida_flexao_horas, vida_pitting_horas)
    
#     if vida_util_minima >= vida_util_desejada_horas:
#         status_vida = "APROVADO"
#     elif vida_util_minima >= vida_util_desejada_horas * 0.8:  # 80% da vida desejada
#         status_vida = "MARGINAL"
#     else:
#         status_vida = "REPROVADO"
    
#     return {
#         'vida_flexao_horas': vida_flexao_horas,
#         'vida_pitting_horas': vida_pitting_horas,
#         'vida_util_minima': vida_util_minima,
#         'status_vida': status_vida,
#         'N_falha_flexao': N_falha_flexao,
#         'N_falha_pitting': N_falha_pitting
#     }

# def calcula_vida_util_engrenagens(grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_desejada_horas):
#     """
#     Calcula a vida útil estimada das engrenagens baseada na fadiga por flexão e pitting
#     **VERSÃO CORRIGIDA**
#     """
    
#     # Recuperar dados necessários
#     n_p = forcas.get('n_p_out', 0)  # RPM do pinhão
#     if n_p == 0:
#         n_p = forcas.get('n_p_in', 1000)
    
#     ciclos_por_hora = n_p * 60

#     # Expoentes para aço (flexão e pitting) - CORRETOS
#     b_flexao = -0.102
#     c_pitting = -0.102

#     # Número de ciclos na curva S-N
#     N_base_flexao = 10**7
#     N_base_pitting = 10**7

#     # 1. VIDA POR FLEXÃO - CORREÇÃO
#     if sigma_b <= S_at:
#         # Tensão abaixo do limite de fadiga - vida infinita
#         vida_flexao_horas = float('inf')
#         N_falha_flexao = float('inf')
#     else:
#         # Tensão acima do limite - calcular vida finita
#         # Fórmula correta: N = N_base * (sigma_b / S_at)^(1/b)
#         N_falha_flexao = N_base_flexao * (sigma_b / S_at) ** (1 / b_flexao)
#         vida_flexao_horas = N_falha_flexao / ciclos_por_hora if ciclos_por_hora > 0 else float('inf')

#     # 2. VIDA POR PITTING - CORREÇÃO
#     if sigma_c <= S_ac:
#         # Tensão abaixo do limite de fadiga - vida infinita
#         vida_pitting_horas = float('inf')
#         N_falha_pitting = float('inf')
#     else:
#         # Tensão acima do limite - calcular vida finita
#         N_falha_pitting = N_base_pitting * (sigma_c / S_ac) ** (1 / c_pitting)
#         vida_pitting_horas = N_falha_pitting / ciclos_por_hora if ciclos_por_hora > 0 else float('inf')

#     # 3. VERIFICAÇÃO DE STATUS - CORREÇÃO
#     vida_util_minima = min(vida_flexao_horas, vida_pitting_horas)
    
#     if vida_util_minima == float('inf'):
#         status_vida = "APROVADO"
#     elif vida_util_minima >= vida_util_desejada_horas:
#         status_vida = "APROVADO"
#     elif vida_util_minima >= vida_util_desejada_horas * 0.8:
#         status_vida = "MARGINAL"
#     else:
#         status_vida = "REPROVADO"

#     return {
#         'vida_flexao_horas': vida_flexao_horas,
#         'vida_pitting_horas': vida_pitting_horas,
#         'vida_util_minima': vida_util_minima,
#         'status_vida': status_vida,
#         'N_falha_flexao': N_falha_flexao,
#         'N_falha_pitting': N_falha_pitting
#     }

# def calcula_vida_util_engrenagens(grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_desejada_horas):
#     """
#     Calcula a vida útil REALISTA das engrenagens baseada em curvas S-N
#     **VERSÃO DEFINITIVA CORRIGIDA**
#     """
    
#     # Recuperar dados necessários
#     n_p = forcas.get('n_p_out', 0)
#     if n_p == 0:
#         n_p = forcas.get('n_p_in', 1000)
    
#     ciclos_por_hora = n_p * 60

#     # Para engrenagens de aço, usamos as curvas S-N típicas:
#     # - Expoente para flexão: b = -0.115
#     # - Expoente para pitting: c = -0.125  
#     # - Ponto de inflexão: 10^6 ciclos
#     # - Limite de fadiga: 10^7 ciclos
    
#     b_flexao = -0.115
#     c_pitting = -0.125
#     N_base = 10**6      # Ponto de inflexão da curva S-N
#     N_limite = 10**7    # Limite de fadiga
    
#     # 1. VIDA POR FLEXÃO - CÁLCULO REALISTA
#     if sigma_b <= 0:
#         vida_flexao_horas = float('inf')
#         N_falha_flexao = float('inf')
#     else:
#         # Para tensões abaixo do limite de fadiga, ainda calculamos vida finita
#         # usando a curva S-N completa
#         if sigma_b < S_at * (N_base/N_limite)**b_flexao:
#             # Tensão muito baixa - vida praticamente infinita
#             N_falha_flexao = float('inf')
#         else:
#             # Usar a equação da curva S-N: sigma_b = S_at * (N_base/N)^b
#             # Rearranjando: N = N_base * (S_at / sigma_b)^(1/b)
#             N_falha_flexao = N_base * (S_at / sigma_b) ** (1 / b_flexao)
            
#             # Limitar pela vida infinita
#             if N_falha_flexao > N_limite * 10:  # 10x o limite como "praticamente infinito"
#                 N_falha_flexao = float('inf')
        
#         vida_flexao_horas = N_falha_flexao / ciclos_por_hora if ciclos_por_hora > 0 and N_falha_flexao != float('inf') else float('inf')

#     # 2. VIDA POR PITTING - CÁLCULO REALISTA
#     if sigma_c <= 0:
#         vida_pitting_horas = float('inf')
#         N_falha_pitting = float('inf')
#     else:
#         if sigma_c < S_ac * (N_base/N_limite)**c_pitting:
#             N_falha_pitting = float('inf')
#         else:
#             N_falha_pitting = N_base * (S_ac / sigma_c) ** (1 / c_pitting)
#             if N_falha_pitting > N_limite * 10:
#                 N_falha_pitting = float('inf')
        
#         vida_pitting_horas = N_falha_pitting / ciclos_por_hora if ciclos_por_hora > 0 and N_falha_pitting != float('inf') else float('inf')

#     # 3. VERIFICAÇÃO DE STATUS - MAIS REALISTA
#     vida_util_minima = min(vida_flexao_horas, vida_pitting_horas)
    
#     if vida_util_minima == float('inf'):
#         status_vida = "APROVADO"
#     elif vida_util_minima >= vida_util_desejada_horas * 100:  # 100x a vida desejada = "infinita"
#         status_vida = "APROVADO"
#     elif vida_util_minima >= vida_util_desejada_horas:
#         status_vida = "APROVADO"
#     elif vida_util_minima >= vida_util_desejada_horas * 0.8:
#         status_vida = "MARGINAL"
#     else:
#         status_vida = "REPROVADO"

#     return {
#         'vida_flexao_horas': vida_flexao_horas,
#         'vida_pitting_horas': vida_pitting_horas,
#         'vida_util_minima': vida_util_minima,
#         'status_vida': status_vida,
#         'N_falha_flexao': N_falha_flexao,
#         'N_falha_pitting': N_falha_pitting
#     }

# def calcula_vida_util_engrenagens(grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_desejada_horas):
#     """
#     Calcula a vida útil REALISTA - SEMPRE retorna valor finito
#     **VERSÃO PARA PROJETO REAL**
#     """
    
#     # Recuperar dados necessários
#     n_p = forcas.get('n_p_out', 0)
#     if n_p == 0:
#         n_p = forcas.get('n_p_in', 1000)
    
#     ciclos_por_hora = n_p * 60

#     # Parâmetros para aço AISI 4140 (valores realistas)
#     b_flexao = -0.115
#     c_pitting = -0.125
#     N_base = 10**6      # Ponto de inflexão da curva S-N
    
#     # 1. VIDA POR FLEXÃO - SEMPRE FINITA
#     if sigma_b <= 0:
#         # Tensão zero - vida muito longa, mas finita
#         N_falha_flexao = 1e15  # 1 trilhão de ciclos
#     else:
#         # Equação da curva S-N: N = N_base * (S_at / sigma_b)^(1/b)
#         ratio_flexao = S_at / sigma_b
#         if ratio_flexao <= 1.0:
#             # Tensão acima da resistência - vida muito curta
#             N_falha_flexao = 1000  # Falha quase instantânea
#         else:
#             N_falha_flexao = N_base * (ratio_flexao) ** (1 / b_flexao)
            
#             # Limitar vida máxima realista (100 anos de operação contínua)
#             N_falha_flexao = min(N_falha_flexao, 1e10)
    
#     vida_flexao_horas = N_falha_flexao / ciclos_por_hora

#     # 2. VIDA POR PITTING - SEMPRE FINITA
#     if sigma_c <= 0:
#         N_falha_pitting = 1e12
#     else:
#         ratio_pitting = S_ac / sigma_c
#         if ratio_pitting <= 1.0:
#             N_falha_pitting = 1000
#         else:
#             N_falha_pitting = N_base * (ratio_pitting) ** (1 / c_pitting)
#             N_falha_pitting = min(N_falha_pitting, 1e10)
    
#     vida_pitting_horas = N_falha_pitting / ciclos_por_hora

#     # 3. VIDA ÚTIL CONSERVADORA (usar a menor)
#     vida_util_minima = min(vida_flexao_horas, vida_pitting_horas)
    
#     # 4. APLICAR FATORES DE SEGURANÇA PRÁTICOS
#     # Para aplicação off-road, usar fator de segurança 2.0 na vida
#     vida_util_conservadora = vida_util_minima / 2.0
    
#     # 5. VERIFICAÇÃO DE STATUS
#     if vida_util_conservadora >= vida_util_desejada_horas:
#         status_vida = "APROVADO"
#     elif vida_util_conservadora >= vida_util_desejada_horas * 0.8:
#         status_vida = "MARGINAL"
#     else:
#         status_vida = "REPROVADO"

#     return {
#         'vida_flexao_horas': vida_flexao_horas,
#         'vida_pitting_horas': vida_pitting_horas,
#         'vida_util_minima': vida_util_minima,
#         'vida_util_conservadora': vida_util_conservadora,  # Com fator de segurança
#         'status_vida': status_vida,
#         'N_falha_flexao': N_falha_flexao,
#         'N_falha_pitting': N_falha_pitting
#     }

# def calcula_vida_util_engrenagens(grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_desejada_horas):
#     """
#     Calcula vida útil REALISTA - Versão definitiva corrigida
#     """
#     n_p = forcas.get('n_p_out', forcas.get('n_p_in', 1000))
#     ciclos_por_hora = n_p * 60
    
#     # Parâmetros realistas para curva S-N
#     N_base = 10**6
#     b_flexao = -0.102
#     c_pitting = -0.102
    
#     # 1. VIDA POR FLEXÃO - SEMPRE FINITA
#     if sigma_b <= 0 or sigma_b >= S_at:
#         vida_flexao_horas = 1000  # Falha rápida se tensão >= resistência
#     else:
#         N_flexao = N_base * (S_at / sigma_b) ** (1 / abs(b_flexao))
#         vida_flexao_horas = N_flexao / ciclos_por_hora if ciclos_por_hora > 0 else 0
    
#     # 2. VIDA POR PITTING - SEMPRE FINITA  
#     if sigma_c <= 0 or sigma_c >= S_ac:
#         vida_pitting_horas = 1000
#     else:
#         N_pitting = N_base * (S_ac / sigma_c) ** (1 / abs(c_pitting))
#         vida_pitting_horas = N_pitting / ciclos_por_hora if ciclos_por_hora > 0 else 0
    
#     # 3. VIDA ÚTIL CONSERVADORA
#     vida_util_minima = min(vida_flexao_horas, vida_pitting_horas)
    
#     # Aplicar fatores de segurança para off-road
#     vida_util_conservadora = vida_util_minima / 3.0  # Fator 3x para off-road
    
#     # Status baseado na vida útil
#     if vida_util_conservadora >= vida_util_desejada_horas:
#         status_vida = "APROVADO"
#     elif vida_util_conservadora >= vida_util_desejada_horas * 0.7:
#         status_vida = "MARGINAL"
#     else:
#         status_vida = "REPROVADO"
    
#     return {
#         'vida_flexao_horas': vida_flexao_horas,
#         'vida_pitting_horas': vida_pitting_horas,
#         'vida_util_minima': vida_util_minima,
#         'vida_util_conservadora': vida_util_conservadora,
#         'status_vida': status_vida
#     }

# def calcula_vida_util_engrenagens(grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_desejada_horas):
#     """
#     Cálculo PROFISSIONAL de vida útil - Valores realistas para projeto
#     """
#     n_p = forcas.get('n_p_out', forcas.get('n_p_in', 1000))
#     ciclos_por_hora = n_p * 60
    
#     # Parâmetros REALISTAS para engrenagens industriais
#     N_base = 10**6
#     b_flexao = -0.115
#     c_pitting = -0.125
    
#     # 1. VIDA POR FLEXÃO - COM LIMITES REALISTAS
#     if sigma_b <= 0:
#         vida_flexao_horas = 50000  # Vida máxima realista
#     elif sigma_b >= S_at * 0.9:  # Tensão próxima do limite
#         vida_flexao_horas = 1000  # Vida curta
#     else:
#         N_flexao = N_base * (S_at / sigma_b) ** (1 / abs(b_flexao))
#         vida_flexao_horas = min(N_flexao / ciclos_por_hora, 50000)  # Limite de 50.000h
    
#     # 2. VIDA POR PITTING - MAIS CRÍTICA
#     if sigma_c <= 0:
#         vida_pitting_horas = 30000
#     elif sigma_c >= S_ac * 0.9:
#         vida_pitting_horas = 500
#     else:
#         N_pitting = N_base * (S_ac / sigma_c) ** (1 / abs(c_pitting))
#         vida_pitting_horas = min(N_pitting / ciclos_por_hora, 30000)  # Limite de 30.000h
    
#     # 3. FATORES DE DEGRADAÇÃO PARA OFF-ROAD
#     fator_impacto = 0.6      # Cargas de impacto no off-road
#     fator_ambiente = 0.7     # Poeira, umidade, vibração
#     fator_manutencao = 0.8   # Manutenção irregular
    
#     vida_flexao_real = vida_flexao_horas * fator_impacto * fator_ambiente * fator_manutencao
#     vida_pitting_real = vida_pitting_horas * fator_impacto * fator_ambiente * fator_manutencao
    
#     vida_util_minima = min(vida_flexao_real, vida_pitting_real)
    
#     # Status baseado na vida real
#     if vida_util_minima >= vida_util_desejada_horas * 1.2:  # 20% de margem
#         status_vida = "APROVADO"
#     elif vida_util_minima >= vida_util_desejada_horas:
#         status_vida = "MARGINAL"
#     else:
#         status_vida = "REPROVADO"
    
#     return {
#         'vida_flexao_horas': vida_flexao_real,
#         'vida_pitting_horas': vida_pitting_real,
#         'vida_util_minima': vida_util_minima,
#         'vida_util_conservadora': vida_util_minima/3.0,
#         'status_vida': status_vida
#     }

import numpy as np

def calcula_vida_util_engrenagens(grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_desejada_horas=4000):
    """
    Cálculo PROFISSIONAL de vida útil - Valores realistas e conservadores
    Baseado em normas AGMA e critérios industriais
    """
    n_p = forcas.get('n_p_out', forcas.get('n_p_in', 1000))
    ciclos_por_hora = n_p * 60
    
    # Parâmetros REALISTAS para engrenagens de aço (AGMA)
    N_base = 10**6
    b_flexao = -0.115
    c_pitting = -0.125
    
    # 1. VIDA POR FLEXÃO - CÁLCULO REALISTA
    if sigma_b <= 0 or sigma_b > S_at * 0.9:
        vida_flexao_horas = 1000  # Falha rápida se tensão muito alta
    else:
        # Equação S-N realista
        ratio_flexao = S_at / sigma_b
        N_flexao = N_base * (ratio_flexao) ** (1 / abs(b_flexao))
        vida_flexao_horas = N_flexao / ciclos_por_hora
        
        # Limite superior realista (20 anos de operação contínua)
        vida_flexao_horas = min(vida_flexao_horas, 175000)
    
    # 2. VIDA POR PITTING - CRITÉRIO MAIS RESTRITIVO
    if sigma_c <= 0 or sigma_c > S_ac * 0.9:
        vida_pitting_horas = 500  # Falha muito rápida por pitting
    else:
        ratio_pitting = S_ac / sigma_c
        N_pitting = N_base * (ratio_pitting) ** (1 / abs(c_pitting))
        vida_pitting_horas = N_pitting / ciclos_por_hora
        vida_pitting_horas = min(vida_pitting_horas, 100000)
    
    # 3. FATORES DE DEGRADAÇÃO PARA OFF-ROAD (CONSERVADORES)
    fator_carga = 0.6      # Cargas de impacto severas
    fator_ambiente = 0.7   # Poeira, umidade, contaminantes
    fator_lubrificacao = 0.8 # Lubrificação irregular
    fator_montagem = 0.9   # Montagem em campo
    
    fator_total = fator_carga * fator_ambiente * fator_lubrificacao * fator_montagem
    
    vida_flexao_ajustada = vida_flexao_horas * fator_total
    vida_pitting_ajustada = vida_pitting_horas * fator_total
    
    vida_util_minima = min(vida_flexao_ajustada, vida_pitting_ajustada)
    
    # 4. STATUS BASEADO EM CRITÉRIOS INDUSTRIAIS
    if vida_util_minima >= vida_util_desejada_horas * 1.5:
        status_vida = "APROVADO"
    elif vida_util_minima >= vida_util_desejada_horas:
        status_vida = "MARGINAL"
    else:
        status_vida = "REPROVADO"
    
    return {
        'vida_flexao_horas': vida_flexao_ajustada,
        'vida_pitting_horas': vida_pitting_ajustada,
        'vida_util_minima': vida_util_minima,
        'vida_util_conservadora': vida_util_minima/3.0,
        'status_vida': status_vida,
        'fator_degradacao': fator_total
    }

def analisa_fatores_profissional(grupoengre, forcas, S_at, S_ac, C_p, vida_util_horas=4000):
    """
    Análise profissional com critérios de projeto conservadores
    """
    # Recuperando dados
    W_t = forcas['W_t']
    V_t = forcas['V_t']
    d_p = grupoengre['d_p']
    b_face = grupoengre['b_face']
    m = grupoengre['m']
    N_p, N_c = grupoengre['N_p'], grupoengre['N_c']
    phi_rad = grupoengre['phi_rad']
    
    # Fatores AGMA conservadores para off-road
    K_v = calcula_Kv_AGMA(V_t, 6, phi_rad)
    K_m = calcula_Km_AGMA(b_face, d_p, 'comercial', 'Q6')
    
    # Fatores de segurança aumentados para off-road
    K_s = 1.1
    K_R = 1.5  # Confiabilidade para off-road
    C_R = 1.5
    K_T = 1.0
    J = 0.32
    
    # Análise de Flexão
    sigma_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
    sigma_b_adm = S_at / (K_T * K_R)
    FS_flexao = sigma_b_adm / sigma_b if sigma_b > 0 else float('inf')

    # Análise de Pitting
    mg = N_c / N_p
    I_geo = (np.sin(phi_rad) * np.cos(phi_rad) / 2) * (mg / (mg + 1))
    pitting = (W_t * K_s * K_m) / (b_face * d_p * I_geo * K_v)
    sigma_c = C_p * np.sqrt(pitting) if pitting > 0 else 0
    sigma_c_adm = S_ac / C_R
    FS_pitting = (sigma_c_adm / sigma_c)**2 if sigma_c > 0 else float('inf')
    
    # Vida útil profissional
    vida_util = calcula_vida_util_profissional(
        grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_horas
    )
    
    # Critério de aprovação profissional
    FS_min_flexao = 2.0
    FS_min_pitting = 2.0
    
    if (FS_flexao >= FS_min_flexao and 
        FS_pitting >= FS_min_pitting and
        vida_util['status_vida'] == "APROVADO" and
        vida_util['vida_util_minima'] >= vida_util_horas):
        status_completo = "APROVADO"
    elif (FS_flexao >= 1.7 and 
          FS_pitting >= 1.7 and
          vida_util['vida_util_minima'] >= vida_util_horas * 0.8):
        status_completo = "MARGINAL"
    else:
        status_completo = "REPROVADO"
    
    return {
        "FS_flexao": FS_flexao,
        "FS_pitting": FS_pitting,
        'sigma_b': sigma_b,
        'sigma_c': sigma_c,
        'K_m': K_m,
        'K_v': K_v,
        'status_completo': status_completo,
        **vida_util
    }

def debug_estagio1(guardaengrenagem1):
    """Analisa por que nenhuma configuração foi aprovada no Estágio 1"""
    
    if not guardaengrenagem1:
        print("DEBUG: Nenhuma engrenagem foi testada no Estágio 1")
        return
    
    # Analisar os primeiros 5 resultados para entender o problema
    print("\n=== DEBUG ESTÁGIO 1 ===")
    print("Analisando as primeiras 5 configurações testadas:")
    
    for i, config in enumerate(guardaengrenagem1[:5]):
        print(f"\nConfiguração {i+1}:")
        print(f"  Np: {config['N_p']}, Nc: {config['N_c']}, m: {config['m']}mm")
        print(f"  FS Flexão: {config['FS_flexao']:.2f}")
        print(f"  FS Pitting: {config['FS_pitting']:.2f}")
        print(f"  Vida Flexão: {config['vida_flexao_horas']:.0f}h")
        print(f"  Vida Pitting: {config['vida_pitting_horas']:.0f}h")
        print(f"  Status: {config['status_completo']}")
        
        # Identificar qual critério está falhando
        problemas = []
        if config['FS_flexao'] < 1.5:
            problemas.append(f"FS Flexão baixo ({config['FS_flexao']:.2f})")
        if config['FS_pitting'] < 1.5:
            problemas.append(f"FS Pitting baixo ({config['FS_pitting']:.2f})")
        if config['vida_flexao_horas'] < 4000:
            problemas.append(f"Vida flexão insuficiente ({config['vida_flexao_horas']:.0f}h)")
        if config['vida_pitting_horas'] < 4000:
            problemas.append(f"Vida pitting insuficiente ({config['vida_pitting_horas']:.0f}h)")
        
        if problemas:
            print(f"  PROBLEMAS: {', '.join(problemas)}")

# def analisa_fatores_com_vida(grupoengre, forcas, S_at, S_ac, C_p, vida_util_horas=4000, montagem_tipo='comercial', precisao_engrenagem='Q6'):
#     """
#     Análise de fatores AGMA com verificação de vida útil
#     """
    
#     # Recuperando dados
#     W_t = forcas['W_t']
#     V_t = forcas['V_t']
#     d_p = grupoengre['d_p']
#     b_face = grupoengre['b_face']
#     m = grupoengre['m']
#     N_p, N_c = grupoengre['N_p'], grupoengre['N_c']
#     phi_rad = grupoengre['phi_rad']
    
#     # --- Fatores AGMA CORRETOS ---
    
#     # Fator Dinâmico (K_v) - AGMA
#     K_v = calcula_Kv_AGMA(V_t, 6, phi_rad)  # Qv = 6 para aplicação comercial
    
#     # Fator de Distribuição de Carga (K_m) - AGMA
#     K_m = calcula_Km_AGMA(b_face, d_p, montagem_tipo, precisao_engrenagem)
    
#     # Outros fatores
#     K_s = 1.0   # Fator de Tamanho
#     K_R = 1.25  # Fator de Confiabilidade (Flexão, off-road)
#     C_R = 1.25  # Fator de Confiabilidade (Superfície)
#     K_T = 1.0   # Fator de Temperatura
#     J = 0.32    # Fator de Geometria
    
#     # Analise de Flexão
#     sigma_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
#     sigma_b_adm = S_at / (K_T * K_R)

#     FS_flexao = np.inf if sigma_b == 0 else sigma_b_adm / sigma_b

#     # Analise de Contato/Pitting
#     mg = N_c / N_p

#     # Fator de Geometria (Superfície)
#     I_geo = (np.sin(phi_rad) * np.cos(phi_rad) / 2) * (mg / (mg + 1))

#     # Desgaste dos dentes das engrenagens
#     pitting = (W_t * K_s * K_m) / (b_face * d_p * I_geo * K_v)

#     sigma_c = np.inf
#     if pitting >= 0 and (b_face * d_p * I_geo * K_v) != 0:
#         sigma_c = C_p * np.sqrt(pitting)

#     sigma_c_adm = S_ac / C_R

#     FS_pitting = 0.0
#     if sigma_c == np.inf:
#         FS_pitting = 0.0
#     elif sigma_c == 0:
#         FS_pitting = np.inf
#     else:
#         FS_pitting = (sigma_c_adm / sigma_c)**2
    
#     # VERIFICAÇÃO DE VIDA ÚTIL
#     vida_util = calcula_vida_util_engrenagens(
#         grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_horas
#     )
    
#     # Status combinado (FS + Vida Útil)
#     if (FS_flexao > 1.5 and FS_pitting > 1.5 and 
#         vida_util['status_vida'] == "APROVADO"):
#         status_completo = "APROVADO"
#     elif (FS_flexao > 1.2 and FS_pitting > 1.2 and 
#           vida_util['status_vida'] in ["APROVADO", "MARGINAL"]):
#         status_completo = "MARGINAL"
#     else:
#         status_completo = "REPROVADO"

#     return {
#         "FS_flexao": FS_flexao, 
#         "FS_pitting": FS_pitting,
#         'sigma_b': sigma_b,
#         'sigma_c': sigma_c,
#         'K_m': K_m,
#         'K_v': K_v,
#         'status_completo': status_completo,
#         **vida_util  # Inclui todos os dados de vida útil
#     }

# def analisa_fatores_com_vida(grupoengre, forcas, S_at, S_ac, C_p, vida_util_horas=4000, montagem_tipo='comercial', precisao_engrenagem='Q6'):
#     """
#     Análise de fatores AGMA com verificação de vida útil CORRIGIDA
#     """
    
#     # Recuperando dados
#     W_t = forcas['W_t']
#     V_t = forcas['V_t']
#     d_p = grupoengre['d_p']
#     b_face = grupoengre['b_face']
#     m = grupoengre['m']
#     N_p, N_c = grupoengre['N_p'], grupoengre['N_c']
#     phi_rad = grupoengre['phi_rad']
    
#     # --- Fatores AGMA CORRETOS ---
#     K_v = calcula_Kv_AGMA(V_t, 6, phi_rad)
#     K_m = calcula_Km_AGMA(b_face, d_p, montagem_tipo, precisao_engrenagem)
    
#     # Outros fatores
#     K_s = 1.0
#     K_R = 1.25
#     C_R = 1.25
#     K_T = 1.0
#     J = 0.32
    
#     # Analise de Flexão
#     sigma_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
#     sigma_b_adm = S_at / (K_T * K_R)
#     FS_flexao = np.inf if sigma_b == 0 else sigma_b_adm / sigma_b

#     # Analise de Contato/Pitting
#     mg = N_c / N_p
#     I_geo = (np.sin(phi_rad) * np.cos(phi_rad) / 2) * (mg / (mg + 1))
#     pitting = (W_t * K_s * K_m) / (b_face * d_p * I_geo * K_v)

#     sigma_c = np.inf
#     if pitting >= 0 and (b_face * d_p * I_geo * K_v) != 0:
#         sigma_c = C_p * np.sqrt(pitting)

#     sigma_c_adm = S_ac / C_R
#     FS_pitting = 0.0
#     if sigma_c == np.inf:
#         FS_pitting = 0.0
#     elif sigma_c == 0:
#         FS_pitting = np.inf
#     else:
#         FS_pitting = (sigma_c_adm / sigma_c)**2
    
#     # VERIFICAÇÃO DE VIDA ÚTIL CORRIGIDA
#     vida_util = calcula_vida_util_engrenagens(
#         grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_horas
#     )
    
#     # Status combinado (FS + Vida Útil)
#     if (FS_flexao > 1.5 and FS_pitting > 1.5 and 
#         vida_util['status_vida'] == "APROVADO"):
#         status_completo = "APROVADO"
#     elif (FS_flexao > 1.2 and FS_pitting > 1.2 and 
#           vida_util['status_vida'] in ["APROVADO", "MARGINAL"]):
#         status_completo = "MARGINAL"
#     else:
#         status_completo = "REPROVADO"

#     return {
#         "FS_flexao": FS_flexao, 
#         "FS_pitting": FS_pitting,
#         'sigma_b': sigma_b,
#         'sigma_c': sigma_c,
#         'K_m': K_m,
#         'K_v': K_v,
#         'status_completo': status_completo,
#         **vida_util
#     }

def analisa_fatores_com_vida(grupoengre, forcas, S_at, S_ac, C_p, vida_util_horas=4000, montagem_tipo='comercial', precisao_engrenagem='Q6'):
    """
    Análise com vida útil REALISTA e fatores de segurança apropriados para off-road
    """
    
    # Recuperando dados
    W_t = forcas['W_t']
    V_t = forcas['V_t']
    d_p = grupoengre['d_p']
    b_face = grupoengre['b_face']
    m = grupoengre['m']
    N_p, N_c = grupoengre['N_p'], grupoengre['N_c']
    phi_rad = grupoengre['phi_rad']
    
    # --- Fatores AGMA COM VALORES MAIS CONSERVADORES PARA OFF-ROAD ---
    K_v = calcula_Kv_AGMA(V_t, 6, phi_rad)
    K_m = calcula_Km_AGMA(b_face, d_p, montagem_tipo, precisao_engrenagem)
    
    # Fatores mais conservadores para aplicação off-road
    K_s = 1.2   # Aumentado para off-road
    K_R = 1.25   # Aumentado para off-road (era 1.25)
    C_R = 1.25   # Aumentado para off-road (era 1.25)
    K_T = 1.0
    J = 0.32
    
    # Analise de Flexão
    sigma_b = (W_t / (b_face * m * J)) * (K_s * K_m) / K_v
    sigma_b_adm = S_at / (K_T * K_R)
    FS_flexao = np.inf if sigma_b == 0 else sigma_b_adm / sigma_b

    # Analise de Contato/Pitting
    mg = N_c / N_p
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
    
    # VERIFICAÇÃO DE VIDA ÚTIL REALISTA
    vida_util = calcula_vida_util_engrenagens(
        grupoengre, forcas, sigma_b, sigma_c, S_at, S_ac, vida_util_horas
    )
    
    # Status combinado MAIS CONSERVADOR
    # Para off-road, exigir FS mínimo de 2.0 e vida com margem
    if (FS_flexao > 2.0 and FS_pitting > 2.0 and 
        vida_util['status_vida'] == "APROVADO" and
        vida_util['vida_util_conservadora'] >= vida_util_horas * 1.5):  # 50% de margem
        status_completo = "APROVADO"
    elif (FS_flexao > 1.5 and FS_pitting > 1.5 and 
          vida_util['vida_util_conservadora'] >= vida_util_horas):
        status_completo = "MARGINAL"
    else:
        status_completo = "REPROVADO"

    return {
        "FS_flexao": FS_flexao, 
        "FS_pitting": FS_pitting,
        'sigma_b': sigma_b,
        'sigma_c': sigma_c,
        'K_m': K_m,
        'K_v': K_v,
        'status_completo': status_completo,
        **vida_util
    }

def dimensionamento_engrenagens_com_vida(T_in, n_in, i_estagio, eta_estagio, S_at, S_ac, C_p, vida_util_horas=4000):
    """
    Dimensionamento automático com verificação de vida útil - CORRIGIDA
    """
    
    # Inicializar todas as variáveis no início
    guardaengrenagem1 = []
    guardaengrenagem2 = []
    modulos_padronizados = [1.25,1.5, 2, 2.5, 3, 4, 5, 6, 8]
    N_pinhao_padronizado = [17, 18, 19, 20, 21, 22, 23, 24, 25]
    resultado_estagio_1 = None
    resultado_estagio_2 = None  # ← INICIALIZAR AQUI
    
    # Estimativas iniciais para K_m e K_v
    K_m_estimado = 1.5
    K_v_estimado = 1.2
    J_estimado = 0.32
    
    # ESTÁGIO 1
    for Np in N_pinhao_padronizado:
        for m in modulos_padronizados:
            # Estimar diâmetro primitivo para cálculo da força tangencial
            d_p_estimado = m * Np
            W_t_estimado = estima_W_t_inicial(T_in, d_p_estimado)
            
            # Calcular grupo de engrenagens com largura automática
            grupo_engrenagens = calcula_grupo_engrenagens_auto(
                i_estagio, m, Np, W_t_estimado, S_at, 
                K_m_estimado, K_v_estimado, J_estimado
            )
            
            # Análise de forças e fatores com verificação de vida útil
            forcas = calcula_forcas(T_in, n_in, grupo_engrenagens, eta_estagio)
            forcas['n_p_in'] = n_in
            
            fatores = analisa_fatores_com_vida(
                grupo_engrenagens, forcas, S_at, S_ac, C_p, vida_util_horas
            )
            
            dados = {**grupo_engrenagens, **forcas, **fatores}
            dados['erro_estagio'] = abs((dados['i_efetiva'] - i_estagio) / i_estagio)
            
            if dados['status_completo'] == "APROVADO":
                dados['status_calc'] = "Aprovado"
                guardaengrenagem1.append(dados)
            elif dados['status_completo'] == "MARGINAL":
                dados['status_calc'] = "Marginal"
                guardaengrenagem1.append(dados)
            else:
                dados['status_calc'] = "Reprovado"
                guardaengrenagem1.append(dados)
    
    # Seleção do melhor par (priorizar "Aprovado")
    EngrenagemAprovada = [p for p in guardaengrenagem1 if p['status_calc'] == "Aprovado"]
    EngrenagemMarginal = [p for p in guardaengrenagem1 if p['status_calc'] == "Marginal"]
    
    if EngrenagemAprovada:
        EngrenagemAprovada.sort(key=lambda x: x['erro_estagio'])
        resultado_estagio_1 = EngrenagemAprovada[0]
    elif EngrenagemMarginal:
        EngrenagemMarginal.sort(key=lambda x: x['erro_estagio'])
        resultado_estagio_1 = EngrenagemMarginal[0]
        print("AVISO: Estágio 1 usando configuração MARGINAL")
    else:
        resultado_estagio_1 = None
        print("AVISO: Nenhuma configuração adequada para Estágio 1")
        # ← AGORA RETORNA resultado_estagio_2 = None JÁ DEFINIDO
    
    # Se Estágio 1 falhou, retorne imediatamente
    if resultado_estagio_1 is None:
        return resultado_estagio_1, resultado_estagio_2, guardaengrenagem1, guardaengrenagem2
    
    # ESTÁGIO 2 (só executa se Estágio 1 foi bem-sucedido)
    T_in_2 = resultado_estagio_1['T_p_out']
    n_in_2 = resultado_estagio_1['n_p_out']
    
    for Np in N_pinhao_padronizado:
        for m in modulos_padronizados:
            d_p_estimado = m * Np
            W_t_estimado = estima_W_t_inicial(T_in_2, d_p_estimado)
            
            grupo_engrenagens = calcula_grupo_engrenagens_auto(
                i_estagio, m, Np, W_t_estimado, S_at,
                K_m_estimado, K_v_estimado, J_estimado
            )
            
            forcas = calcula_forcas(T_in_2, n_in_2, grupo_engrenagens, eta_estagio)
            forcas['n_p_in'] = n_in_2
            
            fatores = analisa_fatores_com_vida(
                grupo_engrenagens, forcas, S_at, S_ac, C_p, vida_util_horas
            )
            
            dados = {**grupo_engrenagens, **forcas, **fatores}
            dados['erro_estagio'] = abs((dados['i_efetiva'] - i_estagio) / i_estagio)
            
            if dados['status_completo'] == "APROVADO":
                dados['status_calc'] = "Aprovado"
                guardaengrenagem2.append(dados)
            elif dados['status_completo'] == "MARGINAL":
                dados['status_calc'] = "Marginal"
                guardaengrenagem2.append(dados)
            else:
                dados['status_calc'] = "Reprovado"
                guardaengrenagem2.append(dados)
    
    # Seleção do melhor par estágio 2
    EngrenagemAprovada2 = [p for p in guardaengrenagem2 if p['status_calc'] == "Aprovado"]
    EngrenagemMarginal2 = [p for p in guardaengrenagem2 if p['status_calc'] == "Marginal"]
    
    if EngrenagemAprovada2:
        EngrenagemAprovada2.sort(key=lambda x: x['erro_estagio'])
        resultado_estagio_2 = EngrenagemAprovada2[0]
    elif EngrenagemMarginal2:
        EngrenagemMarginal2.sort(key=lambda x: x['erro_estagio'])
        resultado_estagio_2 = EngrenagemMarginal2[0]
        print("AVISO: Estágio 2 usando configuração MARGINAL")
    else:
        resultado_estagio_2 = None
        print("AVISO: Nenhuma configuração adequada para Estágio 2")
    
    return resultado_estagio_1, resultado_estagio_2, guardaengrenagem1, guardaengrenagem2

def calcula_reações_mancais_eixo1(forcas_engrenagem1, distancias):
    """
    Calcula reações nos mancais do Eixo 1 (entrada)
    """
    W_t = forcas_engrenagem1['W_t']
    W_r = forcas_engrenagem1['W_r']
    
    LA = distancias['LA']  # Distância mancal A até engrenagem
    AB = distancias['AB']  # Distância entre mancais A e B
    BC = distancias['BC']  # Distância engrenagem até mancal B
    
    # Reações no plano vertical (forças radiais)
    R_Av = (W_r * BC) / AB
    R_Bv = (W_r * LA) / AB
    
    # Reações no plano horizontal (forças tangenciais)  
    R_Ah = (W_t * BC) / AB
    R_Bh = (W_t * LA) / AB
    
    # Reações resultantes
    R_A = np.sqrt(R_Av**2 + R_Ah**2)
    R_B = np.sqrt(R_Bv**2 + R_Bh**2)
    
    return {
        'R_A_vertical': R_Av,
        'R_A_horizontal': R_Ah, 
        'R_A_resultante': R_A,
        'R_B_vertical': R_Bv,
        'R_B_horizontal': R_Bh,
        'R_B_resultante': R_B
    }

def calcula_reações_eixo2(forcas_coroa1, forcas_pinhao2, distancias):
    """
    Calcula reações nos mancais do Eixo 2 (intermediário)
    - Recebe forças da Coroa 1 e Pinhão 2
    """
    W_t1 = forcas_coroa1['W_t']
    W_r1 = forcas_coroa1['W_r']
    W_t2 = forcas_pinhao2['W_t'] 
    W_r2 = forcas_pinhao2['W_r']
    
    LA = distancias['LA']  # Mancal A até Coroa 1
    L1 = distancias['L1']  # Coroa 1 até Pinhão 2
    L2 = distancias['L2']  # Pinhão 2 até Mancal B
    
    AB = LA + L1 + L2  # Distância total entre mancais
    
    # Reações no plano vertical
    R_Av = (W_r1 * (L1 + L2) + W_r2 * L2) / AB
    R_Bv = (W_r1 * LA + W_r2 * (LA + L1)) / AB
    
    # Reações no plano horizontal
    R_Ah = (W_t1 * (L1 + L2) + W_t2 * L2) / AB
    R_Bh = (W_t1 * LA + W_t2 * (LA + L1)) / AB
    
    # Reações resultantes
    R_A = np.sqrt(R_Av**2 + R_Ah**2)
    R_B = np.sqrt(R_Bv**2 + R_Bh**2)
    
    return {
        'R_A_vertical': R_Av, 'R_A_horizontal': R_Ah, 'R_A_resultante': R_A,
        'R_B_vertical': R_Bv, 'R_B_horizontal': R_Bh, 'R_B_resultante': R_B
    }
import numpy as np
import matplotlib.pyplot as plt

def calcula_diagramas_eixo_simples(forcas_engrenagem, distancias, torque, nome_eixo):
    """
    Calcula diagramas completos para eixo com uma engrenagem
    Retorna valores críticos e plota diagramas
    """
    W_t = forcas_engrenagem['W_t']
    W_r = forcas_engrenagem['W_r']
    
    LA = distancias['LA']
    AB = distancias['AB']
    BC = distancias['BC']
    
    # 1. CÁLCULO DAS REAÇÕES
    R_Av = (W_r * BC) / AB
    R_Bv = (W_r * LA) / AB
    R_Ah = (W_t * BC) / AB  
    R_Bh = (W_t * LA) / AB
    
    # 2. DIAGRAMAS DE ESFORÇO CORTANTE
    # Pontos ao longo do eixo (de 0 a AB)
    x = np.linspace(0, AB, 100)
    
    # Cortante Vertical
    V_v = np.zeros_like(x)
    V_v[x <= LA] = R_Av
    V_v[x > LA] = -R_Bv
    
    # Cortante Horizontal  
    V_h = np.zeros_like(x)
    V_h[x <= LA] = R_Ah
    V_h[x > LA] = -R_Bh
    
    # Cortante Resultante
    V_resultante = np.sqrt(V_v**2 + V_h**2)
    
    # 3. DIAGRAMAS DE MOMENTO FLETOR
    # Momento Vertical
    M_v = np.zeros_like(x)
    M_v[x <= LA] = R_Av * x[x <= LA]
    M_v[x > LA] = R_Av * x[x > LA] - W_r * (x[x > LA] - LA)
    
    # Momento Horizontal
    M_h = np.zeros_like(x)  
    M_h[x <= LA] = R_Ah * x[x <= LA]
    M_h[x > LA] = R_Ah * x[x > LA] - W_t * (x[x > LA] - LA)
    
    # Momento Resultante
    M_resultante = np.sqrt(M_v**2 + M_h**2)
    
    # 4. VALORES CRÍTICOS
    M_max = np.max(M_resultante)
    V_max = np.max(V_resultante)
    posicao_M_max = x[np.argmax(M_resultante)]
    
    # 5. PLOTAR DIAGRAMAS
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.plot(x, V_v, 'b-', linewidth=2, label='Cortante Vertical')
    plt.plot(x, V_h, 'r-', linewidth=2, label='Cortante Horizontal')
    plt.plot(x, V_resultante, 'k--', linewidth=2, label='Cortante Resultante')
    plt.axvline(x=LA, color='g', linestyle=':', alpha=0.7, label='Engrenagem')
    plt.grid(True, alpha=0.3)
    plt.title(f'{nome_eixo} - Esforço Cortante')
    plt.xlabel('Posição (mm)')
    plt.ylabel('Força (N)')
    plt.legend()
    
    plt.subplot(2, 2, 2)
    plt.plot(x, M_v, 'b-', linewidth=2, label='Momento Vertical')
    plt.plot(x, M_h, 'r-', linewidth=2, label='Momento Horizontal') 
    plt.plot(x, M_resultante, 'k-', linewidth=3, label='Momento Resultante')
    plt.axvline(x=LA, color='g', linestyle=':', alpha=0.7, label='Engrenagem')
    plt.axvline(x=posicao_M_max, color='red', linestyle='--', alpha=0.7, 
                label=f'M_max = {M_max:.1f} N.mm')
    plt.grid(True, alpha=0.3)
    plt.title(f'{nome_eixo} - Momento Fletor')
    plt.xlabel('Posição (mm)')
    plt.ylabel('Momento (N.mm)')
    plt.legend()
    
    plt.subplot(2, 2, 3)
    # Diagrama de Torque
    torque_Nmm = torque * 1000  # Converter N.m para N.mm
    T_plot = np.full_like(x, torque_Nmm)
    plt.plot(x, T_plot, 'purple', linewidth=3, label=f'Torque = {torque:.1f} N.m')
    plt.axvline(x=LA, color='g', linestyle=':', alpha=0.7, label='Engrenagem')
    plt.grid(True, alpha=0.3)
    plt.title(f'{nome_eixo} - Torque')
    plt.xlabel('Posição (mm)')
    plt.ylabel('Torque (N.mm)')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'diagramas_{nome_eixo.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return {
        'M_max': M_max,
        'V_max': V_max,
        'posicao_M_max': posicao_M_max,
        'R_A_vertical': R_Av,
        'R_A_horizontal': R_Ah,
        'R_A_resultante': np.sqrt(R_Av**2 + R_Ah**2),
        'R_B_vertical': R_Bv,
        'R_B_horizontal': R_Bh, 
        'R_B_resultante': np.sqrt(R_Bv**2 + R_Bh**2)
    }

def calcula_diagramas_eixo_duplo(forcas_engrenagem1, forcas_engrenagem2, distancias, torque, nome_eixo):
    """
    Calcula diagramas completos para eixo com duas engrenagens (Eixo 2)
    """
    W_t1 = forcas_engrenagem1['W_t']
    W_r1 = forcas_engrenagem1['W_r']
    W_t2 = forcas_engrenagem2['W_t']
    W_r2 = forcas_engrenagem2['W_r']
    
    LA = distancias['LA']  # Mancal A até Engrenagem 1
    L1 = distancias['L1']  # Engrenagem 1 até Engrenagem 2  
    L2 = distancias['L2']  # Engrenagem 2 até Mancal B
    AB = LA + L1 + L2
    
    # 1. CÁLCULO DAS REAÇÕES
    R_Av = (W_r1 * (L1 + L2) + W_r2 * L2) / AB
    R_Bv = (W_r1 * LA + W_r2 * (LA + L1)) / AB
    R_Ah = (W_t1 * (L1 + L2) + W_t2 * L2) / AB
    R_Bh = (W_t1 * LA + W_t2 * (LA + L1)) / AB
    
    # 2. DIAGRAMAS DE ESFORÇO CORTANTE
    x = np.linspace(0, AB, 200)
    
    # Cortante Vertical
    V_v = np.zeros_like(x)
    mask1 = x <= LA
    mask2 = (x > LA) & (x <= LA + L1)
    mask3 = x > LA + L1
    
    V_v[mask1] = R_Av
    V_v[mask2] = R_Av - W_r1
    V_v[mask3] = R_Av - W_r1 - W_r2
    
    # Cortante Horizontal
    V_h = np.zeros_like(x)
    V_h[mask1] = R_Ah
    V_h[mask2] = R_Ah - W_t1  
    V_h[mask3] = R_Ah - W_t1 - W_t2
    
    V_resultante = np.sqrt(V_v**2 + V_h**2)
    
    # 3. DIAGRAMAS DE MOMENTO FLETOR
    M_v = np.zeros_like(x)
    M_v[mask1] = R_Av * x[mask1]
    M_v[mask2] = R_Av * x[mask2] - W_r1 * (x[mask2] - LA)
    M_v[mask3] = R_Av * x[mask3] - W_r1 * (x[mask3] - LA) - W_r2 * (x[mask3] - LA - L1)
    
    M_h = np.zeros_like(x)
    M_h[mask1] = R_Ah * x[mask1]
    M_h[mask2] = R_Ah * x[mask2] - W_t1 * (x[mask2] - LA)
    M_h[mask3] = R_Ah * x[mask3] - W_t1 * (x[mask3] - LA) - W_t2 * (x[mask3] - LA - L1)
    
    M_resultante = np.sqrt(M_v**2 + M_h**2)
    
    # 4. VALORES CRÍTICOS
    M_max = np.max(M_resultante)
    V_max = np.max(V_resultante)
    posicao_M_max = x[np.argmax(M_resultante)]
    
    # 5. PLOTAR DIAGRAMAS
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.plot(x, V_v, 'b-', linewidth=2, label='Cortante Vertical')
    plt.plot(x, V_h, 'r-', linewidth=2, label='Cortante Horizontal')
    plt.plot(x, V_resultante, 'k--', linewidth=2, label='Cortante Resultante')
    plt.axvline(x=LA, color='g', linestyle=':', alpha=0.7, label='Engrenagem 1')
    plt.axvline(x=LA+L1, color='orange', linestyle=':', alpha=0.7, label='Engrenagem 2')
    plt.grid(True, alpha=0.3)
    plt.title(f'{nome_eixo} - Esforço Cortante')
    plt.xlabel('Posição (mm)')
    plt.ylabel('Força (N)')
    plt.legend()
    
    plt.subplot(2, 2, 2)
    plt.plot(x, M_v, 'b-', linewidth=2, label='Momento Vertical')
    plt.plot(x, M_h, 'r-', linewidth=2, label='Momento Horizontal')
    plt.plot(x, M_resultante, 'k-', linewidth=3, label='Momento Resultante')
    plt.axvline(x=LA, color='g', linestyle=':', alpha=0.7, label='Engrenagem 1')
    plt.axvline(x=LA+L1, color='orange', linestyle=':', alpha=0.7, label='Engrenagem 2')
    plt.axvline(x=posicao_M_max, color='red', linestyle='--', alpha=0.7,
                label=f'M_max = {M_max:.1f} N.mm')
    plt.grid(True, alpha=0.3)
    plt.title(f'{nome_eixo} - Momento Fletor')
    plt.xlabel('Posição (mm)')
    plt.ylabel('Momento (N.mm)')
    plt.legend()
    
    plt.subplot(2, 2, 3)
    # Diagrama de Torque
    torque_Nmm = torque * 1000
    T_plot = np.full_like(x, torque_Nmm)
    plt.plot(x, T_plot, 'purple', linewidth=3, label=f'Torque = {torque:.1f} N.m')
    plt.axvline(x=LA, color='g', linestyle=':', alpha=0.7, label='Engrenagem 1')
    plt.axvline(x=LA+L1, color='orange', linestyle=':', alpha=0.7, label='Engrenagem 2')
    plt.grid(True, alpha=0.3)
    plt.title(f'{nome_eixo} - Torque')
    plt.xlabel('Posição (mm)')
    plt.ylabel('Torque (N.mm)')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(f'diagramas_{nome_eixo.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return {
        'M_max': M_max,
        'V_max': V_max, 
        'posicao_M_max': posicao_M_max,
        'R_A_vertical': R_Av,
        'R_A_horizontal': R_Ah,
        'R_A_resultante': np.sqrt(R_Av**2 + R_Ah**2),
        'R_B_vertical': R_Bv,
        'R_B_horizontal': R_Bh,
        'R_B_resultante': np.sqrt(R_Bv**2 + R_Bh**2)
    }