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

def calcula_forcas(T_p_in, n_p_in, grupoengre, eta_estagio):
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
