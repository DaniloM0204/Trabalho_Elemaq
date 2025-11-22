import numpy as np
import matplotlib.pyplot as plt

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
        if fator > 1.0:
            fator = 1.0

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
    plt.savefig(f'Plots/diagramas_{nome_eixo.replace(" ", "_")}_momento.png', dpi=300, bbox_inches='tight')

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
    plt.plot(x, T_plot, 'purple', linewidth=3, label=f'Torque = {torque:.2f} N.m')
    plt.axvline(x=LA, color='g', linestyle=':', alpha=0.7, label='Engrenagem 1')
    plt.axvline(x=LA+L1, color='orange', linestyle=':', alpha=0.7, label='Engrenagem 2')
    plt.grid(True, alpha=0.3)
    plt.title(f'{nome_eixo} - Torque')
    plt.xlabel('Posição (mm)')
    plt.ylabel('Torque (N.mm)')
    plt.legend()

    plt.tight_layout()
    plt.savefig(f'Plots/diagramas_{nome_eixo.replace(" ", "_")}.png', dpi=300, bbox_inches='tight')

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

def dimensiona_eixo_por_fadiga(M_max, Tm, Sut, Sy, Se_linha, vida_util_horas=4000, Nf=2.0, tipo_eixo="simples"):
    """
    Dimensiona eixo por fadiga usando critério de Goodman modificado
    Considera fatores de concentração de tensão para chavetas
    """

    # Converter unidades
    M_max_Nmm = M_max  # Já está em N.mm
    Tm_Nmm = Tm * 1000  # Converter N.m para N.mm

    # Fatores de concentração de tensão para chavetas (valores conservadores)
    if tipo_eixo == "simples":
        Kt_flexao = 2.0  # Para entalhe de chaveta
        Kt_torcao = 1.7  # Para entalhe de chaveta
    else:
        Kt_flexao = 2.2  # Para eixo com múltiplas chavetas
        Kt_torcao = 1.8

    # Sensibilidade ao entalhe (para aço com Sut = 469 MPa)
    q = 0.85  # Sensibilidade média
    Kf_flexao = 1 + q * (Kt_flexao - 1)
    Kf_torcao = 1 + q * (Kt_torcao - 1)

    # Fatores de correção do limite de fadiga
    C_carreg = 1.0      # Flexão
    C_temp = 1.0        # Temperatura ambiente
    C_conf = 0.814      # 99% confiabilidade

    # Loop iterativo para encontrar diâmetro
    d_guess = 20.0  # mm - chute inicial
    tolerance = 0.01
    max_iter = 50
    iter_count = 0

    print(f"Dimensionando eixo {tipo_eixo}...")
    print(f"M_max = {M_max:.0f} N.mm, Tm = {Tm:.1f} N.m")

    while iter_count < max_iter:
        # Fator de tamanho
        if d_guess <= 8:
            C_tam = 1.0
        elif d_guess <= 250:
            C_tam = 1.189 * (d_guess ** -0.097)
        else:
            C_tam = 0.6

        # Fator de superfície
        A = 4.51  # Usinado
        b = -0.265
        C_sup = A * (Sut ** b)
        C_sup = min(C_sup, 1.0)

        # Limite de fadiga corrigido
        Se = Se_linha * C_carreg * C_tam * C_sup * C_temp * C_conf

        # Tensões nominais
        sigma_a_nominal = (32 * M_max_Nmm) / (np.pi * d_guess**3)
        tau_m_nominal = (16 * Tm_Nmm) / (np.pi * d_guess**3)

        # Tensões efetivas com concentradores
        sigma_a = Kf_flexao * sigma_a_nominal
        tau_m = Kf_torcao * tau_m_nominal

        # Critério de Goodman modificado para eixos
        # (σ_a / Se) + (τ_m / Sut) ≤ 1/Nf
        left_side = (sigma_a / Se) + (tau_m / Sut)

        # Verificar critério
        if abs(left_side - 1/Nf) < tolerance:
            break

        # Atualizar diâmetro
        if left_side > 1/Nf:
            # Muito estressado - aumentar diâmetro
            d_guess *= 1.05
        else:
            # Superdimensionado - diminuir diâmetro
            d_guess *= 0.98

        iter_count += 1

    # Tensões finais
    sigma_a_final = Kf_flexao * (32 * M_max_Nmm) / (np.pi * d_guess**3)
    tau_m_final = Kf_torcao * (16 * Tm_Nmm) / (np.pi * d_guess**3)

    # Fator de segurança final
    FS_final = 1 / ((sigma_a_final / Se) + (tau_m_final / Sut))

    # Verificação de escoamento estático
    sigma_max = (32 * M_max_Nmm) / (np.pi * d_guess**3) + (16 * Tm_Nmm) / (np.pi * d_guess**3)
    FS_escoamento = Sy / sigma_max if sigma_max > 0 else float('inf')

    return {
        'diametro_minimo': d_guess,
        'FS_fadiga': FS_final,
        'FS_escoamento': FS_escoamento,
        'sigma_a_efetiva': sigma_a_final,
        'tau_m_efetiva': tau_m_final,
        'Se_corrigido': Se,
        'iteracoes': iter_count
    }

def analisa_deflexao_eixo(M_max, comprimento, diametro, material='aco'):
    """
    Analisa deflexão máxima do eixo
    """
    E = 200000  # Módulo de elasticidade do aço (MPa)

    # Momento de inércia
    Inercia = (np.pi * diametro**4) / 64

    # Deflexão máxima aproximada para viga bi-apoiada com carga no centro
    # δ_max = (F * L^3) / (48 * E * I)
    # Mas temos momento fletor, então usamos relação M = F * L / 4 => F = 4 * M / L
    F_equivalente = (4 * M_max) / comprimento  # M_max em N.mm, comprimento em mm
    deflexao_max = (F_equivalente * comprimento**3) / (48 * E * Inercia)

    # Critério: deflexão ≤ 0.001 * comprimento (critério industrial)
    deflexao_admissivel = 0.001 * comprimento
    status_deflexao = "APROVADO" if deflexao_max <= deflexao_admissivel else "REPROVADO"

    return {
        'deflexao_max_mm': deflexao_max,
        'deflexao_admissivel_mm': deflexao_admissivel,
        'status_deflexao': status_deflexao,
        'relacao_deflexao': deflexao_max / comprimento
    }

def calcula_velocidade_critica(diametro, comprimento, material='aco'):
    """
    Calcula velocidade crítica do eixo (primeira frequência natural)
    """
    E = 200000  # MPa
    rho = 7850  # kg/m³ - densidade do aço

    # Momento de inércia
    diametro_m = diametro / 1000  # Converter para metros
    Inercia = (np.pi * diametro_m**4) / 64

    # Área da seção transversal
    A = (np.pi * diametro_m**2) / 4

    # Comprimento em metros
    L = comprimento / 1000

    # Massa por unidade de comprimento
    m = rho * A

    # Frequência natural (Hz) - viga bi-apoiada
    # ω_n = (π² / L²) * √(EI / m)
    omega_n = (np.pi**2 / L**2) * np.sqrt(E * 1e6 * Inercia / m)  # E em Pa

    # Velocidade crítica em RPM
    N_critica = omega_n * 60 / (2 * np.pi)

    return {
        'velocidade_critica_rpm': N_critica,
        'frequencia_natural_Hz': omega_n / (2 * np.pi)
    }
