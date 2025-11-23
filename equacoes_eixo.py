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
        fator = A * (S_ut ** b)
        if fator > 1.0: fator = 1.0
        c_lista.append(fator)

    maior_valor = max(c_lista)
    indice_maior = c_lista.index(maior_valor)
    melhor_acabamento = nomes[indice_maior]

    return maior_valor, melhor_acabamento

def calcula_fator_tamanho(d):
    """Calcula C_tamanho segundo Norton para Flexão/Torção"""
    if d <= 8: return 1.0
    elif d <= 250: return 1.189 * (d ** -0.097)
    else: return 0.6

def dimensiona_eixo_por_fadiga(M_max, Tm, Sut, Sy, Se_linha, vida_util_horas=4000, Nf=2.0, tipo_eixo="simples"):
    """
    Dimensiona eixo por fadiga E escoamento simultaneamente.
    """
    M_max_Nmm = M_max
    Tm_Nmm = Tm * 1000

    if tipo_eixo == "simples":
        Kt_flexao = 2.0
        Kt_torcao = 1.7
    else:
        Kt_flexao = 2.2
        Kt_torcao = 1.8

    q = 0.75
    Kf_flexao = 1 + q * (Kt_flexao - 1)
    Kf_torcao = 1 + q * (Kt_torcao - 1)

    C_carreg = 1.0
    C_temp = 1.0
    C_conf = 0.814

    d_guess = 20.0
    tolerance = 0.01
    max_iter = 100
    iter_count = 0

    print(f"Dimensionando eixo {tipo_eixo} (M={M_max:.0f} N.mm)...")

    while iter_count < max_iter:
        C_tam = calcula_fator_tamanho(d_guess)
        A_sup, b_sup = 4.51, -0.265
        C_sup = min(A_sup * (Sut ** b_sup), 1.0)
        Se = Se_linha * C_carreg * C_tam * C_sup * C_temp * C_conf

        sigma_nominal = (32 * M_max_Nmm) / (np.pi * d_guess**3)
        tau_nominal = (16 * Tm_Nmm) / (np.pi * d_guess**3)

        sigma_a = Kf_flexao * sigma_nominal
        tau_m = 1.0 * tau_nominal

        utilizacao_fadiga = (sigma_a / Se) + (tau_m / Sut)
        sigma_von_mises = np.sqrt(sigma_nominal**2 + 3 * tau_nominal**2)
        utilizacao_escoamento = sigma_von_mises / Sy

        utilizacao_maxima = max(utilizacao_fadiga, utilizacao_escoamento)
        objetivo = 1.0 / Nf

        if abs(utilizacao_maxima - objetivo) < (tolerance * objetivo):
            break

        if utilizacao_maxima > objetivo:
            d_guess *= 1.02
        else:
            d_guess *= 0.98

        iter_count += 1

    FS_final_fadiga = 1 / ((sigma_a / Se) + (tau_m / Sut))
    sigma_vm_final = np.sqrt(((32 * M_max_Nmm) / (np.pi * d_guess**3))**2 + 3 * ((16 * Tm_Nmm) / (np.pi * d_guess**3))**2)
    FS_final_escoamento = Sy / sigma_vm_final

    return {
        'diametro_minimo': d_guess,
        'FS_fadiga': FS_final_fadiga,
        'FS_escoamento': FS_final_escoamento,
        'Se_corrigido': Se,
        'iteracoes': iter_count
    }

def analisa_deflexao_eixo(M_max, comprimento, diametro, material='aco'):
    E = 200000
    Inercia = (np.pi * diametro**4) / 64
    deflexao_max = (M_max * comprimento**2) / (10 * E * Inercia)
    deflexao_admissivel = 0.0005 * comprimento
    status_deflexao = "APROVADO" if deflexao_max <= deflexao_admissivel else "REPROVADO"

    return {
        'deflexao_max_mm': deflexao_max,
        'deflexao_admissivel_mm': deflexao_admissivel,
        'status_deflexao': status_deflexao,
    }

def calcula_velocidade_critica(diametro, comprimento, material='aco'):
    E = 200000
    rho = 7850
    diametro_m = diametro / 1000
    Inercia = (np.pi * diametro_m**4) / 64
    A = (np.pi * diametro_m**2) / 4
    L = comprimento / 1000
    m = rho * A
    omega_n = (np.pi**2 / L**2) * np.sqrt(E * 1e6 * Inercia / m)
    N_critica = omega_n * 60 / (2 * np.pi)
    return {'velocidade_critica_rpm': N_critica}

# --- FUNÇÕES DE CÁLCULO DE ESFORÇOS ---

def calcula_reações_mancais_eixo1(forcas_engrenagem1, distancias):
    W_t, W_r = forcas_engrenagem1['W_t'], forcas_engrenagem1['W_r']
    LA, AB, BC = distancias['LA'], distancias['AB'], distancias['BC']

    R_Av = (W_r * BC) / AB
    R_Bv = (W_r * LA) / AB
    R_Ah = (W_t * BC) / AB
    R_Bh = (W_t * LA) / AB

    return {'R_A_resultante': np.sqrt(R_Av**2 + R_Ah**2), 'R_B_resultante': np.sqrt(R_Bv**2 + R_Bh**2)}

def calcula_reações_eixo2(forcas_coroa1, forcas_pinhao2, distancias):
    W_t1, W_r1 = forcas_coroa1['W_t'], forcas_coroa1['W_r']
    W_t2, W_r2 = forcas_pinhao2['W_t'], forcas_pinhao2['W_r']
    LA, L1, L2 = distancias['LA'], distancias['L1'], distancias['L2']
    AB = LA + L1 + L2

    R_Av = (W_r1 * (L1 + L2) + W_r2 * L2) / AB
    R_Bv = (W_r1 * LA + W_r2 * (LA + L1)) / AB
    R_Ah = (W_t1 * (L1 + L2) + W_t2 * L2) / AB
    R_Bh = (W_t1 * LA + W_t2 * (LA + L1)) / AB

    return {'R_A_resultante': np.sqrt(R_Av**2 + R_Ah**2), 'R_B_resultante': np.sqrt(R_Bv**2 + R_Bh**2)}

def calcula_reações_mancais_eixo3(forcas_engrenagem, distancias):
    # Para o eixo 3 (saída) que só tem 1 engrenagem
    return calcula_reações_mancais_eixo1(forcas_engrenagem, distancias)

def calcula_diagramas_eixo_simples(forcas_engrenagem, distancias, torque, nome_eixo):
    W_t, W_r = forcas_engrenagem['W_t'], forcas_engrenagem['W_r']
    LA, AB, BC = distancias['LA'], distancias['AB'], distancias['BC']

    R_Av = (W_r * BC) / AB
    R_Ah = (W_t * BC) / AB

    # Momento Máximo ocorre na engrenagem (LA)
    M_v_max = R_Av * LA
    M_h_max = R_Ah * LA
    M_resultante = np.sqrt(M_v_max**2 + M_h_max**2)

    return {'M_max': M_resultante, 'V_max': 0, 'posicao_M_max': LA,
            'R_A_resultante': 0, 'R_B_resultante': 0}

def calcula_diagramas_eixo_duplo(forcas_engrenagem1, forcas_engrenagem2, distancias, torque, nome_eixo):
    W_t1, W_r1 = forcas_engrenagem1['W_t'], forcas_engrenagem1['W_r']
    W_t2, W_r2 = forcas_engrenagem2['W_t'], forcas_engrenagem2['W_r']
    LA, L1, L2 = distancias['LA'], distancias['L1'], distancias['L2']
    AB = LA + L1 + L2

    R_Av = (W_r1 * (L1 + L2) + W_r2 * L2) / AB
    R_Ah = (W_t1 * (L1 + L2) + W_t2 * L2) / AB

    # Momentos nos pontos de aplicação
    # Ponto 1 (Coroa 1): x = LA
    M_v1 = R_Av * LA
    M_h1 = R_Ah * LA
    M_res1 = np.sqrt(M_v1**2 + M_h1**2)

    # Ponto 2 (Pinhão 2): x = LA + L1
    # M = R_A * x - F1 * (x - LA)
    M_v2 = R_Av * (LA + L1) - W_r1 * L1
    M_h2 = R_Ah * (LA + L1) - W_t1 * L1
    M_res2 = np.sqrt(M_v2**2 + M_h2**2)

    M_max = max(M_res1, M_res2)

    return {'M_max': M_max, 'V_max': 0, 'posicao_M_max': 0,
            'R_A_resultante': 0, 'R_B_resultante': 0}

def plotar_diagramas_eixo(L_total, posicoes_elementos, forcas, torque, nome_eixo):
    """
    Gera e salva diagramas de Cortante, Momento e Torque.
    Plota componentes Vertical, Horizontal e a Resultante.
    """
    # Cria diretório se não existir
    import os
    if not os.path.exists('Outputs'):
        os.makedirs('Outputs')

    x = np.linspace(0, L_total, 500)
    xA = posicoes_elementos['A']
    xB = posicoes_elementos['B']

    # 1. CÁLCULO DAS REAÇÕES DE APOIO (Equilíbrio Estático)
    # Cargas externas passadas na lista 'forcas'
    cargas = forcas

    # Somatório de Momentos em A para achar Rb
    M_A_v = sum(Fv * (pos - xA) for pos, Fv, Fh in cargas)
    M_A_h = sum(Fh * (pos - xA) for pos, Fv, Fh in cargas)

    R_Bv = M_A_v / (xB - xA)
    R_Bh = M_A_h / (xB - xA)

    # Somatório de Forças para achar Ra
    R_Av = sum(c[1] for c in cargas) - R_Bv
    R_Ah = sum(c[2] for c in cargas) - R_Bh

    # Lista completa de forças para o método das seções
    # Reações entram com sinal oposto às cargas externas para fechar o diagrama
    todas_forcas = cargas.copy()
    todas_forcas.append((xA, -R_Av, -R_Ah))
    todas_forcas.append((xB, -R_Bv, -R_Bh))
    todas_forcas.sort(key=lambda k: k[0])

    # Arrays para os diagramas
    V_v, V_h, M_v, M_h, T_plot = [np.zeros_like(x) for _ in range(5)]

    # 2. CÁLCULO DOS ESFORÇOS PONTO A PONTO
    for i, xi in enumerate(x):
        for pos, Fv, Fh in todas_forcas:
            if pos <= xi:
                # Lógica de Sinais para Diagrama:
                # Se for reação (A ou B), soma. Se for carga externa, subtrai.
                # (Assumindo convenção padrão de engenharia)
                is_reaction = (abs(pos - xA) < 0.1 or abs(pos - xB) < 0.1)

                # Ajuste de sinal: Reação sobe o diagrama, Carga desce
                val_v = abs(Fv) if is_reaction else -abs(Fv)
                val_h = abs(Fh) if is_reaction else -abs(Fh)

                # Cortante (Soma das forças à esquerda)
                V_v[i] += val_v
                V_h[i] += val_h

                # Momento (Força * Braço)
                M_v[i] += val_v * (xi - pos)
                M_h[i] += val_h * (xi - pos)

        # Lógica do Torque (Mantida)
        elementos_torque = [k for k in posicoes_elementos.keys() if 'Eng' in k or 'Coroa' in k or 'Pinhao' in k]
        if len(elementos_torque) >= 1: # Eixo simples
             if posicoes_elementos['A'] <= xi <= posicoes_elementos[elementos_torque[0]]:
                 T_plot[i] = torque

        if len(elementos_torque) >= 2: # Eixo intermediário
             p_t = sorted([posicoes_elementos[k] for k in elementos_torque])
             if p_t[0] <= xi <= p_t[-1]:
                 T_plot[i] = torque

    # Resultantes (Pitágoras)
    V_res = np.sqrt(V_v**2 + V_h**2)
    M_res = np.sqrt(M_v**2 + M_h**2)

    # 3. PLOTAGEM
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

    # --- Gráfico 1: Esforço Cortante (V) ---
    ax1.plot(x, V_v, label='Vertical', color='blue', alpha=0.5, linestyle='--')
    ax1.plot(x, V_h, label='Horizontal', color='green', alpha=0.5, linestyle='--')
    ax1.plot(x, V_res, label='Resultante', color='black', linewidth=2)
    ax1.set_ylabel('Cortante (N)')
    ax1.set_title(f'{nome_eixo} - Esforço Cortante')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)

    # --- Gráfico 2: Momento Fletor (M) ---
    ax2.plot(x, M_v, label='Vertical', color='blue', alpha=0.5, linestyle='--')
    ax2.plot(x, M_h, label='Horizontal', color='green', alpha=0.5, linestyle='--')
    ax2.plot(x, M_res, label='Resultante', color='red', linewidth=2)
    ax2.set_ylabel('Momento Fletor (N.mm)')
    ax2.set_title(f'{nome_eixo} - Momento Fletor')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)

    # Anotação do Máximo
    max_M = np.max(M_res)
    pos_max = x[np.argmax(M_res)]
    ax2.annotate(f'Mmax: {max_M:.0f}', xy=(pos_max, max_M),
                 xytext=(pos_max, max_M*1.1),
                 arrowprops=dict(facecolor='black', arrowstyle='->'))

    # --- Gráfico 3: Torque (T) ---
    ax3.plot(x, T_plot/1000, 'purple', lw=2, label='Torque')
    ax3.fill_between(x, T_plot/1000, color='purple', alpha=0.1)
    ax3.set_ylabel('Torque (N.m)')
    ax3.set_xlabel('Posição no Eixo (mm)')
    ax3.set_title(f'{nome_eixo} - Torque')
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()

    nome_arquivo = f'Plots/Diagrama_{nome_eixo.replace(" ", "_")}.png'
    plt.savefig(nome_arquivo, dpi=300)
    plt.close()
