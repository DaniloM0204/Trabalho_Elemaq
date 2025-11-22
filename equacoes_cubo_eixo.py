import numpy as np

# ============================================================
# PASSO 3.1 - DIMENSIONAR CHAVETAS PARA CADA ENGRENAGEM NOS EIXOS
# ============================================================

def dimensiona_chavetas(diametro_eixo, torque, largura_engrenagem, nome_chaveta):
    """
    Dimensiona chavetas para engrenagens nos eixos conforme norma ABNT
    """
    
    # Tabela de chavetas paralelas conforme ABNT
    # [diâmetro_min, diâmetro_max, largura_b, altura_h, profundidade_eixo_t1, profundidade_cubo_t2]
    tabela_chavetas = [
        [6, 8, 2, 2, 1.2, 1.0],
        [8, 10, 3, 3, 1.8, 1.4],
        [10, 12, 4, 4, 2.5, 1.8],
        [12, 17, 5, 5, 3.0, 2.3],
        [17, 22, 6, 6, 3.5, 2.8],
        [22, 30, 8, 7, 4.0, 3.3],
        [30, 38, 10, 8, 5.0, 3.3],
        [38, 44, 12, 8, 5.0, 3.3],
        [44, 50, 14, 9, 5.5, 3.8],
        [50, 58, 16, 10, 6.0, 4.3],
        [58, 65, 18, 11, 7.0, 4.4],
        [65, 75, 20, 12, 7.5, 4.9]
    ]
    
    # Selecionar chaveta baseada no diâmetro do eixo
    chaveta_selecionada = None
    for faixa in tabela_chavetas:
        if faixa[0] <= diametro_eixo <= faixa[1]:
            chaveta_selecionada = {
                'b': faixa[2],  # largura
                'h': faixa[3],  # altura
                't1': faixa[4], # profundidade no eixo
                't2': faixa[5]  # profundidade no cubo
            }
            break
    
    if chaveta_selecionada is None:
        raise ValueError(f"Diâmetro do eixo {diametro_eixo} mm fora da faixa da tabela de chavetas")
    
    # Comprimento da chaveta (usar 80% da largura da engrenagem como prática comum)
    comprimento_chaveta = 0.8 * largura_engrenagem
    
    # Força tangencial na chaveta
    raio_eixo = diametro_eixo / 2
    F_t = torque / raio_eixo  # N
    
    return {
        'nome': nome_chaveta,
        'diametro_eixo': diametro_eixo,
        'torque': torque,
        'dimensoes_chaveta': chaveta_selecionada,
        'comprimento': comprimento_chaveta,
        'forca_tangencial': F_t
    }

# ============================================================
# PASSO 3.2 - VERIFICAÇÃO AO CISALHAMENTO DAS CHAVETAS
# ============================================================

def verifica_cisalhamento_chaveta(resultado_chaveta, S_y, C_seg):
    """
    Verifica a chaveta ao cisalhamento
    """
    b = resultado_chaveta['dimensoes_chaveta']['b']
    L = resultado_chaveta['comprimento']
    F_t = resultado_chaveta['forca_tangencial']
    
    # Área resistente ao cisalhamento
    A_cisalhamento = b * L  # mm²
    
    # Tensão de cisalhamento
    tau = F_t / A_cisalhamento  # MPa
    
    # Tensão admissível ao cisalhamento (0.577 * Sy segundo critério de Von Mises)
    tau_adm = (0.577 * S_y) / C_seg
    
    # Fator de segurança ao cisalhamento
    FS_cisalhamento = tau_adm / tau if tau > 0 else float('inf')
    
    resultado_chaveta['tau_cisalhamento'] = tau
    resultado_chaveta['tau_admissivel'] = tau_adm
    resultado_chaveta['FS_cisalhamento'] = FS_cisalhamento
    resultado_chaveta['status_cisalhamento'] = "ATENDE" if FS_cisalhamento >= 1 else "NÃO ATENDE"
    
    return resultado_chaveta

# ============================================================
# PASSO 3.3 - VERIFICAÇÃO AO ESMAGAMENTO DAS CHAVETAS
# ============================================================

def verifica_esmagamento_chaveta(resultado_chaveta, S_y, C_seg):
    """
    Verifica a chaveta ao esmagamento (compressão)
    """
    h = resultado_chaveta['dimensoes_chaveta']['h']
    L = resultado_chaveta['comprimento']
    F_t = resultado_chaveta['forca_tangencial']
    
    # Área resistente ao esmagamento (metade da altura da chaveta)
    A_esmagamento = (h / 2) * L  # mm²
    
    # Tensão de esmagamento
    sigma_esmagamento = F_t / A_esmagamento  # MPa
    
    # Tensão admissível ao esmagamento
    sigma_adm_esmagamento = S_y / C_seg
    
    # Fator de segurança ao esmagamento
    FS_esmagamento = sigma_adm_esmagamento / sigma_esmagamento if sigma_esmagamento > 0 else float('inf')
    
    resultado_chaveta['sigma_esmagamento'] = sigma_esmagamento
    resultado_chaveta['sigma_adm_esmagamento'] = sigma_adm_esmagamento
    resultado_chaveta['FS_esmagamento'] = FS_esmagamento
    resultado_chaveta['status_esmagamento'] = "ATENDE" if FS_esmagamento >= 1 else "NÃO ATENDE"
    
    return resultado_chaveta
