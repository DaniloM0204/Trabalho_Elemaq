def dimensiona_chavetas(diametro_eixo, torque_Nmm, largura_engrenagem, nome_chaveta):
    """
    Dimensiona chavetas conforme norma ABNT/DIN 6885
    """

    # Tabela simplificada DIN 6885 (Aco)
    # [d_min, d_max, b, h]
    tabela_chavetas = [
        [6, 8, 2, 2],
        [8, 10, 3, 3],
        [10, 12, 4, 4],
        [12, 17, 5, 5],
        [17, 22, 6, 6],
        [22, 30, 8, 7],
        [30, 38, 10, 8],
        [38, 44, 12, 8],
        [44, 50, 14, 9],
        [50, 58, 16, 10],
        [58, 65, 18, 11],
        [65, 75, 20, 12]
    ]

    # Selecionar chaveta baseada no diametro do eixo
    sel = None
    for faixa in tabela_chavetas:
        # >= min e < max (ou <= max para o ultimo)
        if diametro_eixo >= faixa[0] and diametro_eixo <= faixa[1]:
            sel = {'b': faixa[2], 'h': faixa[3]}
            break

    # Se diametro for muito grande (usa regra pratica b=d/4)
    if sel is None:
        if diametro_eixo > 75:
            sel = {'b': int(diametro_eixo/4), 'h': int(diametro_eixo/4)}
        else:
            # Diametro muito pequeno, usa a menor
            sel = {'b': 2, 'h': 2}

    comprimento_chaveta = min(1.5 * diametro_eixo, largura_engrenagem)

    # Forca tangencial na chaveta (F = 2*T / d)
    F_t = (2 * torque_Nmm) / diametro_eixo

    return {
        'nome': nome_chaveta,
        'diametro_eixo': diametro_eixo,
        'torque': torque_Nmm,
        'dimensoes_chaveta': sel,
        'comprimento': comprimento_chaveta,
        'forca_tangencial': F_t
    }

def verifica_cisalhamento_chaveta(resultado_chaveta):
    b = resultado_chaveta['dimensoes_chaveta']['b']
    L = resultado_chaveta['comprimento']
    F_t = resultado_chaveta['forca_tangencial']

    # Material Aco AISI 1020
    S_y = 390  # MPa
    C_seg = 1.5

    # Area cisalhamento = b * L
    A_cis = b * L
    tau = F_t / A_cis

    # Von Mises: Sys = 0.577 * Sy
    tau_adm = (0.577 * S_y) / C_seg

    FS_cis = tau_adm / tau if tau > 0 else 999

    resultado_chaveta['tau_cisalhamento'] = tau
    resultado_chaveta['tau_admissivel'] = tau_adm
    resultado_chaveta['FS_cisalhamento'] = FS_cis
    resultado_chaveta['status_cisalhamento'] = "ATENDE" if FS_cis >= 1 else "FALHA"

    return resultado_chaveta


def verifica_esmagamento_chaveta(resultado_chaveta):
    h = resultado_chaveta['dimensoes_chaveta']['h']
    L = resultado_chaveta['comprimento']
    F_t = resultado_chaveta['forca_tangencial']

    S_y = 390
    C_seg = 1.5

    # Area esmagamento = (h/2) * L
    A_esm = (h / 2) * L
    sigma_esm = F_t / A_esm

    # Tensão admissível compressão = Sy (ou 0.9*Sy)
    sigma_adm = S_y / C_seg

    FS_esm = sigma_adm / sigma_esm if sigma_esm > 0 else 999

    resultado_chaveta['sigma_esmagamento'] = sigma_esm
    resultado_chaveta['sigma_adm_esmagamento'] = sigma_adm
    resultado_chaveta['FS_esmagamento'] = FS_esm
    resultado_chaveta['status_esmagamento'] = "ATENDE" if FS_esm >= 1 else "FALHA"

    return resultado_chaveta
