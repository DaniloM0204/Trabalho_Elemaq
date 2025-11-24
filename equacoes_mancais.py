import math


def calcula_carga_equivalente_rolamento(Frad, Fax):
    """
    Cálculo simplificado de P.
    Para engrenagens retas, Fax costuma ser 0, logo P = Frad.
    """
    if Fax <= 0.01:
        return float(Frad)

    # Fatores X e Y padrão para rolamentos rígidos de esferas (aprox.)
    X, Y = 0.56, 1.5
    P = X * Frad + Y * Fax
    return max(P, Frad)


def seleciona_rolamento(Frad, Fax, rpm, vida_util_horas, Kapl=1.0):
    """
    Calcula a Capacidade de Carga Dinâmica Requerida (C_req).
    """
    P = calcula_carga_equivalente_rolamento(Frad, Fax)
    P_proj = P * Kapl

    # L10 em milhões de ciclos
    L10_milhoes = (60 * rpm * vida_util_horas) / 10**6

    # Expoente: a=3 para esferas, a=10/3 para rolos
    a = 3.0

    if L10_milhoes > 0:
        C_req = P_proj * (L10_milhoes ** (1 / a))
    else:
        C_req = P_proj

    return {
        "C_req_N": C_req,
        "C_req_kN": C_req / 1000.0,
        "P_proj_N": P_proj,
        "L10_milhoes": L10_milhoes,
    }


def busca_catalogo(C_req_kN, d_eixo_mm, catalogo):
    """
    Busca o rolamento mais econômico que atenda C_req e tenha o diâmetro d_eixo.
    """
    # Filtra apenas rolamentos com o diâmetro exato do eixo (ou um pouco maior se necessário)
    # Rolamentos tem furo padrão (17, 20, 25, 30, 35, 40...)
    candidatos = [r for r in catalogo if r["d"] == d_eixo_mm]

    if not candidatos:
        # Se não achou exato, tenta achar um com d maior (necessitaria bucha ou degrau)
        # Mas para este projeto, vamos forçar o diametro exato ou retornar None
        return None

    # Ordena por capacidade de carga (C_din)
    candidatos.sort(key=lambda x: x["C_din"])

    for rol in candidatos:
        if rol["C_din"] >= C_req_kN:
            return rol

    return None  # Nenhum rolamento desse diâmetro aguenta a carga
