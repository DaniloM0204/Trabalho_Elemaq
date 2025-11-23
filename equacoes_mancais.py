import numpy as np

def calcula_carga_equivalente_rolamento(Frad, Fax):
    if Fax <= 0.01: return Frad
    X, Y = 0.56, 1.5
    P = X * Frad + Y * Fax
    return max(P, Frad)

def seleciona_rolamento(Frad, Fax, rpm, vida_util_horas, Kapl=1.0):
    P = calcula_carga_equivalente_rolamento(Frad, Fax)
    P_proj = P * Kapl
    L10 = (60 * rpm * vida_util_horas) / 1e6
    a = 3.0

    if L10 > 0: C_req = P_proj * (L10 ** (1/a))
    else: C_req = P_proj

    return {
        'C_req_N': C_req,
        'C_req_kN': C_req / 1000.0,
        'P_proj_N': P_proj
    }

def busca_catalogo(C_req_kN, d_min_mm, catalogo):
    # Filtra e ordena
    candidatos = [r for r in catalogo if r['d'] >= d_min_mm]
    if not candidatos: return None
    candidatos.sort(key=lambda x: (x['d'], x['C_din']))

    for rol in candidatos:
        if rol['C_din'] >= C_req_kN:
            return rol
    return None
