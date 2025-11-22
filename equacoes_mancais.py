def seleciona_rolamento(Frad, Fax, omega, vida_util, Kapl):
    """
    Para selecionar o rolamento vai ser lido um arquivo .csv com os dados dos rolamentos de cargas e etc
    """
    Carga_base = calcula_carga_equivalente_rolamento(Frad, Fax)
    Carga_projeto = Carga_base * Kapl
    # Para milhoes de revolucoes
    P_proj = (60*omega*vida_util)/1e6

    # Para maior confiabilidade, maior vida util portanto
    Kr = 1.0

    P_equivalente = P_proj * Kr

    # Calculo da Carga total
    a = 3 # Fator para rolamentos de esferas

    Carga_requerida = Carga_projeto * (P_equivalente)**(1/a)

    return {
        'Carga_req': Carga_requerida,
        'Carga_reqKn': Carga_requerida/1000,
        'Carga_Constante': P_proj,
        'Carga_projetoN': Carga_projeto
    }


def calcula_vida_util_rolamento(carga_dinamica, carga_equivalente, velocidade_rpm):
    """
    Calcula vida útil do rolamento (norma ISO)
    """
    pass

def calcula_carga_equivalente_rolamento(Frad, Fax, fatores_geometricos=None):
    """
    Pelo fato de termos carga puramente farial, a axial é negativa portanto a carga equivalente é igual a carga radial
    """
    if Fax == 0:
        return Frad

    # Caso haja carga axial para fatores padrão de esfera
    X = 0.56
    Y = 1.5
    P = X * Frad + Y * Fax

    return max(P, Frad)


def busca_catalogo(Carga_requerida_kN,diametro_eixo_min,catalogo):
    """
    Busca no catálogo o rolamento adequado pela capacidade de carga
    e diametro minimo do eixo
    Em que C_catalogo => C_requerida, Diametro_interno >= diametro_eixo_min
    """

    rolamento = [r for r in catalogo if r['d']>=diametro_eixo_min]

    # Organiza do primeiro diametro, depois pela capacidade de carga
    rolamento = sorted(rolamento, key=lambda x: (x['d'], x['C_din']))

    for rol in rolamento:
        if rol['C_din'] >= Carga_requerida_kN:
            return rol

    return None


def verifica_capacidade_carga_estatica_rolamento(carga_radial, carga_axial, capacidade_estatica):
    """
    Verifica capacidade de carga estática
    """
    pass

def calcula_reacoes_mancal(forcas_verticals, forcas_horizontais, distancias):
    """
    Calcula reações nos mancais considerando múltiplas forças
    """
    pass
