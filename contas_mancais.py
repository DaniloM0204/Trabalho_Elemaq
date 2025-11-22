import equacoes_mancais as mancal
import contas_eixo as eixo
import contas_cubo_eixo as ce
import utilities as util
import csv


dados_engre = util.ler_Estagios_Engrenagem('Outputs/Estagios_Engrenagem.txt')
dados_entradas = util.ler_Dados_De_Entrada('Inputs/Dados_De_Entrada.txt')
dados_eixos = util.ler_dimensionamento_eixos('Outputs/dimensionamento_eixos.txt')


vida_util = 4000  # horas


catalogo = []
with open('Inputs/rolamento_sfk.csv', newline='') as csvfile:
    leitor = csv.DictReader(csvfile)
    for linha in leitor:
        catalogo.append({
            'codigo': linha['codigo'],
            'd': float(linha['d']),
            'C_din': float(linha['C_din']),
        })

# Distancia e distancias entre engrenagens que tem que ajustar
L_total = (eixo.L_e1 + eixo.L_e2 + eixo.L_e3)/3
dist_e1 = eixo.distancias_eixo1
dist_e2 = eixo.distancias_eixo2
dist_e3 = eixo.distancias_eixo3

# Diametros minimos
diametro_minimo1 = ce.diametro_eixo1
diametro_minimo2 = ce.diametro_eixo2
diametro_minimo3 = ce.diametro_eixo3

# Reações nos mancais
ReacoesEixo1 = eixo.reações_eixo1
ReacoesEixo2 = eixo.reações_eixo2
ReacoesEixo3 = eixo.reações_eixo3

# Lista para os mancais
# Cada tupla: (nome, ReacaoA, ReacaoB, pos_rel_A, largura, diametro_minimo, vida_util, rotacao_rpm, catalogo)
mancais = [
    ('Mancal 1', ReacoesEixo1['R_A_resultante'], ReacoesEixo1['R_B_resultante'], dist_e1['LA'], dist_e1['AB'], diametro_minimo1, vida_util, dados_eixos.get('n_eixo1', 0), catalogo),
    ('Mancal 2', ReacoesEixo1['R_B_resultante'], 0.0, dist_e1['AB'], L_total - dist_e1['LA'] - dist_e1['AB'], diametro_minimo1, vida_util, dados_eixos.get('n_eixo1', 0), catalogo),
    ('Mancal 3', ReacoesEixo2['R_A_resultante'], ReacoesEixo2['R_B_resultante'], dist_e2['LA'], dist_e2['L1'], diametro_minimo2, vida_util, dados_eixos.get('n_eixo2', 0), catalogo),
    ('Mancal 4', ReacoesEixo2['R_B_resultante'], 0.0, dist_e2['L1'], L_total - dist_e2['LA'] - dist_e2['L1'], diametro_minimo2, vida_util, dados_eixos.get('n_eixo2', 0), catalogo),
    ('Mancal 5', ReacoesEixo3['R_A_resultante'], ReacoesEixo3['R_B_resultante'], dist_e3['LA'], dist_e3['AB'], diametro_minimo3, vida_util, dados_eixos.get('n_eixo3', 0), catalogo),
    ('Mancal 6', ReacoesEixo3['R_B_resultante'], 0.0, dist_e3['AB'], L_total - dist_e3['LA'] - dist_e3['AB'], diametro_minimo3, vida_util, dados_eixos.get('n_eixo3', 0), catalogo),
]

# Seleção dos mancais e registo dos resultados
with open('Outputs/Resultados_Mancais.txt', 'w') as arquivo_resultados:
    # Calculo da carga equivalente
    for (nome, R_A, R_B, posA, largura, diametro_minimo, vida_item, omega, catalogo_item) in mancais:
        # Carga radial no mancal: usar resultante máxima entre reações locais
        Frad = max(float(R_A or 0.0), float(R_B or 0.0))

        # Chama função compatível: equacoes_mancais.calcula_carga_equivalente_rolamento(Frad, Fax)
        P = mancal.calcula_carga_equivalente_rolamento(Frad, 0.0)

        # Capacidade dinâmica requerida (funcionalidade pode estar em outro módulo)
        try:
            seleciona = mancal.calcula_capacidade_dinamica_requerida(P=P, vida_util_horas=vida_item, omega=omega)
            C_req = seleciona.get('C_reqKn', None)
        except AttributeError:
            # Função não encontrada no módulo mancal — calcular C_req simples conhecendo P
            C_req = P / 1000.0

        # Busca no catalogo
        rolamento = mancal.busca_catalogo(Carga_requerida_kN=C_req, diametro_eixo_min=diametro_minimo, catalogo=catalogo)

        arquivo_resultados.write(f"{nome}:\n")
        arquivo_resultados.write(f"  Carga Radial (N): {Frad:.2f}\n")
        arquivo_resultados.write(f"  Velocidade (rpm): {omega:.2f}\n")
        arquivo_resultados.write(f"  Diametro minimo do eixo (mm): {diametro_minimo:.2f}\n")
        arquivo_resultados.write(f"  Carga Equivalente (N): {P:.2f}\n")

        if rolamento:
            msg = f'Escolhido o rolamento SFK {rolamento["codigo"]} com capacidade de carga dinamica de {rolamento["C_din"]} kN\n'
        else:
            msg = "Nenhum rolamento adequado encontrado no catálogo.\n"

        arquivo_resultados.write(msg)
        arquivo_resultados.write("\n")
