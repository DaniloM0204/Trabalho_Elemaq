import equacoes_eixo as eixo
import utilities as util
import math

parametros = util.ler_Estagios_Engrenagem('Outputs/Estagios_Engrenagem.txt')

parametros_torque = {
    'Torque eixo1': 10.41,
    'Torque eixo2': 35.18,
    'Torque eixo3': 118.91
}

# Material Aco AISI 1020
S_ut = 469
S_y = 390
Se_linha = 0.5 * S_ut
C_seg = 2.5

# Configuração Geométrica
distancias_eixo1 = {'LA': 50, 'AB': 200, 'BC': 150}
L_e1 = 200
distancias_eixo2 = {'LA': 40, 'L1': 120, 'L2': 40} # LA=40, L1=120, L2=40 -> Total 200
L_e2 = 200
distancias_eixo3 = {'LA': 60, 'AB': 200, 'BC': 160}
L_e3 = 200

def check_force(val):
    return val * 1000 if val < 100 else val

Wt1 = check_force(parametros['parametros_eixos']['forca_tangencial_pinhao1'])
Wr1 = check_force(parametros['parametros_eixos']['forca_radial_pinhao1'])
Wtc1 = check_force(parametros['parametros_eixos']['forca_tangencial_coroa1'])
Wrc1 = check_force(parametros['parametros_eixos']['forca_radial_coroa1'])
Wtp2 = check_force(parametros['parametros_eixos']['forca_tangencial_pinhao2'])
Wrp2 = check_force(parametros['parametros_eixos']['forca_radial_pinhao2'])
Wtc2 = check_force(parametros['parametros_eixos']['forca_tangencial_coroa2'])
Wrc2 = check_force(parametros['parametros_eixos']['forca_radial_coroa2'])

# Cálculo das reações nos mancais
reações_eixo1 = eixo.calcula_reações_mancais_eixo1({'W_t': Wt1, 'W_r': Wr1}, distancias_eixo1)
reações_eixo2 = eixo.calcula_reações_eixo2({'W_t': Wtc1, 'W_r': Wrc1}, {'W_t': Wtp2, 'W_r': Wrp2}, distancias_eixo2)
reações_eixo3 = eixo.calcula_reações_mancais_eixo3({'W_t': Wtc2, 'W_r': Wrc2}, distancias_eixo3)

res_eixo1 = eixo.calcula_diagramas_eixo_simples({'W_t': Wt1, 'W_r': Wr1}, distancias_eixo1, parametros_torque['Torque eixo1']*1000, "Eixo 1")
res_eixo2 = eixo.calcula_diagramas_eixo_duplo({'W_t': Wtc1, 'W_r': Wrc1}, {'W_t': Wtp2, 'W_r': Wrp2}, distancias_eixo2, parametros_torque['Torque eixo2']*1000, "Eixo 2")
res_eixo3 = eixo.calcula_diagramas_eixo_simples({'W_t': Wtc2, 'W_r': Wrc2}, distancias_eixo3, parametros_torque['Torque eixo3']*1000, "Eixo 3")

resultados_diagramas = {'eixo1': res_eixo1, 'eixo2': res_eixo2, 'eixo3': res_eixo3}

# Eixo 1: A(0) -> Eng(50) -> B(200)
eixo.plotar_diagramas_eixo(
    L_total=L_e1,
    posicoes_elementos={'A': 0, 'Eng1': 50, 'B': 200},
    forcas=[(50, Wr1, Wt1)],
    torque=parametros_torque['Torque eixo1']*1000,
    nome_eixo="Eixo 1 - Entrada"
)

# Eixo 2: A(0) -> Coroa1(40) -> Pinhao2(160) -> B(200)
# Nota: Pinhao2 fica em 40 + 120 = 160
eixo.plotar_diagramas_eixo(
    L_total=L_e2,
    posicoes_elementos={'A': 0, 'Coroa1': 40, 'Pinhao2': 160, 'B': 200},
    forcas=[(40, Wrc1, Wtc1), (160, Wrp2, Wtp2)],
    torque=parametros_torque['Torque eixo2']*1000,
    nome_eixo="Eixo 2 - Intermediario"
)

# Eixo 3: A(0) -> Coroa2(60) -> B(200)
eixo.plotar_diagramas_eixo(
    L_total=L_e3,
    posicoes_elementos={'A': 0, 'Coroa2': 60, 'B': 200},
    forcas=[(60, Wrc2, Wtc2)],
    torque=parametros_torque['Torque eixo3']*1000,
    nome_eixo="Eixo 3 - Saida"
)

with open("Outputs/resultados_diagramas.txt", "w") as f:
    for ex, dados in resultados_diagramas.items():
        f.write(f"{ex.upper()}:\n")
        f.write(f"  Momento maximo: {dados['M_max']:.0f} N.mm\n")
        f.write(f"  Reacao A: {dados['R_A_resultante']:.0f} N\n")
        f.write(f"  Reacao B: {dados['R_B_resultante']:.0f} N\n\n")

# Dimensionamento por fadiga e analise de deflexao
resultado_e1 = eixo.dimensiona_eixo_por_fadiga(resultados_diagramas['eixo1']['M_max'], parametros_torque['Torque eixo1'], S_ut, S_y, Se_linha, Nf=C_seg, tipo_eixo="simples")
deflexao_e1 = eixo.analisa_deflexao_eixo(resultados_diagramas['eixo1']['M_max'], L_e1, resultado_e1['diametro_minimo'])

resultado_e2 = eixo.dimensiona_eixo_por_fadiga(resultados_diagramas['eixo2']['M_max'], parametros_torque['Torque eixo2'], S_ut, S_y, Se_linha, Nf=C_seg, tipo_eixo="duplo")
deflexao_e2 = eixo.analisa_deflexao_eixo(resultados_diagramas['eixo2']['M_max'], L_e2, resultado_e2['diametro_minimo'])

resultado_e3 = eixo.dimensiona_eixo_por_fadiga(resultados_diagramas['eixo3']['M_max'], parametros_torque['Torque eixo3'], S_ut, S_y, Se_linha, Nf=C_seg, tipo_eixo="simples")
deflexao_e3 = eixo.analisa_deflexao_eixo(resultados_diagramas['eixo3']['M_max'], L_e3, resultado_e3['diametro_minimo'])

resumo_eixos = {
    "Eixo 1": {**resultado_e1, **deflexao_e1},
    "Eixo 2": {**resultado_e2, **deflexao_e2},
    "Eixo 3": {**resultado_e3, **deflexao_e3}
}

d1_com = max(20, math.ceil(resumo_eixos["Eixo 1"]['diametro_minimo']/5)*5)
d2_com = max(25, math.ceil(resumo_eixos["Eixo 2"]['diametro_minimo']/5)*5)
d3_com = max(30, math.ceil(resumo_eixos["Eixo 3"]['diametro_minimo']/5)*5)

with open("Outputs/dimensionamento_eixos.txt", "w") as f:
    f.write("="*60 + "\n\n")
    f.write("Material: Aco AISI 1020\n")
    f.write(f"Fator de seguranca: {C_seg}\n\n")

    for nome, d in resumo_eixos.items():
        f.write(f"{nome}:\n")
        f.write(f"  Diametro minimo: {d['diametro_minimo']:.2f} mm\n")
        f.write(f"  FS Fadiga: {d['FS_fadiga']:.2f}\n")
        f.write(f"  Deflexao maxima: {d['deflexao_max_mm']:.3f} mm\n\n")

    f.write("RECOMENDACOES:\n")
    f.write(f"- Diametros comerciais: {d1_com}mm, {d2_com}mm, {d3_com}mm\n")
