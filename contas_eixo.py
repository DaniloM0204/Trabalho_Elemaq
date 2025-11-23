import equacoes_eixo as eixo
import utilities as util
import math

parametros = util.ler_Estagios_Engrenagem("Outputs/Estagios_Engrenagem.txt")

parametros_torque = {
    "Torque eixo1": 10.41,
    "Torque eixo2": 35.18,
    "Torque eixo3": 118.91,
}

# Material Aco AISI 1020
S_ut = 469
S_y = 390
Se_linha = 0.5 * S_ut
C_seg = 2.5

# Configuracao Geométrica
posicoes_eixo1 = {"A": 0, "pinhao1": 30, "B": 60}
L_e1 = 60
posicoes_eixo2 = {"A": 0, "coroa1": 30, "pinhao2": 100, "B": 130}
L_e2 = 130
posicoes_eixo3 = {"A": 0, "coroa2": 30, "B": 60}
L_e3 = 60


def check_force(val):
    return val * 1000 if val < 100 else val


Wt1 = check_force(parametros["parametros_eixos"]["forca_tangencial_pinhao1"])
Wr1 = check_force(parametros["parametros_eixos"]["forca_radial_pinhao1"])
Wtc1 = check_force(parametros["parametros_eixos"]["forca_tangencial_coroa1"])
Wrc1 = check_force(parametros["parametros_eixos"]["forca_radial_coroa1"])
Wtp2 = check_force(parametros["parametros_eixos"]["forca_tangencial_pinhao2"])
Wrp2 = check_force(parametros["parametros_eixos"]["forca_radial_pinhao2"])
Wtc2 = check_force(parametros["parametros_eixos"]["forca_tangencial_coroa2"])
Wrc2 = check_force(parametros["parametros_eixos"]["forca_radial_coroa2"])

# ============================================================
# CaLCULO DOS ESFORcOS USANDO NOVA FUNcaO
# ============================================================

# Eixo 1 - Entrada
forcas_eixo1 = [
    (posicoes_eixo1["pinhao1"], Wr1, Wt1)
]  # (pos, F_vertical, F_horizontal)
resultados_eixo1 = eixo.calcula_esforcos_analitico(
    L_total=L_e1,
    posicoes=posicoes_eixo1,
    forcas=forcas_eixo1,
    torque=parametros_torque["Torque eixo1"] * 1000,
    nome_eixo="Eixo 1 - Entrada",
)

# Eixo 2 - Intermediario
forcas_eixo2 = [
    (posicoes_eixo2["coroa1"], Wrc1, Wtc1),
    (posicoes_eixo2["pinhao2"], Wrp2, Wtp2),
]
resultados_eixo2 = eixo.calcula_esforcos_analitico(
    L_total=L_e2,
    posicoes=posicoes_eixo2,
    forcas=forcas_eixo2,
    torque=parametros_torque["Torque eixo2"] * 1000,
    nome_eixo="Eixo 2 - Intermediario",
)

# Eixo 3 - Saida
forcas_eixo3 = [(posicoes_eixo3["coroa2"], Wrc2, Wtc2)]
resultados_eixo3 = eixo.calcula_esforcos_analitico(
    L_total=L_e3,
    posicoes=posicoes_eixo3,
    forcas=forcas_eixo3,
    torque=parametros_torque["Torque eixo3"] * 1000,
    nome_eixo="Eixo 3 - Saida",
)

# Salvar resultados
resultados_diagramas = {
    "eixo1": resultados_eixo1,
    "eixo2": resultados_eixo2,
    "eixo3": resultados_eixo3,
}

with open("Outputs/resultados_diagramas.txt", "w") as f:
    f.write("RESULTADOS DOS DIAGRAMAS - NOVO MÉTODO\n")
    f.write("=" * 50 + "\n\n")

    for ex, dados in resultados_diagramas.items():
        f.write(f"{ex.upper()}:\n")
        f.write(f"  Momento maximo: {dados['M_max']:.0f} N.mm\n")
        f.write(f"  Posicao M_max: {dados['posicao_M_max']:.1f} mm\n")
        f.write(f"  Reacao A: {dados['R_A_resultante']:.0f} N\n")
        f.write(f"  Reacao B: {dados['R_B_resultante']:.0f} N\n")
        f.write(f"  Reacao Av: {dados['R_Av']:.0f} N\n")
        f.write(f"  Reacao Ah: {dados['R_Ah']:.0f} N\n")
        f.write(f"  Reacao Bv: {dados['R_Bv']:.0f} N\n")
        f.write(f"  Reacao Bh: {dados['R_Bh']:.0f} N\n\n")

# ============================================================
# PLOTAGEM DOS DIAGRAMAS COM NOVA FUNcaO
# ============================================================

# Eixo 1
eixo.plotar_diagramas_completo(
    L_total=L_e1,
    posicoes=posicoes_eixo1,
    forcas=forcas_eixo1,
    torque=parametros_torque["Torque eixo1"] * 1000,
    nome_eixo="Eixo 1 - Entrada",
    d_estimado=30,
)

# Eixo 2
eixo.plotar_diagramas_completo(
    L_total=L_e2,
    posicoes=posicoes_eixo2,
    forcas=forcas_eixo2,
    torque=parametros_torque["Torque eixo2"] * 1000,
    nome_eixo="Eixo 2 - Intermediario",
    d_estimado=35,
)

# Eixo 3
eixo.plotar_diagramas_completo(
    L_total=L_e3,
    posicoes=posicoes_eixo3,
    forcas=forcas_eixo3,
    torque=parametros_torque["Torque eixo3"] * 1000,
    nome_eixo="Eixo 3 - Saida",
    d_estimado=40,
)

# ============================================================
# DIMENSIONAMENTO POR FADIGA E ANaLISE DE DEFLEXaO
# ============================================================

# Eixo 1
resultado_e1 = eixo.dimensiona_eixo_por_fadiga(
    resultados_diagramas["eixo1"]["M_max"],
    parametros_torque["Torque eixo1"],
    S_ut,
    S_y,
    Se_linha,
    vida_util_horas=4000,
    Nf=C_seg,
    tipo_eixo="simples",
)
deflexao_e1 = eixo.analisa_deflexao_eixo(
    resultados_diagramas["eixo1"]["M_max"],
    L_e1,
    resultado_e1["diametro_minimo"],
)
velocidade_critica_e1 = eixo.calcula_velocidade_critica(
    resultado_e1["diametro_minimo"], L_e1
)

# Eixo 2
resultado_e2 = eixo.dimensiona_eixo_por_fadiga(
    resultados_diagramas["eixo2"]["M_max"],
    parametros_torque["Torque eixo2"],
    S_ut,
    S_y,
    Se_linha,
    vida_util_horas=4000,
    Nf=C_seg,
    tipo_eixo="duplo",
)
deflexao_e2 = eixo.analisa_deflexao_eixo(
    resultados_diagramas["eixo2"]["M_max"],
    L_e2,
    resultado_e2["diametro_minimo"],
)
velocidade_critica_e2 = eixo.calcula_velocidade_critica(
    resultado_e2["diametro_minimo"], L_e2
)

# Eixo 3
resultado_e3 = eixo.dimensiona_eixo_por_fadiga(
    resultados_diagramas["eixo3"]["M_max"],
    parametros_torque["Torque eixo3"],
    S_ut,
    S_y,
    Se_linha,
    vida_util_horas=4000,
    Nf=C_seg,
    tipo_eixo="simples",
)
deflexao_e3 = eixo.analisa_deflexao_eixo(
    resultados_diagramas["eixo3"]["M_max"],
    L_e3,
    resultado_e3["diametro_minimo"],
)
velocidade_critica_e3 = eixo.calcula_velocidade_critica(
    resultado_e3["diametro_minimo"], L_e3
)

# ============================================================
# RESUMO COMPLETO
# ============================================================

resumo_eixos = {
    "Eixo 1 (Entrada)": {
        **resultado_e1,
        **deflexao_e1,
        "velocidade_critica_rpm": velocidade_critica_e1[
            "velocidade_critica_rpm"
        ],
    },
    "Eixo 2 (Intermediario)": {
        **resultado_e2,
        **deflexao_e2,
        "velocidade_critica_rpm": velocidade_critica_e2[
            "velocidade_critica_rpm"
        ],
    },
    "Eixo 3 (Saida)": {
        **resultado_e3,
        **deflexao_e3,
        "velocidade_critica_rpm": velocidade_critica_e3[
            "velocidade_critica_rpm"
        ],
    },
}

# Diametros comerciais
d1_com = max(
    20, math.ceil(resumo_eixos["Eixo 1 (Entrada)"]["diametro_minimo"] / 5) * 5
)
d2_com = max(
    25,
    math.ceil(resumo_eixos["Eixo 2 (Intermediario)"]["diametro_minimo"] / 5)
    * 5,
)
d3_com = max(
    30, math.ceil(resumo_eixos["Eixo 3 (Saida)"]["diametro_minimo"] / 5) * 5
)

with open("Outputs/dimensionamento_eixos.txt", "w") as f:
    f.write("DIMENSIONAMENTO DOS EIXOS - RESULTADOS FINAIS\n")
    f.write("=" * 60 + "\n\n")
    f.write("Material: Aco AISI 1020\n")
    f.write(f"S_ut = {S_ut} MPa, S_y = {S_y} MPa\n")
    f.write(f"Fator de seguranca: {C_seg}\n")
    f.write(f"Vida util: 4000 horas\n\n")

    for nome, dados in resumo_eixos.items():
        f.write(f"{nome}:\n")
        f.write(f"  Diametro minimo: {dados['diametro_minimo']:.1f} mm\n")
        f.write(f"  FS Fadiga: {dados['FS_fadiga']:.2f}\n")
        f.write(f"  FS Escoamento: {dados['FS_escoamento']:.2f}\n")
        f.write(
            f"  Limite de fadiga corrigido: {dados['Se_corrigido']:.0f} MPa\n"
        )
        f.write(f"  Deflexao maxima: {dados['deflexao_max_mm']:.3f} mm\n")
        f.write(
            f"  Deflexao admissivel: {dados['deflexao_admissivel_mm']:.3f} mm\n"
        )
        f.write(f"  Status deflexao: {dados['status_deflexao']}\n")
        f.write(
            f"  Velocidade critica: {dados['velocidade_critica_rpm']:.0f} rpm\n"
        )
        f.write(f"  Iteracoes: {dados['iteracoes']}\n\n")

    f.write("RECOMENDACOES:\n")
    f.write(
        f"- Diametros comerciais sugeridos: {d1_com}mm, {d2_com}mm, {d3_com}mm\n"
    )
    f.write("- Implementar raios de concordancia nos degraus\n")
    f.write("- Verificar compatibilidade com rolamentos selecionados\n")

print("Analise completa concluida!")
print(
    f"Diametros minimos: E1={resultado_e1['diametro_minimo']:.1f}mm, E2={resultado_e2['diametro_minimo']:.1f}mm, E3={resultado_e3['diametro_minimo']:.1f}mm"
)
print("Arquivos salvos em: Outputs/")
