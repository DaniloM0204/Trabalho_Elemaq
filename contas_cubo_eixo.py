import equacoes_cubo_eixo as ce
import contas_eixo as eixo
import utilities as util

param = util.ler_Estagios_Engrenagem('Outputs/Estagios_Engrenagem.txt')



# VARIÁVEIS NECESSÁRIAS PARA O DIMENSIONAMENTO DE CHAVETAS
# (Valores extraídos das fases anteriores do código)

# Diâmetros dos eixos dimensionados
diametro_eixo1 = eixo.resultado_e1['diametro_minimo']
diametro_eixo2 = eixo.resultado_e2['diametro_minimo']
diametro_eixo3 = eixo.resultado_e3['diametro_minimo']

# Torques nos eixos
torque_eixo1 = param['parametros_eixos']['Torque eixo1']
torque_eixo2 = param['parametros_eixos']['Torque eixo2']
torque_eixo3 = param['parametros_eixos']['Torque eixo3']

# Propriedades do material (Aço AISI 1020)
S_y = 390  # MPa - Tensão de escoamento
S_ut = 469  # MPa - Tensão de ruptura

# Fator de segurança
C_seg = 1.5

estagio_1 = param['estagio_1']
estagio_2 = param['estagio_2']


# Larguras das engrenagens (assumidas - seriam calculadas no PASSO 1.2)
largura_engrenagem_eixo1 = estagio_1['b1']  # mm
largura_engrenagem_eixo2_coroa1 =estagio_1['b2']  # mm
largura_engrenagem_eixo2_pinhao2 = estagio_2['b1']  # mm
largura_engrenagem_eixo3 = estagio_2['b2']  # mm

# ============================================================
# APLICAÇÃO PARA TODAS AS ENGRENAGENS
# ============================================================

# Dimensionar chavetas para cada engrenagem
chaveta_eixo1_pinhao1 = ce.dimensiona_chavetas(
    diametro_eixo1, torque_eixo1, largura_engrenagem_eixo1, "Chaveta Eixo1-Pinhão1"
)

chaveta_eixo2_coroa1 = ce.dimensiona_chavetas(
    diametro_eixo2, torque_eixo2, largura_engrenagem_eixo2_coroa1, "Chaveta Eixo2-Coroa1"
)

chaveta_eixo2_pinhao2 = ce.dimensiona_chavetas(
    diametro_eixo2, torque_eixo2, largura_engrenagem_eixo2_pinhao2, "Chaveta Eixo2-Pinhão2"
)

chaveta_eixo3_coroa2 = ce.dimensiona_chavetas(
    diametro_eixo3, torque_eixo3, largura_engrenagem_eixo3, "Chaveta Eixo3-Coroa2"
)

# Verificar todas as chavetas ao cisalhamento e esmagamento
chavetas = [
    ce.verifica_esmagamento_chaveta(ce.verifica_cisalhamento_chaveta(chaveta_eixo1_pinhao1)),
    ce.verifica_esmagamento_chaveta(ce.verifica_cisalhamento_chaveta(chaveta_eixo2_coroa1)),
    ce.verifica_esmagamento_chaveta(ce.verifica_cisalhamento_chaveta(chaveta_eixo2_pinhao2)),
    ce.verifica_esmagamento_chaveta(ce.verifica_cisalhamento_chaveta(chaveta_eixo3_coroa2))
]

# ============================================================
# RELATÓRIO DAS CHAVETAS
# ============================================================

with open("dimensionamento_chavetas.txt", "w") as f:
    f.write("DIMENSIONAMENTO DE CHAVETAS - RESULTADOS\n")
    f.write("="*60 + "\n\n")

    f.write("Material: Aço AISI 1020\n")
    f.write(f"S_y = {S_y} MPa, S_ut = {S_ut} MPa\n")
    f.write(f"Fator de segurança: {C_seg}\n\n")

    for chaveta in chavetas:
        f.write(f"{chaveta['nome']}:\n")
        f.write(f"  Diâmetro do eixo: {chaveta['diametro_eixo']:.1f} mm\n")
        f.write(f"  Torque transmitido: {chaveta['torque']:.0f} N.mm\n")
        f.write(f"  Dimensões chaveta: {chaveta['dimensoes_chaveta']['b']} x {chaveta['dimensoes_chaveta']['h']} mm\n")
        f.write(f"  Comprimento chaveta: {chaveta['comprimento']:.1f} mm\n")
        f.write(f"  Força tangencial: {chaveta['forca_tangencial']:.0f} N\n")
        f.write(f"  Tensão cisalhamento: {chaveta['tau_cisalhamento']:.1f} MPa\n")
        f.write(f"  Tensão admissível cisalhamento: {chaveta['tau_admissivel']:.1f} MPa\n")
        f.write(f"  FS cisalhamento: {chaveta['FS_cisalhamento']:.2f} ({chaveta['status_cisalhamento']})\n")
        f.write(f"  Tensão esmagamento: {chaveta['sigma_esmagamento']:.1f} MPa\n")
        f.write(f"  Tensão admissível esmagamento: {chaveta['sigma_adm_esmagamento']:.1f} MPa\n")
        f.write(f"  FS esmagamento: {chaveta['FS_esmagamento']:.2f} ({chaveta['status_esmagamento']})\n\n")

    # Verificação geral
    todas_atendem = all(chaveta['status_cisalhamento'] == "ATENDE" and
                       chaveta['status_esmagamento'] == "ATENDE" for chaveta in chavetas)

    f.write("VERIFICAÇÃO GERAL:\n")
    if todas_atendem:
        f.write("✓ TODAS AS CHAVETAS ATENDEM AOS CRITÉRIOS DE DIMENSIONAMENTO\n")
    else:
        f.write("✗ ALGUMAS CHAVETAS NÃO ATENDEM - REVISAR DIMENSIONAMENTO\n")

    f.write("\nRECOMENDAÇÕES:\n")
    f.write("- Utilizar chavetas paralelas ABNT\n")
    f.write("- Ajustar comprimento da chaveta conforme largura real da engrenagem\n")
    f.write("- Verificar folgas de montagem entre chaveta e rasgo\n")
