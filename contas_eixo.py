import equacoes_eixo as eixo
import utilities as util
import numpy as np

parametros = util.ler_Estagios_Engrenagem('Estagios_Engrenagem.txt')

# Aço AISI 1020 Normalizado
S_ut = 469  # MPa - Tensão de ruptura
S_y = 390   # MPa - Tensão de escoamento
Se_linha = 0.5 * S_ut  # Limite de fadiga não corrigido

# Coeficientes
C_seg = 1.5  # Fator de segurança
C_superf, acabamento_superf = eixo.calcula_c_superf(S_ut) #fator de acabamento superficial e o acabamento superficial definido
C_conf = 0.814  # Fator de confiabilidade de 99%
C_temp = 1.0   # Fator de temperatura para T < 450C
C_carr = 1.0   # Fator de carregamento para carga moderada

# Defina distâncias estimadas (em mm)
distancias_eixo1 = {'LA': 50, 'AB': 200, 'BC': 150}
L_e1 = 200
distancias_eixo2 = {'LA': 40, 'L1': 120, 'L2': 40}  # AB = 220mm
L_e2 = 220
distancias_eixo3 = {'LA': 60, 'AB': 200, 'BC': 160}
L_e3 = 200

# Calcule as reações
reações_eixo1 = eixo.calcula_reações_mancais_eixo1(
    {'W_t': parametros['parametros_eixos']['forca_tangencial_pinhao1'], 'W_r': parametros['parametros_eixos']['forca_radial_pinhao1']},  # Pinhão 1
    distancias_eixo1
)

reações_eixo2 = eixo.calcula_reações_eixo2(
    {'W_t': parametros['parametros_eixos']['forca_tangencial_coroa1'], 'W_r': parametros['parametros_eixos']['forca_radial_coroa1']},  # Coroa 1
    {'W_t': parametros['parametros_eixos']['forca_tangencial_pinhao2'], 'W_r': parametros['parametros_eixos']['forca_radial_pinhao2']},  # Pinhão 2  
    distancias_eixo2
)

reações_eixo3 = eixo.calcula_reações_mancais_eixo1(
    {'W_t': parametros['parametros_eixos']['forca_tangencial_coroa2'], 'W_r': parametros['parametros_eixos']['forca_radial_coroa2']},  # Coroa 2
    distancias_eixo3
)
# EIXO 1 - ENTRADA
resultados_eixo1 = eixo.calcula_diagramas_eixo_simples(
    forcas_engrenagem={'W_t': parametros['parametros_eixos']['forca_tangencial_pinhao1'], 'W_r': parametros['parametros_eixos']['forca_radial_pinhao1']},
    distancias=distancias_eixo1,
    torque=parametros['parametros_eixos']['torque_eixo1'],
    nome_eixo="Eixo 1 - Entrada"
)

# EIXO 2 - INTERMEDIÁRIO  
resultados_eixo2 = eixo.calcula_diagramas_eixo_duplo(
    forcas_engrenagem1={'W_t': parametros['parametros_eixos']['forca_tangencial_coroa1'], 'W_r': parametros['parametros_eixos']['forca_radial_coroa1']},  # Coroa 1
    forcas_engrenagem2={'W_t': parametros['parametros_eixos']['forca_tangencial_pinhao2'], 'W_r': parametros['parametros_eixos']['forca_radial_pinhao2']},  # Pinhão 2
    distancias=distancias_eixo2,
    torque=parametros['parametros_eixos']['torque_eixo2'],
    nome_eixo="Eixo 2 - Intermediário"
)

# EIXO 3 - SAÍDA
resultados_eixo3 = eixo.calcula_diagramas_eixo_simples(
    forcas_engrenagem={'W_t': parametros['parametros_eixos']['forca_tangencial_coroa2'], 'W_r': parametros['parametros_eixos']['forca_radial_coroa2']},
    distancias=distancias_eixo3, 
    torque=parametros['parametros_eixos']['torque_eixo3'],
    nome_eixo="Eixo 3 - Saída"
)

# Salvar resultados para uso no dimensionamento
resultados_diagramas = {
    'eixo1': resultados_eixo1,
    'eixo2': resultados_eixo2, 
    'eixo3': resultados_eixo3
}

# Salvar em arquivo
with open("resultados_diagramas.txt", "w") as f:
    f.write("RESULTADOS DOS DIAGRAMAS - VALORES CRÍTICOS\n")
    f.write("="*50 + "\n\n")
    
    for eixo_variavel, dados in resultados_diagramas.items():
        f.write(f"{eixo_variavel.upper()}:\n")
        f.write(f"  Momento máximo: {dados['M_max']:.0f} N.mm\n")
        f.write(f"  Cortante máximo: {dados['V_max']:.0f} N\n")
        f.write(f"  Reação A: {dados['R_A_resultante']:.0f} N\n")
        f.write(f"  Reação B: {dados['R_B_resultante']:.0f} N\n")
        f.write(f"  Posição M_max: {dados['posicao_M_max']:.1f} mm\n\n")

# Dimensionamento do EIXO 1 usando a função definida localmente
resultado_e1 = eixo.dimensiona_eixo_por_fadiga(
    resultados_diagramas['eixo1']['M_max'],
    parametros['parametros_eixos']['torque_eixo1'], 
    S_ut, S_y, Se_linha, 
    tipo_eixo="simples"
)

deflexao_e1 = eixo.analisa_deflexao_eixo(
    resultados_diagramas['eixo1']['M_max'],
    L_e1,
    resultado_e1['diametro_minimo']
)

velocidade_critica_e1 = eixo.calcula_velocidade_critica(
    resultado_e1['diametro_minimo'], 
    L_e1
)

# Dimensionamento do EIXO 2 
resultado_e2 = eixo.dimensiona_eixo_por_fadiga(
    resultados_diagramas['eixo2']['M_max'],
    parametros['parametros_eixos']['torque_eixo2'],
    S_ut, S_y, Se_linha,
    tipo_eixo="duplo"
)

deflexao_e2 = eixo.analisa_deflexao_eixo(
    resultados_diagramas['eixo2']['M_max'],
    L_e2,
    resultado_e2['diametro_minimo']
)

velocidade_critica_e2 = eixo.calcula_velocidade_critica(
    resultado_e2['diametro_minimo'], 
    L_e2
)

# Dimensionamento do EIXO 3
resultado_e3 = eixo.dimensiona_eixo_por_fadiga(
    resultados_diagramas['eixo3']['M_max'],
    parametros['parametros_eixos']['torque_eixo3'],
    S_ut, S_y, Se_linha,
    tipo_eixo="simples"
)

deflexao_e3 = eixo.analisa_deflexao_eixo(
    resultados_diagramas['eixo3']['M_max'],
    L_e3,
    resultado_e3['diametro_minimo']
)

velocidade_critica_e3 = eixo.calcula_velocidade_critica(
    resultado_e3['diametro_minimo'], 
    L_e3
)

resumo_eixos = {
    "Eixo 1 (Entrada)": {
        "Diâmetro": resultado_e1['diametro_minimo'],
        "FS Fadiga": resultado_e1['FS_fadiga'],
        "FS Escoamento": resultado_e1['FS_escoamento'],
        "Deflexão": deflexao_e1['deflexao_max_mm'],
        "Status Deflexão": deflexao_e1['status_deflexao'],
        "N_critica": velocidade_critica_e1['velocidade_critica_rpm']
    },
    "Eixo 2 (Intermediário)": {
        "Diâmetro": resultado_e2['diametro_minimo'], 
        "FS Fadiga": resultado_e2['FS_fadiga'],
        "FS Escoamento": resultado_e2['FS_escoamento'],
        "Deflexão": deflexao_e2['deflexao_max_mm'],
        "Status Deflexão": deflexao_e2['status_deflexao'],
        "N_critica": velocidade_critica_e2['velocidade_critica_rpm']
    },
    "Eixo 3 (Saída)": {
        "Diâmetro": resultado_e3['diametro_minimo'],
        "FS Fadiga": resultado_e3['FS_fadiga'],
        "FS Escoamento": resultado_e3['FS_escoamento'],
        "Deflexão": deflexao_e3['deflexao_max_mm'],
        "Status Deflexão": deflexao_e3['status_deflexao'],
        "N_critica": velocidade_critica_e3['velocidade_critica_rpm']
    }
}
fadiga1 = resultado_e1['FS_fadiga']
fadiga2 = resultado_e2['FS_fadiga']
fadiga3 = resultado_e3['FS_fadiga']
valor_min = min(fadiga1,fadiga2,fadiga3)

with open("dimensionamento_eixos.txt", "w") as f:
    f.write("DIMENSIONAMENTO DOS EIXOS - RESULTADOS FINAIS\n")
    f.write("="*60 + "\n\n")
    
    f.write(f"Material: Aço AISI 1020\n")
    f.write(f"S_ut = {S_ut} MPa, S_y = {S_y} MPa\n")
    f.write(f"Fator de segurança mínimo: {valor_min:2f}\n\n")
    
    for eixo, dados in resumo_eixos.items():
        f.write(f"{eixo}:\n")
        f.write(f"  Diâmetro mínimo: {dados['Diâmetro']:.1f} mm\n")
        f.write(f"  FS Fadiga: {dados['FS Fadiga']:.2f}\n")
        f.write(f"  FS Escoamento: {dados['FS Escoamento']:.2f}\n")
        f.write(f"  Deflexão máxima: {dados['Deflexão']:.3f} mm\n")
        f.write(f"  Status Deflexão: {dados['Status Deflexão']}\n")
        f.write(f"  Velocidade crítica: {dados['N_critica']:.0f} rpm\n\n")
    
    f.write("RECOMENDAÇÕES:\n")
    f.write("- Utilizar diâmetros comerciais (25, 30, 35, 40 mm)\n")
    f.write("- Considerar acabamento retificado nas seções críticas\n")
    f.write("- Implementar raios de concordância generosos nos degraus\n")