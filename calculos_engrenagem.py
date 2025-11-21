import numpy as np
import bib_equacoes as beta

#Dados de Entrada
P_in = 1580      # Potência de entrada (W)
n_in = 1450   # Rotacao de entrada (rpm)
F_cabo = 1988    # Forca de tracao (N)
D_tambor = 0.12  # Diametro do tambor (m)
eta_total = 0.85   # Eficiencia total do redutor
vida_util = 4000 # Vida util em nº de horas considerando Pmax

# Material Engrenagens
"""Confesso que pedi pro gemini escolher um material pq eu tava com preguiça de ficar olhando cada um, dai ele escolheu esse que a gente vai padronizar pra todas as engrenagens"""

# Aço AISI 4340 Temperado e Revenido (Material premium para off-road)
S_at = 655  # MPa (aumentado de 655)
S_ac = 1020  # MPa (aumentado de 1020) 
C_p = 191    # MPa^0.5

# Parametros iniciais
"""Nessa seção são calculados os valores que irão compor o dimensionamento dos estágios dos pares de engrenagens, obtendo aqui uma velocidade angular de entrada para o calculo do Torque de entrada, e saída, para poder calcular a razão de transmissão que vai nos dar quanto cada estágio tem que reduzir de velocidade."""

omega_in = n_in * (2*np.pi / 60) # Em rad/s mesmo
T_in = P_in / omega_in
T_out = F_cabo * (D_tambor / 2)

# Calculo da razao de transmissao total
i_total = T_out / (T_in * eta_total)
i_estagio = np.sqrt(i_total)

# Eficiencias dos estagios
eta_estagio = np.sqrt(eta_total)

"""Foram definidos modulos padronizados da tabela 12.2 do Norton assim como o Numero de dentes dos pinhões, e a largura das faces."""
# modulos_padronizados = [1.25,1.5, 2, 2.5, 3, 4, 5, 6, 8]
# N_pinhao_padronizado = [17, 18, 19, 20, 21, 22, 23, 24, 25]
b_Face = 10

resultado_estagio_1 = None
resultado_estagio_2 = None

# guardaengrenagem1 = []
# guardaengrenagem2 = []

# # Dimensionamento por partes

# for Np in N_pinhao_padronizado:
#     for m in modulos_padronizados:
#         grupo_engrenagens = beta.calcula_grupo_engrenagens(i_estagio, m, Np, b_Face)
#         forcas = beta.calcula_forcas(T_in, n_in, grupo_engrenagens, eta_estagio)
#         fatores = beta.analisa_fatores(grupo_engrenagens, forcas, S_at, S_ac, C_p)

#         dados = {**grupo_engrenagens, **forcas, **fatores}
#         # Calcula erro
#         dados['erro_estagio'] = abs((dados['i_efetiva'] - i_estagio) / i_estagio)

#         if fatores['FS_flexao'] > 1.5 and fatores['FS_pitting'] > 1.5:
#             dados['status_calc'] = "Aprovado"
#             guardaengrenagem1.append(dados)
#             break
#         else:
#             dados['status_calc'] = "Reprovado"
#             guardaengrenagem1.append(dados)


# # Escolhe a engrenagem com menor erro
# EngrenagemAprovada = [p for p in guardaengrenagem1 if p['status_calc'] == "Aprovado"]

# if EngrenagemAprovada:
#     EngrenagemAprovada.sort(key=lambda x: x['erro_estagio'])
#     resultado_estagio_1 = EngrenagemAprovada[0]
# else:
#     resultado_estagio_1 = guardaengrenagem1[-1]
#     print("Nenhuma engrenagem passou no Estagio 1")

# # Estagio 2
# if resultado_estagio_1 is not None:
#     T_in_2 = resultado_estagio_1['T_p_out']
#     n_in_2 = resultado_estagio_1['n_p_out']

#     for Np in N_pinhao_padronizado:
#         for m in modulos_padronizados:
#             grupo_engrenagens = beta.calcula_grupo_engrenagens(i_estagio, m, Np, b_Face)
#             forcas = beta.calcula_forcas(T_in_2, n_in_2, grupo_engrenagens, eta_estagio)
#             fatores = beta.analisa_fatores(grupo_engrenagens, forcas, S_at, S_ac, C_p)

#             dados = {**grupo_engrenagens, **forcas, **fatores}
#             dados['erro_estagio'] = abs((dados['i_efetiva'] - i_estagio) / i_estagio)

#             if fatores['FS_flexao'] > 1.5 and fatores['FS_pitting'] > 1.5:
#                 dados['status_calc'] = "Aprovado"
#                 guardaengrenagem2.append(dados)
#                 break
#             else:
#                 dados['status_calc'] = "Reprovado"
#                 guardaengrenagem2.append(dados)


# # Escolhe o melhor do Estagio 2
# EngrenagemAprovada2 = [p for p in guardaengrenagem2 if p['status_calc'] == "Aprovado"]

# if EngrenagemAprovada2:
#     EngrenagemAprovada2.sort(key=lambda x: x['erro_estagio'])
#     resultado_estagio_2 = EngrenagemAprovada2[0]
# else:
#     resultado_estagio_2 = guardaengrenagem2[-1]
#     print("Nenhuma engrenagem passou no Estagio 2")

# resultado_estagio_1, resultado_estagio_2, guardaengrenagem1, guardaengrenagem2 = beta.dimensionamento_engrenagens_automatico(T_in, n_in, i_estagio, eta_estagio, S_at, S_ac, C_p)
# No seu arquivo calculos_engrenagem.py, após a chamada da função:

resultado_estagio_1, resultado_estagio_2, guardaengrenagem1, guardaengrenagem2 = beta.dimensionamento_engrenagens_com_vida(
    T_in, n_in, i_estagio, eta_estagio, S_at, S_ac, C_p, vida_util_horas=4000
)

# DEBUG: Verificar por que Estágio 1 falhou
if resultado_estagio_1 is None:
    print("\n ESTÁGIO 1 FALHOU - ANALISANDO CAUSA...")
    beta.debug_estagio1(guardaengrenagem1)
    
    # Salvar todos os resultados para análise detalhada
    with open("debug_estagio1_detalhado.txt", "w") as f:
        f.write("Np,Nc,m,b_face,FS_flexao,FS_pitting,vida_flexao_horas,vida_pitting_horas,status\n")
        for config in guardaengrenagem1:
            f.write(f"{config['N_p']},{config['N_c']},{config['m']},{config['b_face']:.1f},{config['FS_flexao']:.2f},{config['FS_pitting']:.2f},{config['vida_flexao_horas']:.0f},{config['vida_pitting_horas']:.0f},{config['status_completo']}\n")
    
    print("Arquivo 'debug_estagio1_detalhado.txt' criado para análise.")
    exit(1)  # Parar execução para corrigir o problema

#Engrenagem
# def escreve_estagio(file, nome_estagio, res):
#     if res is None:
#         file.write("Nenhum modulo padrao foi suficiente.\n")
#         return

#     file.write(f"\n--- {nome_estagio}: aprovado ---\n")
#     file.write(f"  Modulo (m): {res['m']:.1f} mm\n")
#     file.write(f"  Largura da Face (b): {res['b_face']:.1f} mm\n")
#     file.write(f"  Angulo de Pressao: {20.0} graus\n")
#     file.write(f"  Razao de Contato (mp): {res['mp']:.3f}\n")
#     file.write("\n  Pinhao:\n")
#     file.write(f"    N de Dentes (Np): {res['N_p']}\n")
#     file.write(f"    Diametro Primitivo (dp): {res['d_p']:.2f} mm\n")
#     file.write(f"    Diametro Externo (dep): {res['d_ext_p']:.2f} mm\n")
#     file.write("  Coroa:\n")
#     file.write(f"    N de Dentes (Nc): {res['N_c']}\n")
#     file.write(f"    Diametro Primitivo (dc): {res['d_c']:.2f} mm\n")
#     file.write(f"    Diametro Externo (dec): {res['d_ext_c']:.2f} mm\n")
#     file.write("  Conjunto:\n")
#     file.write(f"    Razao de Transmissao (i): {res['i_efetiva']:.3f}\n")
#     file.write(f"    Distancia entre Centros (C): {res['C']:.2f} mm\n")
#     file.write("  Fatores de Seguranca:\n")
#     file.write(f"    FS (Flexao): {res['FS_flexao']:.2f}\n")
#     file.write(f"    FS (Superficie): {res['FS_pitting']:.2f}\n")
def escreve_estagio(file, nome_estagio, res):
    """
    Escreve resultados incluindo informações de vida útil
    """
    if res is None:
        file.write("Nenhum modulo padrao foi suficiente.\n")
        return

    file.write(f"\n--- {nome_estagio}: {res['status_completo']} ---\n")
    file.write(f"  Modulo (m): {res['m']:.1f} mm\n")
    file.write(f"  Largura da Face (b): {res['b_face']:.1f} mm\n")
    file.write(f"  Angulo de Pressao: {20.0} graus\n")
    file.write(f"  Razao de Contato (mp): {res['mp']:.3f}\n")
    file.write("\n  Pinhao:\n")
    file.write(f"    N de Dentes (Np): {res['N_p']}\n")
    file.write(f"    Diametro Primitivo (dp): {res['d_p']:.2f} mm\n")
    file.write(f"    Diametro Externo (dep): {res['d_ext_p']:.2f} mm\n")
    file.write("  Coroa:\n")
    file.write(f"    N de Dentes (Nc): {res['N_c']}\n")
    file.write(f"    Diametro Primitivo (dc): {res['d_c']:.2f} mm\n")
    file.write(f"    Diametro Externo (dec): {res['d_ext_c']:.2f} mm\n")
    file.write("  Conjunto:\n")
    file.write(f"    Razao de Transmissao (i): {res['i_efetiva']:.3f}\n")
    file.write(f"    Distancia entre Centros (C): {res['C']:.2f} mm\n")
    file.write("  Fatores de Seguranca:\n")
    file.write(f"    FS (Flexao): {res['FS_flexao']:.2f}\n")
    file.write(f"    FS (Superficie): {res['FS_pitting']:.2f}\n")
    file.write("  Vida Util Estimada:\n")
    file.write(f"    Vida por Flexao: {res['vida_flexao_horas']:.0f} horas\n")
    file.write(f"    Vida por Pitting: {res['vida_pitting_horas']:.0f} horas\n")
    file.write(f"    Vida Minima: {res['vida_util_minima']:.0f} horas\n")
    file.write(f"    Status Vida: {res['status_vida']}\n")

with open("resultados_engrenagem.txt", "w") as f:
    if resultado_estagio_1 and resultado_estagio_2:
        i_total_efetiva_final = resultado_estagio_1['i_efetiva'] * resultado_estagio_2['i_efetiva']
        erro_final = (i_total_efetiva_final - i_total) / i_total

    # Trecho de comparação entre pares de engrenagens
    f.write("Comparativo pares estagio 1:\n")
    f.write(f"{'Np':<4} | {'Nc':<4} | {'Mod':<4} | {'C (mm)':<10} | {'i_efetiva':<10} |{'i_estagio':<10} | {'Erro (%)':<10}\n")
    f.write("-" * 80 + "\n")

    guardaengrenagem1.sort(key=lambda x: x['N_p'])

    for par in guardaengrenagem1[:10]:
        status = "Melhor Par" if par == resultado_estagio_1 else "Possivel Par"
        erro_perc = par['erro_estagio'] * 100

        f.write(f"{par['N_p']:<4} | {par['N_c']:<4} | {par['m']:<4.1f} | {par['C']:<10.2f} | {par['i_efetiva']:<10.4f} |{i_estagio:<10.4f} |{erro_perc:<10.2f}\n")

    if guardaengrenagem2:
        f.write("\n")
        f.write("Comparativo pares estagio 2:\n")
        f.write(f"{'Np':<4} | {'Nc':<4} | {'Mod':<4} | {'C (mm)':<10} | {'i_efetiva':<10} |{'i_estagio':<10} | {'Erro (%)':<10}\n")
        f.write("-" * 80 + "\n")

        guardaengrenagem2.sort(key=lambda x: x['N_p'])
        for par in guardaengrenagem2[:10]:
            status = "Melhor Par" if par == resultado_estagio_2 else "Possivel Par"
            erro_perc = par['erro_estagio'] * 100
            f.write(f"{par['N_p']:<4} | {par['N_c']:<4} | {par['m']:<4.1f} | {par['C']:<10.2f} | {par['i_efetiva']:<10.4f} |{i_estagio:<10.4f}| {erro_perc:<10.2f} |\n")

# Trecho de comparação de módulos e tensões
    f.write("\nComparativo modulos estagio 1:\n")
    if guardaengrenagem1:
        f.write("\n")
        f.write(f"{'Modulo':<8} | {'sigma_b (MPa)':<13} | {'FS_flex':<8} | {'sigma_c (MPa)':<13} | {'FS_pit':<8} | {'Status'}\n")
        f.write("-" * 80 + "\n")

        for par in guardaengrenagem1:
            if par['N_p'] == resultado_estagio_1['N_p']:
                status = "Nao Falha" if par['status_calc'] == "Aprovado" else "Falha"
                f.write(f"{par['m']:<8.1f} | {par['sigma_b']:<13.2f} | {par['FS_flexao']:<8.2f} | {par['sigma_c']:<13.2f} | {par['FS_pitting']:<8.2f} | {status}\n")

    f.write("\nComparativo modulos estagio 2:\n")
    if guardaengrenagem2:
        f.write("\n")
        f.write(f"{'Modulo':<8} | {'sigma_b (MPa)':<13} | {'FS_flex':<8} | {'sigma_c (MPa)':<13} | {'FS_pit':<8} | {'Status'}\n")
        f.write("-" * 80 + "\n")

        for par in guardaengrenagem2:
            if par['N_p'] == resultado_estagio_2['N_p']:
                status = "Nao Falha" if par['status_calc'] == "Aprovado" else "Falha"
                f.write(f"{par['m']:<8.1f} | {par['sigma_b']:<13.2f} | {par['FS_flexao']:<8.2f} | {par['sigma_c']:<13.2f} | {par['FS_pitting']:<8.2f} | {status}\n")

with open("Estagios_Engrenagem.txt", "w") as f:
    f.write(f"Torque de Entrada: {T_in:.2f} N.m\n Torque de Saida: {T_out:.2f} N.m\n")
    f.write("\n")

    escreve_estagio(f, "Estagio 1 (Entrada)", resultado_estagio_1)
    escreve_estagio(f, "Estagio 2 (Saida)", resultado_estagio_2)

    f.write(f"Razao de Transmissao Total Efetiva: {i_total_efetiva_final:.3f}\n")
    f.write(f"Razao de Transmissao Total Desejada: {i_total:.3f}\n")
    f.write(f"Erro Percentual: {erro_final * 100:.2f}%\n")

    if abs(erro_final) <= 0.05:
            f.write("Dentro do limite de+/- 5%. Projeto aprovado.\n")
    else:
            f.write("Fora do limite de +/- 5%. Projeto com falha.\n")
            f.write("Dimensionamento Falhou\n")


    # ---  Parametros para outros setores ---
    if resultado_estagio_1 and resultado_estagio_2:
        # Eixo 1 Entrada
        f.write("\n")
        f.write("Parametros dos Eixos:\n")

        f.write(f"  Rotacao (n_eixo1): {n_in:.1f} rpm\n")
        f.write(f"  Torque (T_eixo1): {T_in:.2f} N.m\n")
        f.write("  Forcas no Pinhao 1 (Eixo 1):\n")
        f.write(f"    Forca Tangencial (W_t1): {resultado_estagio_1['W_t']:.2f} N\n")
        f.write(f"    Forca Radial (W_r1): {resultado_estagio_1['W_r']:.2f} N\n")

        # Eixo 2
        f.write(f"  Rotacao (n_eixo2): {resultado_estagio_1['n_p_out']:.1f} rpm\n")
        f.write(f"  Torque (T_eixo2): {resultado_estagio_1['T_p_out']:.2f} N.m\n")
        f.write("  Forcas da Coroa 1 (Eixo 2):\n")
        f.write(f"    Forca Tangencial (W_t_c1): {resultado_estagio_1['W_t']:.2f} N\n")
        f.write(f"    Forca Radial (W_r_c1): {resultado_estagio_1['W_r']:.2f} N\n")
        f.write("  Forcas no Pinhao 2 (Eixo 2):\n")
        f.write(f"    Forca Tangencial (W_t_p2): {resultado_estagio_2['W_t']:.2f} N\n")
        f.write(f"    Forca Radial (W_r_p2): {resultado_estagio_2['W_r']:.2f} N\n")

        # Eixo 3 Saida
        f.write(f"  Rotacao (n_eixo3): {resultado_estagio_2['n_p_out']:.1f} rpm\n")
        f.write(f"  Torque (T_eixo3): {resultado_estagio_2['T_p_out']:.2f} N.m\n")
        f.write("  Forcas da Coroa 2 (Eixo 3):\n")
        f.write(f"    Forca Tangencial (W_t_c2): {resultado_estagio_2['W_t']:.2f} N\n")
        f.write(f"    Forca Radial (W_r_c2): {resultado_estagio_2['W_r']:.2f} N\n")

    else:
        f.write("\nSem parametros adicionais para outros setores devido a falha no dimensionamento.\n")

# Eixos!

# Aço AISI 1020 Normalizado
S_ut = 469  # MPa - Tensão de ruptura
S_y = 390   # MPa - Tensão de escoamento

# Coeficientes
C_seg = 1.5  # Fator de segurança
C_superf, acabamento_superf = beta.calcula_c_superf(S_ut) #fator de acabamento superficial e o acabamento superficial definido
C_conf = 0.814  # Fator de confiabilidade de 99%
C_temp = 1.0   # Fator de temperatura para T < 450C
C_carr = 1.0   # Fator de carregamento para carga moderada

# Defina distâncias estimadas (em mm)
distancias_eixo1 = {'LA': 50, 'AB': 200, 'BC': 150}
distancias_eixo2 = {'LA': 40, 'L1': 120, 'L2': 40}  # AB = 200mm
distancias_eixo3 = {'LA': 60, 'AB': 220, 'BC': 160}

# Calcule as reações
reações_eixo1 = beta.calcula_reações_mancais_eixo1(
    {'W_t': 462.46, 'W_r': 168.32},  # Pinhão 1
    distancias_eixo1
)

reações_eixo2 = beta.calcula_reações_eixo2(
    {'W_t': 462.46, 'W_r': 168.32},  # Coroa 1
    {'W_t': 977.10, 'W_r': 355.64},  # Pinhão 2  
    distancias_eixo2
)

reações_eixo3 = beta.calcula_reações_mancais_eixo1(
    {'W_t': 977.10, 'W_r': 355.64},  # Coroa 2
    distancias_eixo3
)

print("=== REAÇÕES NOS MANCAIS ===")
print(f"Eixo 1 - Mancal A: {reações_eixo1['R_A_resultante']:.1f} N")
print(f"Eixo 1 - Mancal B: {reações_eixo1['R_B_resultante']:.1f} N")
print(f"Eixo 2 - Mancal A: {reações_eixo2['R_A_resultante']:.1f} N") 
print(f"Eixo 2 - Mancal B: {reações_eixo2['R_B_resultante']:.1f} N")
print(f"Eixo 3 - Mancal A: {reações_eixo3['R_A_resultante']:.1f} N")
print(f"Eixo 3 - Mancal B: {reações_eixo3['R_B_resultante']:.1f} N")

# Dados dos eixos baseados nos seus resultados
distancias_eixo1 = {'LA': 50, 'AB': 200, 'BC': 150}
distancias_eixo2 = {'LA': 40, 'L1': 120, 'L2': 40}  # AB = 200mm
distancias_eixo3 = {'LA': 60, 'AB': 220, 'BC': 160}

# EIXO 1 - ENTRADA
print("=== EIXO 1 - DIAGRAMAS ===")
resultados_eixo1 = beta.calcula_diagramas_eixo_simples(
    forcas_engrenagem={'W_t': 462.46, 'W_r': 168.32},
    distancias=distancias_eixo1,
    torque=10.41,
    nome_eixo="Eixo 1 - Entrada"
)

# EIXO 2 - INTERMEDIÁRIO  
print("=== EIXO 2 - DIAGRAMAS ===")
resultados_eixo2 = beta.calcula_diagramas_eixo_duplo(
    forcas_engrenagem1={'W_t': 462.46, 'W_r': 168.32},  # Coroa 1
    forcas_engrenagem2={'W_t': 977.10, 'W_r': 355.64},  # Pinhão 2
    distancias=distancias_eixo2,
    torque=35.18,
    nome_eixo="Eixo 2 - Intermediário"
)

# EIXO 3 - SAÍDA
print("=== EIXO 3 - DIAGRAMAS ===")
resultados_eixo3 = beta.calcula_diagramas_eixo_simples(
    forcas_engrenagem={'W_t': 977.10, 'W_r': 355.64},
    distancias=distancias_eixo3, 
    torque=118.91,
    nome_eixo="Eixo 3 - Saída"
)

# RESUMO DOS VALORES CRÍTICOS
print("\n" + "="*50)
print("RESUMO - VALORES CRÍTICOS PARA DIMENSIONAMENTO")
print("="*50)
print(f"EIXO 1 - Momento máximo: {resultados_eixo1['M_max']:.0f} N.mm")
print(f"EIXO 1 - Cortante máximo: {resultados_eixo1['V_max']:.0f} N")
print(f"EIXO 1 - Reação A: {resultados_eixo1['R_A_resultante']:.0f} N")
print(f"EIXO 1 - Reação B: {resultados_eixo1['R_B_resultante']:.0f} N")

print(f"\nEIXO 2 - Momento máximo: {resultados_eixo2['M_max']:.0f} N.mm")
print(f"EIXO 2 - Cortante máximo: {resultados_eixo2['V_max']:.0f} N")  
print(f"EIXO 2 - Reação A: {resultados_eixo2['R_A_resultante']:.0f} N")
print(f"EIXO 2 - Reação B: {resultados_eixo2['R_B_resultante']:.0f} N")

print(f"\nEIXO 3 - Momento máximo: {resultados_eixo3['M_max']:.0f} N.mm")
print(f"EIXO 3 - Cortante máximo: {resultados_eixo3['V_max']:.0f} N")
print(f"EIXO 3 - Reação A: {resultados_eixo3['R_A_resultante']:.0f} N")
print(f"EIXO 3 - Reação B: {resultados_eixo3['R_B_resultante']:.0f} N")

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
    
    for eixo, dados in resultados_diagramas.items():
        f.write(f"{eixo.upper()}:\n")
        f.write(f"  Momento máximo: {dados['M_max']:.0f} N.mm\n")
        f.write(f"  Cortante máximo: {dados['V_max']:.0f} N\n")
        f.write(f"  Reação A: {dados['R_A_resultante']:.0f} N\n")
        f.write(f"  Reação B: {dados['R_B_resultante']:.0f} N\n")
        f.write(f"  Posição M_max: {dados['posicao_M_max']:.1f} mm\n\n")