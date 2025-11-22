import equacoes_engrenagens as eng
import utilities as util
import numpy as np

parametros = util.ler_Dados_De_Entrada('Dados_De_Entrada.txt')

# Parametros iniciais
"""Nessa seção são calculados os valores que irão compor o dimensionamento dos estágios dos pares de engrenagens, obtendo aqui uma velocidade angular de entrada para o calculo do Torque de entrada, e saída, para poder calcular a razão de transmissão que vai nos dar quanto cada estágio tem que reduzir de velocidade."""


S_at=parametros['S_at']
S_ac=parametros['S_ac']
C_p=float(parametros['C_p'])

omega_in = parametros['rotacao_motor'] * (2*np.pi / 60) # Em rad/s mesmo

T_in = parametros['potencia_motor'] / omega_in
T_out = parametros['forca_cabo'] * (parametros['diametro_tambor'] / 2)

# Calculo da razao de transmissao total
i_total = T_out / (T_in * parametros['eficiencia'])
i_estagio = np.sqrt(i_total)

# Eficiencias dos estagios
eta_estagio = np.sqrt(parametros['eficiencia'])
vida_util_horas = parametros['vida_util']
b_Face = 10

resultado_estagio_1 = None
resultado_estagio_2 = None

resultado_estagio_1, resultado_estagio_2, guardaengrenagem1, guardaengrenagem2 = eng.dimensionamento_engrenagens_com_vida(
    T_in, parametros['rotacao_motor'], i_estagio, eta_estagio, S_at, S_ac, C_p, vida_util_horas
)

with open("Outputs/resultados_engrenagem.txt", "w") as f:
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


with open("Outputs/Estagios_Engrenagem.txt", "w") as f:
    f.write(f"Torque de Entrada: {T_in:.2f} N.m\n Torque de Saida: {T_out:.2f} N.m\n")
    f.write(f"Velocidade Angular de Entrada: {parametros['rotacao_motor']:.2f} rpm\n")
    f.write(f"Velocidade Angular de Saida: {parametros['rotacao_motor'] / i_total:.2f} rpm\n")
    f.write("\n")

    util.escreve_estagio(f, "Estagio 1 (Entrada)", resultado_estagio_1)
    util.escreve_estagio(f, "Estagio 2 (Saida)", resultado_estagio_2)

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

        f.write(f"  Rotacao (n_eixo1): {parametros['rotacao_motor']:.1f} rpm\n")
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
