import equacoes_engrenagens as eng
import utilities as util
import numpy as np

parametros = util.ler_Dados_De_Entrada("Inputs/Dados_De_Entrada.txt")

S_at = parametros["S_at"]
S_ac = parametros["S_ac"]
C_p = float(parametros["C_p"])

omega_in = parametros["rotacao_motor"] * (2 * np.pi / 60)

# Torque Entrada
T_in_Nm = parametros["potencia_motor"] / omega_in
T_in_Nmm = T_in_Nm * 1000

# Torque Saída
raio_tambor_mm = parametros["diametro_tambor"] / 2
T_out_Nmm = parametros["forca_cabo"] * raio_tambor_mm

# Razão de Transmissão
i_total = T_out_Nmm / (T_in_Nmm * parametros["eficiencia"])
i_estagio = np.sqrt(i_total)

eta_estagio = np.sqrt(parametros["eficiencia"])
vida_util_horas = parametros["vida_util"]

# Chamada da Função de Dimensionamento
(
    resultado_estagio_1,
    resultado_estagio_2,
    guardaengrenagem1,
    guardaengrenagem2,
) = eng.dimensionamento_engrenagens_com_vida(
    T_in_Nmm,
    parametros["rotacao_motor"],
    i_estagio,
    eta_estagio,
    S_at,
    S_ac,
    C_p,
    vida_util_horas,
)

# Escrita dos Resultados
with open("Outputs/resultados_engrenagem.txt", "w") as f:

    if resultado_estagio_1 and resultado_estagio_2:
        i_total_efetiva_final = (
            resultado_estagio_1["i_efetiva"] * resultado_estagio_2["i_efetiva"]
        )
        erro_final = (i_total_efetiva_final - i_total) / i_total
    else:
        i_total_efetiva_final = 0
        erro_final = 1.0
        f.write("ERRO: FALHA NO DIMENSIONAMENTO DAS ENGRENAGENS\n")

    # Estágio 1
    f.write("Comparativo pares estagio 1:\n")
    f.write(
        f"{'Np':<4} | {'Nc':<4} | {'Mod':<4} | {'C (mm)':<10} | {'i_efetiva':<10} |{'i_estagio':<10} | {'Erro (%)':<10}\n"
    )
    f.write("-" * 80 + "\n")

    if guardaengrenagem1:
        guardaengrenagem1.sort(key=lambda x: x["N_p"])
        for par in guardaengrenagem1[:10]:
            erro_perc = par["erro_estagio"] * 100
            f.write(
                f"{par['N_p']:<4} | {par['N_c']:<4} | {par['m']:<4.1f} | {par['C']:<10.2f} | {par['i_efetiva']:<10.4f} |{i_estagio:<10.4f} |{erro_perc:<10.2f}\n"
            )

    # Estágio 2
    if guardaengrenagem2:
        f.write("\nComparativo pares estagio 2:\n")
        f.write(
            f"{'Np':<4} | {'Nc':<4} | {'Mod':<4} | {'C (mm)':<10} | {'i_efetiva':<10} |{'i_estagio':<10} | {'Erro (%)':<10}\n"
        )
        f.write("-" * 80 + "\n")

        guardaengrenagem2.sort(key=lambda x: x["N_p"])
        for par in guardaengrenagem2[:10]:
            erro_perc = par["erro_estagio"] * 100
            f.write(
                f"{par['N_p']:<4} | {par['N_c']:<4} | {par['m']:<4.1f} | {par['C']:<10.2f} | {par['i_efetiva']:<10.4f} |{i_estagio:<10.4f}| {erro_perc:<10.2f} |\n"
            )

    if resultado_estagio_1:
        f.write("\nComparativo modulos estagio 1 (Para o Pinhão Escolhido):\n")
        f.write("\n")
        f.write(
            f"{'Modulo':<8} | {'sigma_b (MPa)':<13} | {'FS_flex':<8} | {'sigma_c (MPa)':<13} | {'FS_pit':<8} | {'Status'}\n"
        )
        f.write("-" * 80 + "\n")

        for par in guardaengrenagem1:
            if par["N_p"] == resultado_estagio_1["N_p"]:
                status = (
                    "Nao Falha" if par["status_calc"] == "Aprovado" else "Falha"
                )
                # Garante uso da chave correta FS_pitting ou FS_superficie conforme equacoes
                f.write(
                    f"{par['m']:<8.1f} | {par['sigma_b']:<13.2f} | {par['FS_flexao']:<8.2f} | {par['sigma_c']:<13.2f} | {par.get('FS_pitting', par.get('FS_superficie', 0)):<8.2f} | {status}\n"
                )

    if resultado_estagio_2:
        f.write("\nComparativo modulos estagio 2 (Para o Pinhão Escolhido):\n")
        f.write("\n")
        f.write(
            f"{'Modulo':<8} | {'sigma_b (MPa)':<13} | {'FS_flex':<8} | {'sigma_c (MPa)':<13} | {'FS_pit':<8} | {'Status'}\n"
        )
        f.write("-" * 80 + "\n")

        for par in guardaengrenagem2:
            if par["N_p"] == resultado_estagio_2["N_p"]:
                status = (
                    "Nao Falha" if par["status_calc"] == "Aprovado" else "Falha"
                )
                f.write(
                    f"{par['m']:<8.1f} | {par['sigma_b']:<13.2f} | {par['FS_flexao']:<8.2f} | {par['sigma_c']:<13.2f} | {par.get('FS_pitting', par.get('FS_superficie', 0)):<8.2f} | {status}\n"
                )


with open("Outputs/Estagios_Engrenagem.txt", "w") as f:
    f.write(
        f"Torque de Entrada: {T_in_Nmm/1000:.2f} N.m\n Torque de Saida: {T_out_Nmm/1000:.2f} N.m\n"
    )
    f.write(
        f"Velocidade Angular de Entrada: {parametros['rotacao_motor']:.2f} rpm\n"
    )
    if i_total > 0:
        f.write(
            f"Velocidade Angular de Saida: {parametros['rotacao_motor'] / i_total:.2f} rpm\n"
        )
    f.write("\n")

    if resultado_estagio_1:
        util.escreve_estagio(f, "Estagio 1 (Entrada)", resultado_estagio_1)
    else:
        f.write("Estagio 1 (Entrada): FALHA NO DIMENSIONAMENTO\n")

    if resultado_estagio_2:
        util.escreve_estagio(f, "Estagio 2 (Saida)", resultado_estagio_2)
    else:
        f.write("Estagio 2 (Saida): FALHA NO DIMENSIONAMENTO\n")

    if resultado_estagio_1 and resultado_estagio_2:
        f.write(
            f"Razao de Transmissao Total Efetiva: {i_total_efetiva_final:.3f}\n"
        )
        f.write(f"Razao de Transmissao Total Desejada: {i_total:.3f}\n")
        f.write(f"Erro Percentual: {erro_final * 100:.2f}%\n")

        if abs(erro_final) <= 0.05:
            f.write("Dentro do limite de+/- 5%. Projeto aprovado.\n")
        else:
            f.write("Dimensionamento Falhou\n")

        # Parametros para outros setores
        f.write("\n")
        f.write("Parametros dos Eixos:\n")

        f.write(f"  Rotacao eixo1: {parametros['rotacao_motor']:.1f} rpm\n")
        f.write(f"  Torque eixo1: {T_in_Nmm/1000:.2f} N.m\n")
        f.write(
            f"    Forca Tangencial (W_t1): {resultado_estagio_1['W_t']:.2f} N\n"
        )
        f.write(
            f"    Forca Radial (W_r1): {resultado_estagio_1['W_r']:.2f} N\n"
        )

        f.write(f"  Rotacao eixo2: {resultado_estagio_1['n_p_out']:.1f} rpm\n")

        f.write(
            f"  Torque eixo2: {resultado_estagio_1['T_p_out']/1000:.2f} N.m\n"
        )
        f.write(
            f"    Forca Tangencial (W_t_c1): {resultado_estagio_1['W_t']:.2f} N\n"
        )
        f.write(
            f"    Forca Radial (W_r_c1): {resultado_estagio_1['W_r']:.2f} N\n"
        )
        f.write("  Forcas no Pinhao 2 (Eixo 2):\n")
        f.write(
            f"    Forca Tangencial (W_t_p2): {resultado_estagio_2['W_t']:.2f} N\n"
        )
        f.write(
            f"    Forca Radial (W_r_p2): {resultado_estagio_2['W_r']:.2f} N\n"
        )

        f.write(f"  Rotacao eixo3: {resultado_estagio_2['n_p_out']:.1f} rpm\n")
        f.write(
            f"  Torque eixo3: {resultado_estagio_2['T_p_out']/1000:.2f} N.m\n"
        )
        f.write(
            f"    Forca Tangencial (W_t_c2): {resultado_estagio_2['W_t']:.2f} N\n"
        )
        f.write(
            f"    Forca Radial (W_r_c2): {resultado_estagio_2['W_r']:.2f} N\n"
        )

    else:
        f.write(
            "\nSem parametros adicionais para outros setores devido a falha no dimensionamento.\n"
        )
