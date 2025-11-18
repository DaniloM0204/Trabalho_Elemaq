import numpy as np
import bib_equacoes as β

#Dados de Entrada
P_in = 1580      # Potência de entrada
n_in = 1450   # Rotacao de entrada
F_cabo = 1988    # Forca de tracao
D_tambor = 0.12  # Diametro do tambor
η_total = 0.85   # Eficiencia total do redutor

# Parametros iniciais
ω_in = n_in * (2 * np.pi / 60) * (D_tambor / 2)
T_in = P_in / ω_in
T_out = F_cabo * (D_tambor / 2)
i_total = T_out / (T_in * η_total)
i_alvo_1 = np.sqrt(i_total)
i_alvo_2 = np.sqrt(i_total)
η_1 = np.sqrt(η_total)
η_2 = np.sqrt(η_total)



estagio1_iteracoes = []
estagio2_iteracoes = []

modulos_padronizados = [1.5, 2, 2.5, 3, 4, 5, 6, 8]
Np_1 = 18
b_Face_1 = 10
Np_2 = 18
b_Face_2 = 10

resultado_estagio_1 = None
resultado_estagio_2 = None

#
def escreve_estagio(file, nome_estagio, res):
    if res is None:
        file.write("Nenhum modulo padrao foi suficiente.\n")
        return

    file.write(f"\n--- {nome_estagio}: aprovado ---\n")
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


# Estagio 1
for modulo in modulos_padronizados:
    resultado = β.dimensionar_estagio(T_in, n_in, i_alvo_1, modulo, Np_1, b_Face_1, η_1)

    if resultado['FS_flexao'] > 1.5 and resultado['FS_pitting'] > 1.5:
        resultado_estagio_1 = resultado
        break

# Estágio 2
if resultado_estagio_1 is not None:
    # Saída do estágio 1 é a entrada do estágio 2
    T_p2_in = resultado_estagio_1['T_p_out']
    n_p2_in = resultado_estagio_1['n_p_out']

    for modulo in modulos_padronizados:
        resultado = β.dimensionar_estagio(T_p2_in, n_p2_in, i_alvo_2, modulo, Np_2, b_Face_2, η_2)

        if resultado['FS_flexao'] > 1.5 and resultado['FS_pitting'] > 1.5:
            resultado_estagio_2 = resultado
            break

with open("resultados_engrenagem.txt", "w") as f:
    # Estágio 1
    escreve_estagio(f, "Estagio 1 (Entrada)", resultado_estagio_1)

    # Estágio 2
    escreve_estagio(f, "Estagio 2 (Saida)", resultado_estagio_2)

    if resultado_estagio_1 and resultado_estagio_2:
        i_total_efetiva_final = resultado_estagio_1['i_efetiva'] * resultado_estagio_2['i_efetiva']
        erro_final = (i_total_efetiva_final - i_total) / i_total

        f.write(f"Razao de Transmissao Total Efetiva: {i_total_efetiva_final:.3f}\n")
        f.write(f"Erro Percentual: {erro_final * 100:.2f}%\n")

        if abs(erro_final) <= 0.05:
            f.write("Dentro do limite de+/- 5%. Projeto aprovado.\n")
        else:
            f.write("Fora do limite de +/- 5%. Projeto com falha.\n")
    else:
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
