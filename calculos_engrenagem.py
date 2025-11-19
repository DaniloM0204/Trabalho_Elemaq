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

# Aço AISI 4140 Normalizado
S_at = 240  # MPa
S_ac = 860  # MPa
C_p = 191   # MPa^0.5

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
modulos_padronizados = [1.5, 2, 2.5, 3, 4, 5, 6, 8]
Np = 18
b_Face = 10

resultado_estagio_1 = None
resultado_estagio_2 = None

# Dimensionamento por partes

# Estagio 1
for m in modulos_padronizados:
    # Formato da Engrenagem
    grupo_engrenagens = beta.calcula_grupo_engrenagens(i_estagio, m, Np, b_Face)

    # Forças nas engrenagens
    forcas = beta.calcula_forcas(T_in, n_in, grupo_engrenagens, eta_estagio)

    # Análise dos fatores de segurança
    fatores = beta.analisa_fatores(grupo_engrenagens, forcas, S_at, S_ac, C_p)

    # Verifica se o fator de seguranca é maior que 1.5
    if fatores['FS_flexao'] > 1.5 and fatores['FS_pitting'] > 1.5:
        resultado_estagio_1 = {**grupo_engrenagens, **forcas, **fatores}
        break

# Estagio 2
if resultado_estagio_1 is not None:
    T_in_2 = resultado_estagio_1['T_p_out']
    n_in_2 = resultado_estagio_1['n_p_out']

    for m in modulos_padronizados:
        grupo_engrenagens = beta.calcula_grupo_engrenagens(i_estagio, m, Np, b_Face)
        forcas = beta.calcula_forcas(T_in_2, n_in_2, grupo_engrenagens, eta_estagio)
        fatores = beta.analisa_fatores(grupo_engrenagens, forcas, S_at, S_ac, C_p)

        if fatores['FS_flexao'] > 1.5 and fatores['FS_pitting'] > 1.5:
            resultado_estagio_2 = {**grupo_engrenagens,**forcas ,**fatores}

            break


#Engrenagem
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




with open("resultados_engrenagem.txt", "w") as f:

    f.write(f"Torque de Entrada: {T_in:.2f} N.m\n Torque de Saida: {T_out:.2f} N.m\n")

    escreve_estagio(f, "Estagio 1 (Entrada)", resultado_estagio_1)
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