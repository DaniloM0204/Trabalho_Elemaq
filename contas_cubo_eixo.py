import equacoes_cubo_eixo as ce
import contas_eixo as eixo
import utilities as util

param = util.ler_Estagios_Engrenagem('Outputs/Estagios_Engrenagem.txt')

largura_e1 = 20.0
largura_e2 = 30.0

# Tenta ler larguras reais
try:
    with open('Outputs/resultados_engrenagem.txt', 'r') as f:
        content = f.read()
        largura_e1 = 18.0
        largura_e2 = 30.0
except FileNotFoundError:
    pass

# Torques (N.m)
T1 = 10.41
T2 = 35.18
T3 = 118.91

d1 = eixo.d1_com # ex: 20mm
d2 = eixo.d2_com # ex: 25mm
d3 = eixo.d3_com # ex: 30mm


# Eixo 1 - Pinhao 1
# Torque T1, Diametro d1, Largura b1 (Estagio 1)
chaveta1 = ce.dimensiona_chavetas(d1, T1*1000, largura_e1, "Chaveta Eixo1-Pinhao1")

# Eixo 2 - Coroa 1
# Torque T2, Diametro d2, Largura b1 (Estagio 1)
chaveta2a = ce.dimensiona_chavetas(d2, T2*1000, largura_e1, "Chaveta Eixo2-Coroa1")

# Eixo 2 - Pinhao 2
# Torque T2, Diametro d2, Largura b2 (Estagio 2)
chaveta2b = ce.dimensiona_chavetas(d2, T2*1000, largura_e2, "Chaveta Eixo2-Pinhao2")

# Eixo 3 - Coroa 2
# Torque T3, Diametro d3, Largura b2 (Estagio 2)
chaveta3 = ce.dimensiona_chavetas(d3, T3*1000, largura_e2, "Chaveta Eixo3-Coroa2")

# Lista para verificacao
lista_chavetas_bruta = [chaveta1, chaveta2a, chaveta2b, chaveta3]
lista_verificada = []

for ch in lista_chavetas_bruta:
    res = ce.verifica_cisalhamento_chaveta(ch)
    res = ce.verifica_esmagamento_chaveta(res)
    lista_verificada.append(res)


with open("Outputs/dimensionamento_chavetas.txt", "w") as f:
    f.write("Material: Aco AISI 1020 (Sy=390 MPa)\n")
    f.write("Fator de seguranca: 1.5\n\n")

    todos_ok = True

    for ch in lista_verificada:
        f.write(f"{ch['nome']}:\n")
        f.write(f"  Eixo: {ch['diametro_eixo']:.1f} mm | Torque: {ch['torque']:.0f} N.mm\n")
        f.write(f"  Chaveta: {ch['dimensoes_chaveta']['b']}x{ch['dimensoes_chaveta']['h']}x{ch['comprimento']:.1f} mm\n")
        f.write(f"  Cisalhamento: {ch['tau_cisalhamento']:.1f} MPa (FS={ch['FS_cisalhamento']:.2f}) - {ch['status_cisalhamento']}\n")
        f.write(f"  Esmagamento:  {ch['sigma_esmagamento']:.1f} MPa (FS={ch['FS_esmagamento']:.2f}) - {ch['status_esmagamento']}\n")
        f.write("-" * 40 + "\n")

        if ch['status_cisalhamento'] == "FALHA" or ch['status_esmagamento'] == "FALHA":
            todos_ok = False

    if todos_ok:
        f.write("Chavetas atendem o projeto\n")
    else:
        f.write("Chavetas nao atendem o projeto\n")
