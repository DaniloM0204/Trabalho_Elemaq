import equacoes_cubo_eixo as ce
import utilities as util
import re


# --- FUNÇÃO AUXILIAR PARA LER DADOS REAIS DOS ARQUIVOS ---
def ler_dados_projeto():
    # 1. Ler Diâmetros Comerciais dos Eixos
    diams = {}
    try:
        with open("Outputs/dimensionamento_eixos.txt", "r") as f:
            texto = f.read()
            # Procura por "Diametros comerciais sugeridos: 20mm, 25mm, 30mm"
            match = re.search(
                r"Diametros comerciais sugeridos:\s*(\d+)mm,\s*(\d+)mm,\s*(\d+)mm",
                texto,
            )
            if match:
                diams["d1"] = float(match.group(1))
                diams["d2"] = float(match.group(2))
                diams["d3"] = float(match.group(3))
            else:
                # Fallback se a regex falhar
                print(
                    "AVISO: Não foi possível ler diâmetros comerciais. Usando padrão."
                )
                diams = {"d1": 20.0, "d2": 25.0, "d3": 30.0}
    except FileNotFoundError:
        diams = {"d1": 20.0, "d2": 25.0, "d3": 30.0}

    # 2. Ler Larguras das Engrenagens
    larguras = {"b1": 20.0, "b2": 30.0}  # Valores default seguros
    try:
        with open("Outputs/Estagios_Engrenagem.txt", "r") as f:
            texto = f.read()
            # Procura por "Largura da Face (b): XX.X mm"
            matches = re.findall(
                r"Largura da Face \(b\):\s*([\d\.]+)\s*mm", texto
            )
            if len(matches) >= 2:
                larguras["b1"] = float(matches[0])  # Estágio 1
                larguras["b2"] = float(matches[1])  # Estágio 2
    except FileNotFoundError:
        pass

    return diams, larguras


# --- INÍCIO DO CÁLCULO ---

# Ler parâmetros do sistema
param_eng = util.ler_Estagios_Engrenagem("Outputs/Estagios_Engrenagem.txt")
diams, larguras = ler_dados_projeto()

# Dados extraídos
largura_e1 = larguras["b1"]
largura_e2 = larguras["b2"]
d1 = diams["d1"]
d2 = diams["d2"]
d3 = diams["d3"]

# Torques (N.m) - Pegando do input ou recalculando
T1 = 10.41
T2 = 35.18
T3 = 118.91

print(
    f"Dimensionando Chavetas com: d1={d1}, d2={d2}, d3={d3} | b1={largura_e1}, b2={largura_e2}"
)

# Dimensionamento
# Eixo 1 - Pinhao 1
chaveta1 = ce.dimensiona_chavetas(
    d1, T1 * 1000, largura_e1, "Chaveta Eixo1-Pinhao1"
)

# Eixo 2 - Coroa 1 (Entrada do Eixo 2)
chaveta2a = ce.dimensiona_chavetas(
    d2, T2 * 1000, largura_e1, "Chaveta Eixo2-Coroa1"
)

# Eixo 2 - Pinhao 2 (Saída do Eixo 2)
chaveta2b = ce.dimensiona_chavetas(
    d2, T2 * 1000, largura_e2, "Chaveta Eixo2-Pinhao2"
)

# Eixo 3 - Coroa 2 (Entrada do Eixo 3)
chaveta3 = ce.dimensiona_chavetas(
    d3, T3 * 1000, largura_e2, "Chaveta Eixo3-Coroa2"
)

# Verificação
lista_chavetas = [chaveta1, chaveta2a, chaveta2b, chaveta3]
lista_verificada = []

for ch in lista_chavetas:
    res = ce.verifica_cisalhamento_chaveta(ch)
    res = ce.verifica_esmagamento_chaveta(res)
    lista_verificada.append(res)

# Escrita do Relatório
with open("Outputs/dimensionamento_chavetas.txt", "w") as f:
    f.write("=== DIMENSIONAMENTO DE CHAVETAS (ABNT/DIN 6885) ===\n")
    f.write("Material: Aco AISI 1020 (Sy=390 MPa)\n")
    f.write("Fator de seguranca: 1.5\n")
    f.write(f"Dimensoes Base: d1={d1}mm, d2={d2}mm, d3={d3}mm\n\n")

    todos_ok = True
    for ch in lista_verificada:
        f.write(f"{ch['nome']}:\n")
        f.write(
            f"  Eixo: {ch['diametro_eixo']:.1f} mm | Torque: {ch['torque']:.0f} N.mm\n"
        )

        dims = ch["dimensoes_chaveta"]
        f.write(
            f"  Chaveta Selecionada: {dims['b']} x {dims['h']} x {ch['comprimento']:.1f} mm\n"
        )

        # Formatação bonita dos status
        s_cis = ch["status_cisalhamento"]
        s_esm = ch["status_esmagamento"]
        f.write(
            f"  Cisalhamento: {ch['tau_cisalhamento']:.1f} MPa (Adm: {ch['tau_admissivel']:.1f}) -> [{s_cis}]\n"
        )
        f.write(
            f"  Esmagamento:  {ch['sigma_esmagamento']:.1f} MPa (Adm: {ch['sigma_adm_esmagamento']:.1f}) -> [{s_esm}]\n"
        )
        f.write("-" * 50 + "\n")

        if s_cis == "FALHA" or s_esm == "FALHA":
            todos_ok = False

    if todos_ok:
        f.write(
            "\nCONCLUSAO: Todas as chavetas atendem aos criterios de projeto.\n"
        )
    else:
        f.write(
            "\nCONCLUSAO: ATENCAO! Ha chavetas com falha. Verifique diametros ou aumente o comprimento.\n"
        )

print(
    "Dimensionamento de chavetas concluido. Verifique 'Outputs/dimensionamento_chavetas.txt'"
)
