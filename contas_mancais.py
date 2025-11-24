import equacoes_mancais as mancal
import csv
import os
import re


# --- 1. FUNÇÃO DE LEITURA DE DADOS (REGEX) ---
def ler_dados_para_mancais():
    dados = {
        "n1": 1750.0,
        "n2": 400.0,
        "n3": 100.0,  # Defaults
        "d1": 20,
        "d2": 25,
        "d3": 30,
        "R1A": 0,
        "R1B": 0,
        "R2A": 0,
        "R2B": 0,
        "R3A": 0,
        "R3B": 0,
    }

    # A. Ler RPMs (Estagios_Engrenagem.txt)
    try:
        with open("Outputs/Estagios_Engrenagem.txt", "r") as f:
            txt = f.read()
            m = re.search(r"Rotacao eixo1:\s*([\d\.]+)", txt)
            if m:
                dados["n1"] = float(m.group(1))
            m = re.search(r"Rotacao eixo2:\s*([\d\.]+)", txt)
            if m:
                dados["n2"] = float(m.group(1))
            m = re.search(r"Rotacao eixo3:\s*([\d\.]+)", txt)
            if m:
                dados["n3"] = float(m.group(1))
    except:
        pass

    # B. Ler Diametros Comerciais (dimensionamento_eixos.txt)
    try:
        with open("Outputs/dimensionamento_eixos.txt", "r") as f:
            txt = f.read()
            m = re.search(
                r"Diametros comerciais sugeridos:\s*(\d+)mm,\s*(\d+)mm,\s*(\d+)mm",
                txt,
            )
            if m:
                dados["d1"] = int(m.group(1))
                dados["d2"] = int(m.group(2))
                dados["d3"] = int(m.group(3))
    except:
        pass

    # C. Ler Reações (resultados_diagramas.txt)
    # Esse arquivo tem estrutura: "EIXO1:\n ... Reacao A: 500 N ... Reacao B: 600 N"
    try:
        with open("Outputs/resultados_diagramas.txt", "r") as f:
            content = f.read()
            # Regex esperta para pegar blocos
            blocos = content.split("EIXO")
            for b in blocos:
                if "1" in b[:5]:  # Bloco do Eixo 1
                    ra = re.search(r"Reacao A:\s*([\d\.]+)", b)
                    rb = re.search(r"Reacao B:\s*([\d\.]+)", b)
                    if ra:
                        dados["R1A"] = float(ra.group(1))
                    if rb:
                        dados["R1B"] = float(rb.group(1))
                if "2" in b[:5]:  # Bloco do Eixo 2
                    ra = re.search(r"Reacao A:\s*([\d\.]+)", b)
                    rb = re.search(r"Reacao B:\s*([\d\.]+)", b)
                    if ra:
                        dados["R2A"] = float(ra.group(1))
                    if rb:
                        dados["R2B"] = float(rb.group(1))
                if "3" in b[:5]:  # Bloco do Eixo 3
                    ra = re.search(r"Reacao A:\s*([\d\.]+)", b)
                    rb = re.search(r"Reacao B:\s*([\d\.]+)", b)
                    if ra:
                        dados["R3A"] = float(ra.group(1))
                    if rb:
                        dados["R3B"] = float(rb.group(1))
    except:
        print(
            "AVISO: Não foi possível ler as reações. Verifique 'Outputs/resultados_diagramas.txt'"
        )

    return dados


# --- 2. CATÁLOGO EXPANDIDO (FALLBACK) ---
# Adicionei séries 60, 62 e 63 para diâmetros de 15 a 60mm
catalogo_fallback = [
    # d=15
    {"codigo": "6002", "d": 15, "C_din": 5.85},
    {"codigo": "6202", "d": 15, "C_din": 8.06},
    {"codigo": "6302", "d": 15, "C_din": 11.9},
    # d=17
    {"codigo": "6003", "d": 17, "C_din": 6.37},
    {"codigo": "6203", "d": 17, "C_din": 9.95},
    {"codigo": "6303", "d": 17, "C_din": 14.3},
    # d=20
    {"codigo": "6004", "d": 20, "C_din": 9.95},
    {"codigo": "6204", "d": 20, "C_din": 13.5},
    {"codigo": "6304", "d": 20, "C_din": 16.8},
    # d=25
    {"codigo": "6005", "d": 25, "C_din": 11.9},
    {"codigo": "6205", "d": 25, "C_din": 14.8},
    {"codigo": "6305", "d": 25, "C_din": 23.6},
    # d=30
    {"codigo": "6006", "d": 30, "C_din": 13.8},
    {"codigo": "6206", "d": 30, "C_din": 20.3},
    {"codigo": "6306", "d": 30, "C_din": 29.6},
    # d=35
    {"codigo": "6007", "d": 35, "C_din": 16.8},
    {"codigo": "6207", "d": 35, "C_din": 27.0},
    {"codigo": "6307", "d": 35, "C_din": 35.1},
    # d=40
    {"codigo": "6008", "d": 40, "C_din": 17.8},
    {"codigo": "6208", "d": 40, "C_din": 30.7},
    {"codigo": "6308", "d": 40, "C_din": 42.3},
    # d=45
    {"codigo": "6009", "d": 45, "C_din": 22.1},
    {"codigo": "6209", "d": 45, "C_din": 33.2},
    {"codigo": "6309", "d": 45, "C_din": 55.3},
    # d=50
    {"codigo": "6010", "d": 50, "C_din": 22.9},
    {"codigo": "6210", "d": 50, "C_din": 37.1},
    {"codigo": "6310", "d": 50, "C_din": 65.0},
    # d=55
    {"codigo": "6011", "d": 55, "C_din": 29.6},
    {"codigo": "6211", "d": 55, "C_din": 46.2},
    {"codigo": "6311", "d": 55, "C_din": 74.1},
    # d=60
    {"codigo": "6012", "d": 60, "C_din": 31.9},
    {"codigo": "6212", "d": 60, "C_din": 55.3},
    {"codigo": "6312", "d": 60, "C_din": 85.2},
]

# Tenta carregar CSV se existir, senão usa o fallback
catalogo = catalogo_fallback
if os.path.exists("Inputs/rolamento_sfk.csv"):
    try:
        temp_cat = []
        with open("Inputs/rolamento_sfk.csv", newline="") as f:
            leitor = csv.DictReader(f)
            for row in leitor:
                temp_cat.append(
                    {
                        "codigo": row["codigo"],
                        "d": float(row["d"]),
                        "C_din": float(row["C_din"]),
                    }
                )
        if temp_cat:
            catalogo = temp_cat
    except:
        pass

# --- 3. EXECUÇÃO ---
d = ler_dados_para_mancais()
vida_util = 4000  # horas

# Monta lista de seleção
# (Nome, ForçaRadial, DiametroEixo, Rotacao)
lista_selecao = [
    ("Mancal 1 (Eixo 1 - A)", d["R1A"], d["d1"], d["n1"]),
    ("Mancal 2 (Eixo 1 - B)", d["R1B"], d["d1"], d["n1"]),
    ("Mancal 3 (Eixo 2 - A)", d["R2A"], d["d2"], d["n2"]),
    ("Mancal 4 (Eixo 2 - B)", d["R2B"], d["d2"], d["n2"]),
    ("Mancal 5 (Eixo 3 - A)", d["R3A"], d["d3"], d["n3"]),
    ("Mancal 6 (Eixo 3 - B)", d["R3B"], d["d3"], d["n3"]),
]

print("=== SELECAO DE ROLAMENTOS (SKF) ===")
print(f"Vida Util de Projeto: {vida_util} horas")

with open("Outputs/Resultados_Mancais.txt", "w") as f:
    f.write(f"DIMENSIONAMENTO DE MANCAIS (Vida: {vida_util}h)\n")
    f.write("=" * 60 + "\n\n")

    for nome, Fr, d_eixo, rpm in lista_selecao:
        f.write(f"{nome}:\n")
        f.write(
            f"  Condicoes: Fr={Fr:.1f} N | RPM={rpm:.1f} | d_eixo={d_eixo} mm\n"
        )

        # 1. Calcula C_req
        res = mancal.seleciona_rolamento(Fr, 0.0, rpm, vida_util)
        C_req = res["C_req_kN"]

        f.write(f"  Capacidade Dinamica Requerida (C_req): {C_req:.2f} kN\n")

        # 2. Busca no Catalogo
        rol = mancal.busca_catalogo(C_req, d_eixo, catalogo)

        if rol:
            f.write(f"  > SELECIONADO: Rolamento SKF {rol['codigo']}\n")
            f.write(f"    (d={rol['d']}mm, C_din={rol['C_din']} kN)\n")

            # Verificação extra de vida real
            vida_real = (rol["C_din"] * 1000 / res["P_proj_N"]) ** 3 * (
                1e6 / (60 * rpm)
            )
            f.write(f"    Vida Estimada: {vida_real:.0f} horas\n")
        else:
            f.write(
                f"  > FALHA: Nenhum rolamento serie 60/62/63 com d={d_eixo}mm atende a carga.\n"
            )
            f.write(
                "    Sugestao: Aumentar diametro do eixo ou usar rolamento de rolos.\n"
            )

        f.write("-" * 50 + "\n")

print("Calculo finalizado. Verifique 'Outputs/Resultados_Mancais.txt'.")
