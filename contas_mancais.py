import equacoes_mancais as mancal
import contas_eixo as eixo
import csv
import os

# Definicao manual das rotacoes para evitar erro de leitura de arquivo com acento
# (Copie os valores do seu arquivo de engrenagens se forem diferentes)
n1 = 1450.0
n2 = 395.5
n3 = 107.9

vida_util = 4000

# Catalogo Caso CSV falhe
catalogo = [
    {'codigo': '6004', 'd': 20, 'C_din': 9.95},
    {'codigo': '6005', 'd': 25, 'C_din': 11.9},
    {'codigo': '6006', 'd': 30, 'C_din': 13.3},
    {'codigo': '6205', 'd': 25, 'C_din': 14.0},
    {'codigo': '6206', 'd': 30, 'C_din': 19.5},
    {'codigo': '6306', 'd': 30, 'C_din': 28.1}
]

# Tenta ler CSV
if os.path.exists('Inputs/rolamento_sfk.csv'):
    try:
        with open('Inputs/rolamento_sfk.csv', newline='') as f: # Sem encoding especifico
            leitor = csv.DictReader(f)
            catalogo_csv = []
            for row in leitor:
                catalogo_csv.append({'codigo': row['codigo'], 'd': float(row['d']), 'C_din': float(row['C_din'])})
            if catalogo_csv: catalogo = catalogo_csv
    except: pass

# Puxa os diametros COMERCIAIS calculados no contas_eixo.py
d1 = eixo.d1_com
d2 = eixo.d2_com
d3 = eixo.d3_com

# Puxa as reacoes
R1 = eixo.reações_eixo1
R2 = eixo.reações_eixo2
R3 = eixo.reações_eixo3

lista_mancais = [
    ('Mancal 1 (Eixo 1 - A)', R1['R_A_resultante'], d1, n1),
    ('Mancal 2 (Eixo 1 - B)', R1['R_B_resultante'], d1, n1),
    ('Mancal 3 (Eixo 2 - A)', R2['R_A_resultante'], d2, n2),
    ('Mancal 4 (Eixo 2 - B)', R2['R_B_resultante'], d2, n2),
    ('Mancal 5 (Eixo 3 - A)', R3['R_A_resultante'], d3, n3),
    ('Mancal 6 (Eixo 3 - B)', R3['R_B_resultante'], d3, n3),
]

with open('Outputs/Resultados_Mancais.txt', 'w') as f:
    f.write("Dimensionamento de mancais\n")

    for (nome, Fr, d_min, rpm) in lista_mancais:
        sel = mancal.seleciona_rolamento(Fr, 0.0, rpm, vida_util)
        C_req_kN = sel['C_req_kN']

        rol = mancal.busca_catalogo(C_req_kN, d_min, catalogo)

        f.write(f"{nome}:\n")
        f.write(f"  Carga: {Fr:.1f} N, RPM: {rpm:.1f}, D_min: {d_min:.1f} mm\n")
        f.write(f"  C_req: {C_req_kN:.2f} kN\n")

        if rol:
            f.write(f"  > SELECIONADO: SKF {rol['codigo']} (C={rol['C_din']} kN)\n")
        else:
            f.write("  > FALHA: Nenhum rolamento adequado.\n")
        f.write("-" * 50 + "\n")
