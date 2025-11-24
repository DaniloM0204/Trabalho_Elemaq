from __future__ import print_function
import equacoes_cubo_eixo as ce
import os
import math


def mm(x):
    return x / 1000.0

# (fornecidos manualmente para evitar dependência em leitura de arquivos)
d1 = 30.0
d2 = 35.0
d3 = 40.0

# diâmetros primitivos das engrenagens (mm) — valores usados no relatório
dp_pinh1 = 36.0
dp_coroa1 = 132.0
dp_pinh2 = 45.0
dp_coroa2 = 165.0

# rugosidade assumida (m)
R1 = 6e-6
R2 = 6e-6

# Observação: a pressão mínima de assentamento (`p_min`) é calculada por
# acoplamento usando a equação tangencial (p_min_tangencial). Não é
# necessária uma constante global predefinida.

# Pressões máximas limitantes (valores obtidos a partir da análise - ver Tabela da imagem):
# p_max_eixo = 0.35 GPa para todos os eixos (350 MPa)
# p_max_cubo varia por acoplamento (extraído da Tabela 13 da imagem)
p_max_eixo = 0.35e9   # 0.35 GPa

# p_max do cubo por acoplamento (GPa -> Pa)
p_max_cubo_map = {
    'Pinhão1': 0.16e9,  # 0.16 GPa
    'Coroa1':  0.18e9,  # 0.18 GPa
    'Pinhão2': 0.15e9,  # 0.15 GPa
    'Coroa2':  0.18e9   # 0.18 GPa
}

# Parâmetros para cálculo de p_min (Equação 4.1 da imagem)
# p_min_tangencial = 2 * S_F * T_tangencial / (pi * d^2 * l * mu_et)
# Onde T_tangencial é a força tangencial W_t (N), d é o diâmetro do cubo (m),
# l é o comprimento/face do cubo (m), mu_et é o coeficiente de atrito tangencial,
# S_F é o fator de segurança adotado.
S_F = 1.5
mu_et = 0.2  # coeficiente de atrito tangencial (valor assumido; ajustar conforme dados experimentais)

# Propriedades do material para cálculo de p_max
# p_max_cubo = 0.45 * sigma_e * (r2^2 - r^2) / r2^2
# p_max_eixo = 0.9 * sigma_e
sigma_e = 390e6  # tensão de escoamento AISI 1020 (Pa)

couplings = [
    ("Eixo1 - Pinhão1", d1, dp_pinh1, 'Pinhão1'),
    ("Eixo2 - Coroa1", d2, dp_coroa1, 'Coroa1'),
    ("Eixo2 - Pinhão2", d2, dp_pinh2, 'Pinhão2'),
    ("Eixo3 - Coroa2", d3, dp_coroa2, 'Coroa2'),
]

print("Calculando sobremedida/interferência para acoplamentos (valores em μm):\n")
results = []
for name, shaft_d_mm, hub_d_mm, tag in couplings:
    D1 = mm(shaft_d_mm)
    D2 = mm(hub_d_mm)
    # definimos D nominal como média entre eixo e cubo (aproximação prática)
    D_nom = (D1 + D2) / 2.0

    # --- calcular p_min tangencial usando força tangencial W_t ---
    # força tangencial W_t = 2 * Torque / Dp (T in N.m, Dp in m) -> result in N
    # Para cada acoplamento associamos o torque do eixo: T1, T2, T3 (N.m)
    T1 = 10.41
    T2 = 35.18
    T3 = 118.91
    # escolher torque conforme eixo do acoplamento (Eixo1 uses T1, Eixo2 uses T2, Eixo3 uses T3)
    if name.startswith('Eixo1'):
        T_axis = T1
    elif name.startswith('Eixo2'):
        T_axis = T2
    else:
        T_axis = T3

    # força tangencial (N)
    W_t = (2.0 * T_axis) / D2

    # comprimento do cubo: usar largura de face por estágio (b1=10mm, b2=20mm)
    # associar Pinhão/ Coroa do estágio 1 -> b1; estágio 2 -> b2
    b1 = 0.010  # m
    b2 = 0.020  # m
    # determinar comprimento l: pinhões/coroas do estágio 1 usam b1, do estágio 2 usam b2
    if tag in ('Pinhão1', 'Coroa1'):
        l = b1
    else:
        l = b2

    # diámetro do cubo d (usar D2)
    d_cubo = D2

    # calcular p_min_tangencial (Pa)
    p_min_calc = (2.0 * S_F * W_t) / (math.pi * d_cubo**2 * l * mu_et)

    # --- calcular p_maxes via propriedades (Equações 4.2 e 4.3) ---
    r = D1 / 2.0
    r2 = D2 / 2.0
    # proteger contra divisão por zero
    if r2 <= 0:
        p_max_cubo_calc = float('inf')
    else:
        p_max_cubo_calc = 0.45 * sigma_e * ((r2**2 - r**2) / (r2**2))

    p_max_eixo_calc = 0.9 * sigma_e

    # valor limite adotado
    p_max_lim = min(p_max_eixo_calc, p_max_cubo_calc)

    # usar p_min_calc como p_min for interference calc
    res_min = ce.interferencia_total(D_nom, D1, D2, p_min_calc, R1, R2)
    res_max = ce.interferencia_total(D_nom, D1, D2, p_max_lim, R1, R2)

    delta_dr_um = res_min['Delta_dr_m'] * 1e6
    delta_total_min_um = res_min['Interferencia_total_m'] * 1e6
    delta_total_max_um = res_max['Interferencia_total_m'] * 1e6

    # guardar também p_min_calc e p_max_lim em MPa para saída
    p_min_used_MPa = p_min_calc / 1e6
    p_max_used_MPa = p_max_lim / 1e6

    results.append((name, delta_dr_um, delta_total_min_um, delta_total_max_um, p_min_used_MPa, p_max_used_MPa))

    for r in results:
        name, delta_dr, min_total, max_total, pminm, pmaxm = r
        print(f"{name}: Interf. mínima (rugosidade) = {delta_dr:.2f} µm | Interf. mínima total (p_min) = {min_total:.2f} µm | Interf. máxima (p_max) = {max_total:.2f} µm | p_min={pminm:.2f} MPa | p_max_lim={pmaxm:.2f} MPa")

# também salva em Outputs
outdir = 'Outputs'
if not os.path.exists(outdir):
    os.makedirs(outdir)

with open(os.path.join(outdir, 'interferencia_acoplamentos.txt'), 'w') as f:
    f.write('Acoplamento;Delta_dr_um;Interf_min_total_um;Interf_max_total_um;p_min_MPa;p_max_lim_MPa\n')
    for (name, delta_dr, min_total, max_total, pminm, pmaxm) in results:
        f.write(f"{name};{delta_dr:.4f};{min_total:.4f};{max_total:.4f};{pminm:.3f};{pmaxm:.3f}\n")

print('\nResultados gravados em Outputs/interferencia_acoplamentos.txt')
