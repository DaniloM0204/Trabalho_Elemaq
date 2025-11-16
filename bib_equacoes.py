import numpy as np

constantes_engrenagem = {
        "S_at": 240, "S_ac": 860, "C_p": 191,
        "K_s": 1.0, "K_R": 1.25, "C_R": 1.25,
        "K_T": 1.0, "J": 0.32, "Qv": 6
    }


def geometria_engrenagem(m,N_p,i,b_face,φ=20):
    φ_rad = np.radians(φ)
    N_c = round(N_p * i)
    i_estagio = N_c / N_p
    d_p = m * N_p
    d_c = m * N_c
    r_p = d_p / 2
    r_c = d_c / 2
    b_face = m * b_face
    dist_centros = (d_p + d_c) / 2
    d_ext_p = m * (N_p + 2)
    d_ext_c = m * (N_c + 2)
    r_ext_p = d_ext_p / 2
    r_ext_c = d_ext_c / 2

    return {"m": m,"N_p": N_p,"N_c": N_c,"b_face": b_face,"d_p": d_p,  "d_c": d_c,"d_ext_p": d_ext_p,"d_ext_c": d_ext_c,"dist_centros": dist_centros,"i_estagio": i_estagio,"φ_rad": φ_rad,"r_p": r_p,
        "r_c": r_c,"r_ext_p": r_ext_p,"r_ext_c": r_ext_c,}


def razao_contato(m,r_p,r_c,r_ext_p,r_ext_c,dist_centros,φ_rad):
    p_circ = m * np.pi
    p_base = p_circ * np.cos(φ_rad)
    r_b_p = r_p * np.cos(φ_rad)
    r_b_c = r_c * np.cos(φ_rad)

    Z = (np.sqrt(r_ext_p**2 - r_b_p**2) +
        np.sqrt(r_ext_c**2 - r_b_c**2) -
         dist_centros * np.sin(φ_rad))

    # Razao de Contato
    mp = 0
    if p_base > 0:
        mp = Z / p_base

    return {"p_base": p_base,"Z": Z,"mp": mp,}

def cargas_e_velocidade(T_in,d_p,φ_rad,n_in):
    # Forca Tangencial
    W_t = (T_in * 1000) / (d_p / 2)
    # Forca Radial
    W_r = W_t * np.tan(φ_rad)

    # Velocidade tangencial
    ω_in = n_in * (2 * np.pi / 60)
    V_t = (d_p / 1000) * ω_in / 2

    return {"W_t": W_t,"W_r": W_r,"V_t": V_t,}

def saida_estagio(T_in,i_estagio,eta_estagio,n_in):
    T_out = T_in * i_estagio * eta_estagio
    n_out = n_in / i_estagio

    return {"T_out": T_out,"n_out": n_out,}

def faotr_dinamico(Qv,V_t):
    B = np.pow(12 - Qv, 2/3) / 4
    A = 50 + 56 * (1 - B)
    denominador = A + np.sqrt(V_t * 200)
    if denominador == 0:
        return 1.0
    K_v = (A / denominador)**B
    return K_v

def fatorKm(b_face):
    if b_face <= 50:
        return 1.6
    elif 50 < b_face <= 150:
        return 1.7
    elif 250 < b_face <= 500:
        return 2.0
    else:
        return 2.2

def flexao(W_t, b_face, m, J, K_s, K_m, K_v, S_at, K_T, K_R):
    σ_b_adm = S_at / (K_T * K_R)

    denominador = (b_face * m * J * K_v)
    if denominador == 0:
        σ_b = np.inf
    else:
        σ_b = (W_t / denominador) * (K_s * K_m)

    FS_flexao = np.inf if σ_b == 0 else σ_b_adm / σ_b
    return {"σ_b": σ_b, "σ_b_adm": σ_b_adm, "FS_flexao": FS_flexao}

def pitting(N_c, N_p, φ_rad, W_t, K_s, K_m, b_face, d_p, K_v, C_p, S_ac, C_R):
    σ_c_adm = S_ac / C_R
    mg = N_c / N_p

    I_geo = (np.sin(φ_rad) * np.cos(φ_rad) / 2) * (mg / (mg + 1))

    denominador = (b_face * d_p * I_geo * K_v)

    σ_c = np.inf
    if denominador != 0:
        pitting = (W_t * K_s * K_m) / denominador
        if pitting >= 0:
            σ_c = C_p * np.sqrt(pitting)

    FS_pitting = 0.0
    if σ_c == np.inf:
        FS_pitting = 0.0
    elif σ_c == 0:
        FS_pitting = np.inf
    else:
        FS_pitting = (σ_c_adm / σ_c)**2

    return {"σ_c": σ_c, "σ_c_adm": σ_c_adm, "I_geo": I_geo, "FS_pitting": FS_pitting}
