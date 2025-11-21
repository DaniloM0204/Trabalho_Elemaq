
import re
import numpy as np
import matplotlib.pyplot as plt

def plot_diagramas_cortante_momento(x, V, M, titulo):
    """
    Função genérica para plotar diagramas
    """
    plt.figure(figsize=(10, 6))
    
    plt.subplot(2, 1, 1)
    plt.plot(x, V, 'b-', linewidth=2)
    plt.grid(True, alpha=0.3)
    plt.title(f'{titulo} - Esforço Cortante')
    plt.ylabel('Força (N)')
    
    plt.subplot(2, 1, 2)
    plt.plot(x, M, 'r-', linewidth=2)
    plt.grid(True, alpha=0.3)
    plt.title(f'{titulo} - Momento Fletor')
    plt.xlabel('Posição (mm)')
    plt.ylabel('Momento (N.mm)')
    
    plt.tight_layout()
    return plt

def converte_unidades(valor, unidade_origem, unidade_destino):
    """
    Conversor de unidades comum
    """
    conversoes = {
        ('N.m', 'N.mm'): lambda x: x * 1000,
        ('N.mm', 'N.m'): lambda x: x / 1000,
        ('mm', 'm'): lambda x: x / 1000,
        ('m', 'mm'): lambda x: x * 1000,
        ('RPM', 'rad/s'): lambda x: x * 2 * np.pi / 60,
        ('rad/s', 'RPM'): lambda x: x * 60 / (2 * np.pi)
    }
    
    return conversoes.get((unidade_origem, unidade_destino), lambda x: x)(valor)

def valida_entradas(**kwargs):

    """
    Valida entradas do usuário
    """
    for key, value in kwargs.items():
        if value <= 0:
            raise ValueError(f"Valor de {key} deve ser positivo")
    return True

def escreve_estagio(file, nome_estagio, res):
    """
    Escreve resultados incluindo informações de vida útil
    """
    if res is None:
        file.write("Nenhum modulo padrao foi suficiente.\n")
        return

    file.write(f"\n--- {nome_estagio}: {res['status_completo']} ---\n")
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
    file.write("  Vida Util Estimada:\n")
    file.write(f"    Vida por Flexao: {res['vida_flexao_horas']:.0f} horas\n")
    file.write(f"    Vida por Pitting: {res['vida_pitting_horas']:.0f} horas\n")
    file.write(f"    Vida Minima: {res['vida_util_minima']:.0f} horas\n")
    file.write(f"    Status Vida: {res['status_vida']}\n")


def ler_Dados_De_Entrada(nome_arquivo):
    """
    Versão robusta que garante que todos os valores são números
    """
    parametros = {}
    
    mapeamento_chaves = {
        "Potência de entrada do motor elétrico": "potencia_motor",
        "Velocidade de rotação do motor elétrico": "rotacao_motor", 
        "Força no cabo necessária": "forca_cabo",
        "Diâmetro do tambor do guincho": "diametro_tambor",
        "Vida útil em número de horas": "vida_util",
        "Eficiência mínima do redutor": "eficiencia",
        "S_at": "S_at",
        "S_ac": "S_ac", 
        "C_p": "C_p"
    }
    
    with open(nome_arquivo, 'r', encoding='utf-8') as arquivo:
        for linha in arquivo:
            linha = linha.strip()
            
            if not linha or linha.startswith('#'):
                continue
                
            if ':' in linha:
                chave, valor = linha.split(':', 1)
                chave = chave.strip()
                valor = valor.strip()
                
                # Encontra a chave correspondente no mapeamento
                chave_correspondente = None
                for chave_mapeamento, chave_saida in mapeamento_chaves.items():
                    if chave_mapeamento in chave:
                        chave_correspondente = chave_saida
                        break
                
                if chave_correspondente:
                    # Extrai o primeiro número encontrado
                    import re
                    numeros = re.findall(r'[-+]?\d+[,.]?\d*', valor)
                    if numeros:
                        numero_str = numeros[0].replace(',', '.')
                        try:
                            if '.' in numero_str:
                                numero = float(numero_str)
                            else:
                                numero = int(numero_str)
                            
                            # Conversão especial para eficiência
                            if chave_correspondente == 'eficiencia' and numero > 1:
                                numero = numero / 100
                                
                            parametros[chave_correspondente] = numero
                            
                        except ValueError:
                            print(f"AVISO: Não pude converter '{numero_str}' para número")
    
    return parametros


def ler_Estagios_Engrenagem(nome_arquivo):
    """
    Lê um arquivo de texto com resultados do dimensionamento do redutor
    e extrai todos os parâmetros de forma estruturada.
    
    Args:
        nome_arquivo (str): Caminho para o arquivo de texto
        
    Returns:
        dict: Dicionário estruturado com todos os resultados
    """
    resultados = {
        'torques': {},
        'estagio1': {},
        'estagio2': {},
        'transmissao_total': {},
        'parametros_eixos': {}
    }
    
    try:
        with open(nome_arquivo, 'r', encoding='utf-8') as arquivo:
            linhas = arquivo.readlines()
    except FileNotFoundError:
        print(f"Erro: Arquivo {nome_arquivo} não encontrado.")
        return None
    
    estagio_atual = None
    secao_atual = None
    
    for linha in linhas:
        linha = linha.strip()
        
        # Ignora linhas vazias
        if not linha:
            continue
        
        # Detecta seções principais
        if linha.startswith('Torque de Entrada:'):
            valor = extrair_numero(linha)
            resultados['torques']['entrada'] = valor
            
        elif linha.startswith('Torque de Saida:'):
            valor = extrair_numero(linha)
            resultados['torques']['saida'] = valor
            
        elif '--- Estagio 1' in linha:
            estagio_atual = 'estagio1'
            
        elif '--- Estagio 2' in linha:
            estagio_atual = 'estagio2'
            
        elif linha.startswith('Razao de Transmissao Total Efetiva:'):
            valor = extrair_numero(linha)
            resultados['transmissao_total']['efetiva'] = valor
            
        elif linha.startswith('Razao de Transmissao Total Desejada:'):
            valor = extrair_numero(linha)
            resultados['transmissao_total']['desejada'] = valor
            
        elif linha.startswith('Erro Percentual:'):
            valor = extrair_numero(linha)
            resultados['transmissao_total']['erro_percentual'] = valor
            
        elif linha.startswith('Parametros dos Eixos:'):
            estagio_atual = 'parametros_eixos'
            secao_atual = None
            
        # Processa dados dos estágios
        elif estagio_atual in ['estagio1', 'estagio2']:
            if 'Modulo (m):' in linha:
                resultados[estagio_atual]['modulo'] = extrair_numero(linha)
            elif 'Largura da Face (b):' in linha:
                resultados[estagio_atual]['largura_face'] = extrair_numero(linha)
            elif 'Angulo de Pressao:' in linha:
                resultados[estagio_atual]['angulo_pressao'] = extrair_numero(linha)
            elif 'Razao de Contato (mp):' in linha:
                resultados[estagio_atual]['razao_contato'] = extrair_numero(linha)
            elif 'N de Dentes (Np):' in linha and 'Pinhao:' not in linha:
                resultados[estagio_atual]['pinhao_dentes'] = extrair_numero(linha)
            elif 'Diametro Primitivo (dp):' in linha:
                resultados[estagio_atual]['pinhao_diametro_primitivo'] = extrair_numero(linha)
            elif 'Diametro Externo (dep):' in linha:
                resultados[estagio_atual]['pinhao_diametro_externo'] = extrair_numero(linha)
            elif 'N de Dentes (Nc):' in linha and 'Coroa:' not in linha:
                resultados[estagio_atual]['coroa_dentes'] = extrair_numero(linha)
            elif 'Diametro Primitivo (dc):' in linha:
                resultados[estagio_atual]['coroa_diametro_primitivo'] = extrair_numero(linha)
            elif 'Diametro Externo (dec):' in linha:
                resultados[estagio_atual]['coroa_diametro_externo'] = extrair_numero(linha)
            elif 'Razao de Transmissao (i):' in linha:
                resultados[estagio_atual]['razao_transmissao'] = extrair_numero(linha)
            elif 'Distancia entre Centros (C):' in linha:
                resultados[estagio_atual]['distancia_centros'] = extrair_numero(linha)
            elif 'FS (Flexao):' in linha:
                resultados[estagio_atual]['fs_flexao'] = extrair_numero(linha)
            elif 'FS (Superficie):' in linha:
                resultados[estagio_atual]['fs_superficie'] = extrair_numero(linha)
            elif 'Vida por Flexao:' in linha:
                resultados[estagio_atual]['vida_flexao'] = extrair_numero(linha)
            elif 'Vida por Pitting:' in linha:
                resultados[estagio_atual]['vida_pitting'] = extrair_numero(linha)
            elif 'Vida Minima:' in linha:
                resultados[estagio_atual]['vida_minima'] = extrair_numero(linha)
            elif 'Status Vida:' in linha:
                resultados[estagio_atual]['status_vida'] = linha.split(':')[-1].strip()
                
        # Processa parâmetros dos eixos
        elif estagio_atual == 'parametros_eixos':
            if 'Rotacao (n_eixo1):' in linha:
                resultados['parametros_eixos']['rotacao_eixo1'] = extrair_numero(linha)
            elif 'Torque (T_eixo1):' in linha:
                resultados['parametros_eixos']['torque_eixo1'] = extrair_numero(linha)
            elif 'Forca Tangencial (W_t1):' in linha:
                resultados['parametros_eixos']['forca_tangencial_pinhao1'] = extrair_numero(linha)
            elif 'Forca Radial (W_r1):' in linha:
                resultados['parametros_eixos']['forca_radial_pinhao1'] = extrair_numero(linha)
            elif 'Rotacao (n_eixo2):' in linha:
                resultados['parametros_eixos']['rotacao_eixo2'] = extrair_numero(linha)
            elif 'Torque (T_eixo2):' in linha:
                resultados['parametros_eixos']['torque_eixo2'] = extrair_numero(linha)
            elif 'Forca Tangencial (W_t_c1):' in linha:
                resultados['parametros_eixos']['forca_tangencial_coroa1'] = extrair_numero(linha)
            elif 'Forca Radial (W_r_c1):' in linha:
                resultados['parametros_eixos']['forca_radial_coroa1'] = extrair_numero(linha)
            elif 'Forca Tangencial (W_t_p2):' in linha:
                resultados['parametros_eixos']['forca_tangencial_pinhao2'] = extrair_numero(linha)
            elif 'Forca Radial (W_r_p2):' in linha:
                resultados['parametros_eixos']['forca_radial_pinhao2'] = extrair_numero(linha)
            elif 'Rotacao (n_eixo3):' in linha:
                resultados['parametros_eixos']['rotacao_eixo3'] = extrair_numero(linha)
            elif 'Torque (T_eixo3):' in linha:
                resultados['parametros_eixos']['torque_eixo3'] = extrair_numero(linha)
            elif 'Forca Tangencial (W_t_c2):' in linha:
                resultados['parametros_eixos']['forca_tangencial_coroa2'] = extrair_numero(linha)
            elif 'Forca Radial (W_r_c2):' in linha:
                resultados['parametros_eixos']['forca_radial_coroa2'] = extrair_numero(linha)
    
    return resultados

def extrair_numero(texto):
    """
    Extrai um número de um texto, removendo unidades e convertendo vírgula para ponto
    """
    # Encontra todos os números no texto (incluindo decimais com vírgula)
    padrao = r'[-+]?\d*[,.]?\d+'
    matches = re.findall(padrao, texto)
    
    if matches:
        # Pega o primeiro número encontrado e converte
        numero_str = matches[0].replace(',', '.')
        try:
            # Tenta converter para float primeiro
            return float(numero_str)
        except ValueError:
            # Se falhar, tenta int
            try:
                return int(numero_str)
            except ValueError:
                return None
    return None

# Função auxiliar para exibir os resultados de forma organizada
def exibir_resultados(resultados):
    """
    Exibe os resultados de forma organizada para verificação
    """
    if not resultados:
        print("Nenhum resultado para exibir")
        return
    
    print("=== RESULTADOS DO REDUTOR ===")
    
    # Torques
    print("\n--- TORQUES ---")
    for chave, valor in resultados['torques'].items():
        print(f"{chave}: {valor} N.m")
    
    # Estágio 1
    print("\n--- ESTÁGIO 1 ---")
    for chave, valor in resultados['estagio1'].items():
        print(f"{chave}: {valor}")
    
    # Estágio 2
    print("\n--- ESTÁGIO 2 ---")
    for chave, valor in resultados['estagio2'].items():
        print(f"{chave}: {valor}")
    
    # Transmissão total
    print("\n--- TRANSMISSÃO TOTAL ---")
    for chave, valor in resultados['transmissao_total'].items():
        print(f"{chave}: {valor}")
    
    # Parâmetros dos eixos
    print("\n--- PARÂMETROS DOS EIXOS ---")
    for chave, valor in resultados['parametros_eixos'].items():
        print(f"{chave}: {valor}")
