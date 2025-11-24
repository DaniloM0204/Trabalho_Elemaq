import re
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
    plt.savefig(f'.graficos/{titulo.lower().replace(" ", "_")}_diagramas.png')
    return plt

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
    file.write(f"Fator dinamico (Kv): {res['K_v']:.3f}\n")
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
    Versão atualizada que lê todos os parâmetros do projeto incluindo
    materiais de engrenagens, eixos e comprimentos dos eixos.
    """
    parametros = {}

    mapeamento_chaves = {
        # Parâmetros principais do projeto
        "Potência de entrada do motor elétrico": "potencia_motor",
        "Velocidade de rotação do motor elétrico": "rotacao_motor",
        "Força no cabo necessária": "forca_cabo",
        "Diâmetro do tambor do guincho": "diametro_tambor",
        "Vida útil em número de horas": "vida_util",
        "Eficiência mínima do redutor": "eficiencia",

        # Material das engrenagens
        "S_at": "S_at",
        "S_ac": "S_ac",
        "C_p": "C_p",

        # Material dos eixos
        "S_ut": "S_ut",
        "S_y": "S_y",

        # Coeficientes
        "C_seg": "C_seg",
        "C_conf": "C_conf",
        "C_temp": "C_temp",
        "C_carr": "C_carr"
    }

    with open(nome_arquivo, 'r', encoding='utf-8') as arquivo:
        linhas = arquivo.readlines()

    # Processa cada linha
    for i, linha in enumerate(linhas):
        linha = linha.strip()

        # Ignora linhas vazias e comentários
        if not linha or linha.startswith('#'):
            continue

        # Processa comprimentos dos eixos (formato especial)
        if "comprimento eixo" in linha.lower():
            try:
                # Extrai o número do eixo e o comprimento
                if "eixo 1" in linha:
                    numero = extrair_numero(linha)
                    parametros['comprimento_eixo1'] = numero
                elif "eixo 2" in linha:
                    numero = extrair_numero(linha)
                    parametros['comprimento_eixo2'] = numero
                elif "eixo 3" in linha:
                    numero = extrair_numero(linha)
                    parametros['comprimento_eixo3'] = numero
            except Exception as e:
                print(f"Não foi possível extrair comprimento da linha: {linha}. Erro: {e}")
            continue

        # Processa linhas com separador : ou =
        separadores = [':', '=']
        separador_encontrado = None

        for sep in separadores:
            if sep in linha:
                separador_encontrado = sep
                break

        if separador_encontrado:
            partes = linha.split(separador_encontrado, 1)
            if len(partes) == 2:
                chave = partes[0].strip()
                valor = partes[1].strip()

                # Encontra a chave correspondente no mapeamento
                chave_correspondente = None
                for chave_mapeamento, chave_saida in mapeamento_chaves.items():
                    if chave_mapeamento in chave:
                        chave_correspondente = chave_saida
                        break

                # Se não encontrou no mapeamento, tenta usar a própria chave
                if not chave_correspondente:
                    chave_correspondente = chave.replace(' ', '_').lower()

                # Extrai o número
                numero = extrair_numero(valor)
                if numero is not None:
                    # Conversões especiais
                    if chave_correspondente == 'eficiencia' and numero > 1:
                        numero = numero / 100

                    parametros[chave_correspondente] = numero
                else:
                    # Se não conseguiu extrair número, armazena como string
                    parametros[chave_correspondente] = valor

    # Calcula valores derivados
    if 'S_ut' in parametros and 'Se_linha' not in parametros:
        parametros['Se_linha'] = 0.5 * parametros['S_ut']

    # Define valores padrão para coeficientes se não encontrados
    if 'C_seg' not in parametros:
        parametros['C_seg'] = 1.5
    if 'C_conf' not in parametros:
        parametros['C_conf'] = 0.814
    if 'C_temp' not in parametros:
        parametros['C_temp'] = 1.0
    if 'C_carr' not in parametros:
        parametros['C_carr'] = 1.0

    return parametros

def extrair_numero(texto):
    """
    Extrai um número de um texto, removendo unidades e convertendo vírgula para ponto
    """
    # Remove unidades comuns e espaços
    texto_limpo = re.sub(r'[N\.mm|N\.m|mm|N|MPa|kW|rpm|h]', '', texto).strip()

    # Encontra números (incluindo decimais com vírgula)
    padrao = r'[-+]?\d*[,.]?\d+'
    matches = re.findall(padrao, texto_limpo)

    if matches:
        # Pega o primeiro número encontrado e converte
        numero_str = matches[0].replace(',', '.')
        try:
            return float(numero_str)
        except ValueError:
            try:
                return int(numero_str)
            except ValueError:
                return None
    return None

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

def ler_resultados_diagramas(nome_arquivo='resultados_diagramas.txt'):
    """
    Lê um arquivo de texto com resultados dos diagramas e extrai os valores críticos
    de forma estruturada.

    Args:
        nome_arquivo (str): Caminho para o arquivo de texto

    Returns:
        dict: Dicionário estruturado com os valores críticos de cada eixo
    """
    resultados = {}

    try:
        with open(nome_arquivo, 'r', encoding='utf-8') as arquivo:
            linhas = arquivo.readlines()
    except FileNotFoundError:
        print(f"Erro: Arquivo {nome_arquivo} não encontrado.")
        return None

    eixo_atual = None

    for linha in linhas:
        linha = linha.strip()

        # Ignora linhas vazias e cabeçalhos
        if not linha or "RESULTADOS DOS DIAGRAMAS" in linha or "====" in linha:
            continue

        # Detecta um novo eixo
        if linha.upper().startswith('EIXO'):
            eixo_nome = linha.replace(':', '').strip().lower()
            eixo_atual = eixo_nome
            resultados[eixo_atual] = {}
            continue

        # Processa os dados do eixo atual
        if eixo_atual and ':' in linha:
            chave, valor = linha.split(':', 1)
            chave = chave.strip().lower()
            valor = valor.strip()

            # Extrai o número do valor (remove unidades se houver)
            numero = extrair_numero_diagrama(valor)

            if numero is not None:
                # Mapeia as chaves para nomes padronizados
                if 'momento máximo' in chave:
                    resultados[eixo_atual]['momento_maximo'] = numero
                elif 'cortante máximo' in chave:
                    resultados[eixo_atual]['cortante_maximo'] = numero
                elif 'reação a' in chave:
                    resultados[eixo_atual]['reacao_A'] = numero
                elif 'reação b' in chave:
                    resultados[eixo_atual]['reacao_B'] = numero
                elif 'posição m_max' in chave or 'posicao m_max' in chave:
                    resultados[eixo_atual]['posicao_M_max'] = numero

    return resultados

def extrair_numero_diagrama(texto):
    """
    Extrai um número de um texto, removendo unidades e convertendo vírgula para ponto
    """
    # Remove unidades comuns e espaços
    texto_limpo = re.sub(r'[N\.mm|N\.m|mm|N]', '', texto).strip()

    # Encontra números (incluindo decimais com vírgula)
    padrao = r'[-+]?\d*[,.]?\d+'
    matches = re.findall(padrao, texto_limpo)

    if matches:
        # Pega o primeiro número encontrado e converte
        numero_str = matches[0].replace(',', '.')
        try:
            return float(numero_str)
        except ValueError:
            try:
                return int(numero_str)
            except ValueError:
                return None
    return None

def ler_dimensionamento_eixos(nome_arquivo='dimensionamento_eixos.txt'):
    """
    Lê um arquivo de texto com resultados do dimensionamento dos eixos
    e extrai os valores de forma estruturada.

    Args:
        nome_arquivo (str): Caminho para o arquivo de texto

    Returns:
        dict: Dicionário estruturado com os resultados do dimensionamento
    """
    resultados = {
        'material': {},
        'fs_minimo': None,
        'eixos': {},
        'recomendacoes': []
    }

    try:
        with open(nome_arquivo, 'r', encoding='utf-8') as arquivo:
            linhas = arquivo.readlines()
    except FileNotFoundError:
        print(f"Erro: Arquivo {nome_arquivo} não encontrado.")
        return None

    eixo_atual = None
    secao_atual = None

    for linha in linhas:
        linha = linha.strip()

        # Ignora linhas vazias e cabeçalhos
        if not linha or "DIMENSIONAMENTO DOS EIXOS" in linha or "====" in linha:
            continue

        # Detecta seção de material
        if linha.startswith('Material:'):
            resultados['material']['nome'] = linha.split(':', 1)[1].strip()
            secao_atual = 'material'
            continue

        # Detecta propriedades do material
        if 'S_ut =' in linha and 'S_y =' in linha:
            # Extrai S_ut e S_y da mesma linha
            s_ut_match = re.search(r'S_ut\s*=\s*([\d.,]+)', linha)
            s_y_match = re.search(r'S_y\s*=\s*([\d.,]+)', linha)
            if s_ut_match:
                resultados['material']['S_ut'] = extrair_numero(s_ut_match.group(1))
            if s_y_match:
                resultados['material']['S_y'] = extrair_numero(s_y_match.group(1))
            continue

        # Detecta fator de segurança mínimo
        if 'Fator de segurança mínimo:' in linha:
            resultados['fs_minimo'] = extrair_numero(linha)
            continue

        # Detecta um novo eixo
        if linha.startswith('Eixo'):
            eixo_nome = linha.split(':', 1)[0].strip().lower()
            eixo_atual = eixo_nome
            resultados['eixos'][eixo_atual] = {}
            secao_atual = 'eixo'
            continue

        # Detecta seção de recomendações
        if linha.startswith('RECOMENDAÇÕES:'):
            secao_atual = 'recomendacoes'
            continue

        # Processa dados dos eixos
        if secao_atual == 'eixo' and eixo_atual and ':' in linha:
            chave, valor = linha.split(':', 1)
            chave = chave.strip().lower()
            valor = valor.strip()

            # Extrai o número ou texto
            if any(termo in chave for termo in ['diâmetro', 'fs', 'deflexão', 'velocidade']):
                numero = extrair_numero(valor)
                if numero is not None:
                    # Mapeia as chaves para nomes padronizados
                    if 'diâmetro' in chave:
                        resultados['eixos'][eixo_atual]['diametro_minimo'] = numero
                    elif 'fs fadiga' in chave:
                        resultados['eixos'][eixo_atual]['fs_fadiga'] = numero
                    elif 'fs escoamento' in chave:
                        resultados['eixos'][eixo_atual]['fs_escoamento'] = numero
                    elif 'deflexão' in chave and 'máxima' in chave:
                        resultados['eixos'][eixo_atual]['deflexao_maxima'] = numero
                    elif 'velocidade crítica' in chave:
                        resultados['eixos'][eixo_atual]['velocidade_critica'] = numero
            elif 'status' in chave:
                resultados['eixos'][eixo_atual]['status_deflexao'] = valor

        # Processa recomendações
        if secao_atual == 'recomendacoes' and linha.startswith('-'):
            resultados['recomendacoes'].append(linha[1:].strip())

    return resultados




