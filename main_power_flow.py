from pathlib import Path
import numpy as np
import pandas as pd
import cmath
from leitura_pwf import *



# === Funções auxiliares para cálculo de Fluxo de Potência === #
def calcula_ybarra(nbus, dlin, shunt):
    Y = np.zeros((nbus, nbus), dtype=complex)
    for _, row in dlin.iterrows():
        i, j = int(row.de - 1), int(row.para - 1)
        y = 1 / complex(row.r, row.x)
        tap = row.tap
        defas = row.defasagem

        Y[i, i] += (1 / tap**2) * y + 1j * row.bsh
        Y[j, j] += y + 1j * row.bsh
        Y[i, j] -= y / tap * cmath.exp(-1j * defas)
        Y[j, i] -= y / tap * cmath.exp(1j * defas)

    for k in range(nbus):
        Y[k, k] += 1j * shunt[k]

    return Y

def calcula_residuos(v, teta, ybarra, tipo, pesp, qesp, dbar, CREM):

    # Calcula o vetor de tensões na forma retangular
    x, y = np.cos(teta) * v, np.sin(teta) * v
    vret = x + 1j * y
    
    # Calcula o vetor de correntes injetadas
    i_inj = ybarra @ vret
    
    # Calcula o vetor de potências (S = V x I*)    
    s = vret * np.conj(i_inj)
    pcalc = np.real(s)
    qcalc = np.imag(s)

    # Calcula os resíduos
    delta_p = pesp - pcalc
    delta_q = qesp - qcalc

    # Zera os resíduos de P para swing e Q paar swing e PV
    delta_p[tipo == 2] = 0  # barra swing
    delta_q[tipo == 2] = 0  # barra swing
    delta_q[tipo == 1] = 0  # barra PV

    delta_y = np.concatenate([delta_p, delta_q])

    if CREM:
        delta_v_completo = dbar.v.values - v
        delta_v = delta_v_completo[tipo == 4]

        delta_y = np.hstack([delta_y, delta_v])


    return delta_y, pcalc, qcalc

def calcula_jacobiano(v, teta, pcalc, qcalc, ybarra, tipo, barra, bc, CREM=False):
    nbar = len(v)
    g = np.real(ybarra)
    b = np.imag(ybarra)

    H = np.zeros((nbar, nbar))
    N = np.zeros((nbar, nbar))
    M = np.zeros((nbar, nbar))
    L = np.zeros((nbar, nbar))

    for k in range(nbar):
        for m in range(nbar):
            if k == m:
                H[k, k] = -(qcalc[k] + v[k]**2 * b[k, k])
                N[k, k] = (pcalc[k] + v[k]**2 * g[k, k]) / v[k]
                M[k, k] = pcalc[k] - v[k]**2 * g[k, k]
                L[k, k] = (qcalc[k] - v[k]**2 * b[k, k]) / v[k]

                if tipo[k] == 2:
                    H[k, k] = 1e10
                    L[k, k] = 1e10
                elif tipo[k] == 1:
                    L[k, k] = 1e10
            else:
                delta = teta[k] - teta[m]
                H[k, m] = v[k]*v[m]*(g[k, m]*np.sin(delta) - b[k, m]*np.cos(delta))
                N[k, m] = v[k]*(g[k, m]*np.cos(delta) + b[k, m]*np.sin(delta))
                M[k, m] = -v[k]*v[m]*(g[k, m]*np.cos(delta) + b[k, m]*np.sin(delta))
                L[k, m] = v[k]*(g[k, m]*np.sin(delta) - b[k, m]*np.cos(delta))

    Jac = np.block([[H, N], [M, L]])
        
    if CREM:
        idx_tipo3 = np.where(tipo == 3)[0]
        nbar_tipo3 = idx_tipo3.size

        # Cria linhas adicionais
        M_CREM_linha = np.zeros((nbar_tipo3, nbar))
        L_CREM_linha = np.zeros((nbar_tipo3, nbar))

        for k in range(nbar_tipo3):
            barra_controlada = bc[idx_tipo3[k]]       
            idx_coluna = np.where(barra == barra_controlada)[0][0] 
            L_CREM_linha[k, idx_coluna] = 1

        linha_adicional_CREM = np.hstack([M_CREM_linha, L_CREM_linha])

        # Cria colunas adicionais
        N_CREM_coluna = np.zeros((nbar, nbar_tipo3))
        L_CREM_coluna = np.zeros((nbar, nbar_tipo3))

        for k in range(nbar_tipo3):
            barra_controladora = barra[idx_tipo3[k]]                 
            idx_linha = np.where(barra == barra_controladora)[0][0]
            L_CREM_coluna[idx_linha, k] = -1
        
        coluna_adicional_CREM = np.vstack([N_CREM_coluna, L_CREM_coluna])     
        
        # Cria matriz adicional para CREM (derivadas das tensões em relação à potência reativa gerada ~ 0)
        matriz_adicional_CREM = np.zeros((nbar_tipo3, nbar_tipo3))

        Jac_linhas_adicionais = np.vstack([Jac, linha_adicional_CREM])
        coluna_adicional_com_matriz_adicional = np.vstack([coluna_adicional_CREM, matriz_adicional_CREM])
        Jac_colunas_adicionais = np.hstack([Jac_linhas_adicionais, coluna_adicional_com_matriz_adicional])
        Jac = Jac_colunas_adicionais

    return Jac


# === Função para cálculo de Fluxo de Potência === #
def calcula_fluxo_de_potencia_newt(dbar, dlin, tol=1e-4, iter_max=25, flat_start=True, QLIM=False, CREM=False):
    
    # Variáveis gerais
    nbar = len(dbar)
    barra = dbar.barra.values
    bc = dbar.bc.values
    shunt = dbar.shunt.values
    tipo = dbar.tipo.values.copy()

    # Variáveis de estado
    teta = dbar.teta.values.copy()
    v = dbar.v.values.copy()
    
    # Variáveis de potência
    pg, pl = dbar.pg.values, dbar.pl.values
    qg, ql = dbar.qg.values, dbar.ql.values
    pesp, qesp = pg - pl, qg - ql
    qmin, qmax = dbar.qn.values - ql, dbar.qm.values - ql

    # Cálculo da matriz Ybarra
    ybarra = calcula_ybarra(nbar, dlin, shunt)

    # A barra PV que controla passa a ser P (tipo 3) e barra PQ controlada passa a ser PQV (tipo 4)
    if CREM:
        # Atualizar barras tipo 1 com bc ≠ 0 para tipo 3 (P)
        tipo[(tipo == 1) & (bc != 0) & (bc != barra)] = 3

        # Atualizar barras cujo número está no campo 'bc' de outras barras para tipo 4 (PQV)
        tipo[np.isin(barra, bc[(tipo == 3)])] = 4

        idx_tipo3 = np.where(tipo == 3)[0]
        nbar_tipo3 = idx_tipo3.size


    # Zera ângulos das baras PV e PQ e define como 1,0 pu tensões de barras PQ
    if flat_start:
        teta[tipo == 1] = 0
        teta[tipo == 0] = 0
        v[tipo == 0] = 1.0

    tipo_pre_iteracoes = tipo.copy()

    # Iterações
    for it in range(iter_max):
        delta_residuos, pcalc, qcalc = calcula_residuos(v, teta, ybarra, tipo, pesp, qesp, dbar, CREM)
        
         # Se o QLIM estiver ativado e houver alteração no Q calculado, verifica se a barra pode voltar a ser PV ou não
        # Essa lógica esta incorreta pq a barra foi alterada para tipo 0 lá em cima
        if QLIM and np.any((qesp == qmin) | (qesp == qmax)):
            cond_tensao_tipo_1 = (
                ((v < dbar.v.values) & (qesp == qmin)) |
                ((v > dbar.v.values) & (qesp == qmax))
            )
            
            cond_barras_tipo_1 = (tipo == 0) & (tipo_pre_iteracoes == 1) & cond_tensao_tipo_1
            tipo[cond_barras_tipo_1] = 1

            barras_retornaram_pv = barra[cond_barras_tipo_1]
            print(f"Barras que voltaram a ser PV: {barras_retornaram_pv.tolist()}")

            # Identifica barras que eram tipo 3 e viraram tipo 0
            cond_barras_tipo_3 = (tipo == 0) & (tipo_pre_iteracoes == 3)

            # Pega os códigos das barras controladoras dessas barras
            # barras_controladoras = bc[cond_barras_tipo_3]

            # # Converte vetor 'barra' em índice: encontra onde cada controladora está no vetor 'barra'
            # # Cria um mapeamento: código da barra → índice
            # mapa_barra_para_indice = dict(zip(barra, np.arange(len(barra))))
            # indices_controladoras = np.vectorize(mapa_barra_para_indice.get)(barras_controladoras)

            # # Coleta os valores nas posições certas
            # v_ctrl = v[indices_controladoras]
            # v_ref = dbar.v.values[indices_controladoras]
            # qmin_ctrl = qmin[cond_barras_tipo_3]
            # qmax_ctrl = qmax[cond_barras_tipo_3]
            # qcalc_ctrl = qcalc[cond_barras_tipo_3]

            # # Verifica condição de tensão das controladoras
            # cond_tensao_ctrl = (
            #     ((v_ctrl < v_ref) & (qcalc_ctrl == qmin_ctrl)) |
            #     ((v_ctrl > v_ref) & (qcalc_ctrl == qmax_ctrl))
            # )

            # # Atualiza tipo = 3 nas barras tipo 3 envolvidas
            # idx_tipo3 = np.where(cond_barras_tipo_3)[0]
            # tipo[idx_tipo3[cond_tensao_ctrl]] = 3

            # # Atualiza tipo = 4 nas barras controladoras envolvidas
            # tipo[indices_controladoras[cond_tensao_ctrl]] = 4
            
            # v[cond] = dbar.v.values[cond]



        if QLIM and np.any(((qcalc > qmax) | (qcalc < qmin)) & np.isin(tipo, [1, 3])):
            cond_tipo = (tipo == 1) | (tipo == 3)
            violou_sup = (qcalc > qmax) & cond_tipo 
            violou_inf = (qcalc < qmin) & cond_tipo
            violou_limite = violou_sup | violou_inf

            print(f"Iteração {it:>3}: Q calculado(s) além dos limites informados.")
            print(f"Barras violadas: {barra[violou_limite]}")

            # Identifica índices das barras tipo 3 que violaram limite
            idx_tipo3_violado = np.where((violou_limite) & (tipo == 3))[0]

            # Atualizar tipo para 0 onde houve violação
            tipo = np.where(violou_limite, 0, tipo)
           
            # Barras controladas por essas barras tipo 3
            barras_controladas = bc[idx_tipo3_violado]

            # Atualiza para tipo 0 as barras que estão na lista de barras controladas
            tipo[np.isin(barra, barras_controladas)] = 0

            # Corrigir qesp considerando ambas as condições
            qesp = np.where(violou_sup, qmax, qesp)
            qesp = np.where(violou_inf, qmin, qesp)

            delta_residuos, pcalc, qcalc = calcula_residuos(v, teta, ybarra, tipo, pesp, qesp, dbar, CREM)
        

        erro_max = np.max(np.abs(delta_residuos))
        print(f"Iteração {it:>3}: Erro máximo = {erro_max:.4e}")

        # Verifica se o erro máximo é menor que a tolerância
        if erro_max < tol:
            print(f"Convergência atingida em {it} iterações.")
            return v, teta, pcalc, qcalc, True


        jacobiano = calcula_jacobiano(v, teta, pcalc, qcalc, ybarra, tipo, barra, bc, CREM)
        try:
            delta_var_estado = np.linalg.solve(jacobiano, delta_residuos)
        except np.linalg.LinAlgError:
            print("Erro: Matriz Jacobiana singular.")
            break

        teta += delta_var_estado[:nbar]
        v += delta_var_estado[nbar:nbar*2]

        if CREM and np.any(tipo == 3):
            qesp[idx_tipo3] += delta_var_estado[nbar*2:nbar*2+nbar_tipo3]


    print("Fluxo de potência não convergiu!")
    return v, teta, pcalc, qcalc, False


# === Função para impressão de resultados === #
def imprime_resultados_barras(dbar, v, teta, pcalc, qcalc, pbase=100, casa_decimal=6):
    tipo_map = {0: "PQ", 1: "PV", 2: "Slack"}
    dados = {
        "Barra": dbar.index + 1,
        "Tipo": [tipo_map[t] for t in dbar.tipo],
        "V (pu)": np.round(v, casa_decimal),
        "Teta (graus)": np.round(np.rad2deg(teta), casa_decimal),
        "P Inj (MW)": np.round(pcalc * pbase, casa_decimal),
        "Pg (MW)": np.round((pcalc + dbar.pl.values) * pbase, casa_decimal),
        "Pl (MW)": np.round(dbar.pl.values * pbase, casa_decimal),
        "Q Inj (MVAr)": np.round(qcalc * pbase, casa_decimal),
        "Qg (MVAr)": np.round((qcalc + dbar.ql.values) * pbase, casa_decimal),
        "Ql (MVAr)": np.round(dbar.ql.values * pbase, casa_decimal)
    }

    df = pd.DataFrame(dados)
    print("\n=== Resultados do Fluxo de Potência ===")
    print(df.to_string(index=False))

def imprime_resultados_circuitos(dlin, v, teta, pbase=100, casa_decimal=6):
    print("\n=== Fluxo de Potência nas Linhas (sentido DE → PARA) ===")
    linhas = []

    for idx, row in dlin.iterrows():
        de = int(row.de) - 1
        para = int(row.para) - 1
        r, x, bsh, tap, defas = row.r, row.x, row.bsh, row.tap, row.defasagem

        # Admitância da linha
        z = complex(r, x)
        y = 1 / z
        y_shunt = 1j * bsh

        # Tensões complexas
        vk = v[de] * np.exp(1j * teta[de])
        vm = v[para] * np.exp(1j * teta[para])

        # Corrente de -> para
        ikm = ((vk - vm * np.exp(-1j * defas)) / (tap * z)) + (vk * y_shunt)

        # Potência complexa
        skm = vk * np.conj(ikm)
        pkm = np.real(skm) * pbase
        qkm = np.imag(skm) * pbase

        linhas.append({
            "DE": row.de,
            "PARA": row.para,
            "Pkm (MW)": round(pkm, casa_decimal),
            "Qkm (MVAr)": round(qkm, casa_decimal),
            "R": round(r, casa_decimal),
            "X": round(x, casa_decimal),
            "Bsh": round(bsh, casa_decimal),
            "TAP": round(tap, casa_decimal),
            "Defasagem (rad)": round(defas, casa_decimal),
            "Pmax (MW)": round(row.pkmax, casa_decimal) if "pkmax" in row else None
        })

    df_fluxo = pd.DataFrame(linhas)
    print(df_fluxo.to_string(index=False))

def imprime_balanco_potencia(dbar, pcalc, pbase=100, casa_decimal=6):
    tipo = dbar.tipo.values

    p_slack = pcalc[tipo==2].item() * pbase

    pliq = (dbar.pg.values - dbar.pl.values) * pbase
    pliq[tipo == 2] = 0

    carga = np.sum(np.abs(pliq[pliq < 0]))
    geracao = np.sum(np.abs(pliq[pliq > 0]))
  
    # Perdas = geração - carga (baseado nas potências líquidas injetadas)
    perdas_ativas = p_slack + geracao - carga

    print("\n=== Balanço de Potência ===")
    print(f"→ Carga Total      : {carga:,.{casa_decimal}f} MW")
    print(f"→ Geração Total    : {geracao + p_slack:,.{casa_decimal}f} MW")
    print(f"→ Perdas no Sistema: {perdas_ativas:,.{casa_decimal}f} MW ")


# === Configuração Inicial === #
arquivo_pwf = 'arquivos_do_trabalho/IEEE14_Caso2.pwf'
# arquivo_pwf = 'exemplo_QLIM.pwf'
pbase = 100.
tol = 1e-5	
iter_max = 25

# === Dados de Entrada === #
dbar, dlin = read_pwf(arquivo_pwf)
# imprime_info_dbar(dbar)
# imprime_info_dlin(dlin)
dbar, dlin = inicializa_dbar_dlin(dbar, dlin, pbase)

# === Execução Principal === #
v, teta, pcalc, qcalc, convergiu = calcula_fluxo_de_potencia_newt(
    dbar, dlin, tol, iter_max, flat_start=True, QLIM=True, CREM=False
)

# === Saídas para simples de verificação === #
imprime_resultados_barras(dbar, v, teta, pcalc, qcalc, pbase, 3)
imprime_resultados_circuitos(dlin, v, teta, pbase, 3)
imprime_balanco_potencia(dbar, pcalc, pbase, 3)