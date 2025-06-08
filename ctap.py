from pathlib import Path
import numpy as np
import pandas as pd
import cmath
from leitura_pwf import *

# Passos básicos do controle de tap:
# Medir a tensão atual na barra controlada.

# Comparar a tensão atual com a tensão de referência.

# Se a tensão estiver fora da faixa aceitável (acima ou abaixo da referência com uma tolerância), o tap é ajustado:

# Se tensão > referência + tolerância, reduz o tap (diminui a tensão).

# Se tensão < referência - tolerância, aumenta o tap (aumenta a tensão).

# Atualizar o valor do tap e recalcular a matriz admitância do sistema 
# 𝑌
# Y com o novo tap.

# Refazer o fluxo de potência para encontrar as novas tensões e potências.

# Repetir o processo até que as tensões estejam dentro dos limites ou um número máximo de iterações seja alcançado.

def barras_controladas_tap(df):
    mascara = df['bc'] != 0
    
    # Obtém os valores e as posições como vetores NumPy
    de = np.array(df.loc[mascara, 'de'].tolist())
    para = np.array(df.loc[mascara, 'para'].tolist())
    bc = np.array(df.loc[mascara, 'bc'].tolist())
    tmin = np.array(df.loc[mascara, 'tmin'].tolist())
    tmax = np.array(df.loc[mascara, 'tmax'].tolist())
    tap = np.array(df.loc[mascara, 'tap'].tolist())
    posicoes = np.array(df.index[mascara].tolist())

   
    # Empilha em uma matriz
    resultado = np.vstack((de, para, bc, tmin, tmax, tap, posicoes))
    #print(f"bcs {resultado}")
    return resultado

def verifica_tap(barras_ctrl, tensao_calc, tensao_ref):
    passo_tap = 0.01
    alteracoes_tap = []

    for idx, (de, para, barra_ctrl, tapmin, tapmax, tap, ind) in enumerate(zip(barras_ctrl[0,:], barras_ctrl[1,:], barras_ctrl[2,:], barras_ctrl[3,:], barras_ctrl[4,:], barras_ctrl[5,:],  barras_ctrl[6,:])):
        indice_barra = int(barra_ctrl) - 1

        if tensao_calc[indice_barra] > tensao_ref[indice_barra]  and tap > tapmin:
            tap_new = tap - passo_tap
            alteracoes_tap.append([int(barra_ctrl), float(tap_new), int(de), int(para), int(ind)])
            print(f"Reduzindo tap que controla a barra {barra_ctrl} de {tap} para {tap_new} trafo de {de} para {para}")

        elif tensao_calc[indice_barra] < tensao_ref[indice_barra] and tap < tapmax:
            tap_new = tap + passo_tap
            alteracoes_tap.append([int(barra_ctrl), float(tap_new), int(de), int(para), int(ind)])
            print(f"Aumentando tap que controla a barra {barra_ctrl} de {tap} para {tap_new}")

    return alteracoes_tap

def atualizar_Ybarra(ybarra, dlin, taps_alterados):

    taps_alterados = np.array(taps_alterados)
    matriz_reduzida = taps_alterados[:, [4, 1]] #tirando infos desnecessárias
    taps_alterados = matriz_reduzida[:, [0,1]] #coluna 0: índice do trafo no DLIN e coluna 1: novo valor do tap

    Y = ybarra.copy()

    for linha in taps_alterados:
        idx = int(linha[0])
        novo_tap = linha[1]

        # Obtem a linha no dlin correspondente
        linha_dlin = dlin.iloc[idx]

        # Índices das barras i e j (convertendo para zero-based)
        i = int(linha_dlin.de) - 1
        j = int(linha_dlin.para) - 1
        r = linha_dlin.r
        x = linha_dlin.x
        bsh = linha_dlin.bsh
        defas = linha_dlin.defasagem

        # Calcula a admitância da linha
        y = 1 / complex(r, x)

        # Atualiza o tap com o novo valor
        tap = novo_tap

        # Zera os valores antigos da linha i-j na Ybarra
        # (Importante para não somar repetidamente)
        Y[i, i] = 0
        Y[j, j] = 0
        Y[i, j] = 0
        Y[j, i] = 0

        # Recalcula os elementos da Ybarra com o novo tap
        Y[i, i] += (1 / tap**2) * y + 1j * bsh
        Y[j, j] += y + 1j * bsh
        Y[i, j] -= y / tap * np.exp(-1j * defas)
        Y[j, i] -= y / tap * np.exp(1j * defas)

        print(f"Atualizado ramo {idx}: barras ({i+1}, {j+1}), novo tap = {tap}")

    return Y

def calcula_jacobiano_com_tap(v, teta, pcalc, qcalc, ybarra, tipo, barras_ctrl, taps, dlin):
    nbar = len(v)
    g = np.real(ybarra)
    b = np.imag(ybarra)

    n_ctrl = len(barras_ctrl)

    # Matriz base H N M L
    H = np.zeros((nbar, nbar))
    N = np.zeros((nbar, nbar))
    M = np.zeros((nbar, nbar))
    L = np.zeros((nbar, nbar))

    # Matriz extra com colunas para os taps
    H_tap = np.zeros((nbar, n_ctrl))
    M_tap = np.zeros((nbar, n_ctrl))

    for k in range(nbar):
        for m in range(nbar):
            if k == m:
                H[k, k] = -(qcalc[k] + v[k]**2 * b[k, k])
                N[k, k] = (pcalc[k] + v[k]**2 * g[k, k]) / v[k]
                M[k, k] = pcalc[k] - v[k]**2 * g[k, k]
                L[k, k] = (qcalc[k] - v[k]**2 * b[k, k]) / v[k]

                if tipo[k] == 2:   # Barra Slack
                    H[k, k] = 1e10
                    L[k, k] = 1e10
                elif tipo[k] == 1:  # Barra PV
                    L[k, k] = 1e10
            else:
                delta = teta[k] - teta[m]
                H[k, m] = v[k]*v[m]*(g[k, m]*np.sin(delta) - b[k, m]*np.cos(delta))
                N[k, m] = v[k]*(g[k, m]*np.cos(delta) + b[k, m]*np.sin(delta))
                M[k, m] = -v[k]*v[m]*(g[k, m]*np.cos(delta) + b[k, m]*np.sin(delta))
                L[k, m] = v[k]*(g[k, m]*np.sin(delta) - b[k, m]*np.cos(delta))

    # Derivadas parciais P,Q em relação aos taps controlados
    for idx, barra_ctrl in enumerate(barras_ctrl):
        # Localização da barra controlada
        i = int(barra_ctrl) - 1
        tap = taps[idx]

        # Procurar no dlin a linha correspondente ao trafo que controla essa barra
        for linha_idx in range(len(dlin)):
            linha = dlin.iloc[linha_idx]
            if int(linha.de) -1 == i or int(linha.para) -1 == i:
                j = int(linha.para) -1 if int(linha.de)-1 == i else int(linha.de)-1

                r = linha.r
                x = linha.x
                bsh = linha.bsh
                defas = linha.defasagem

                y = 1 / complex(r, x)
                g_ij = np.real(y)
                b_ij = np.imag(y)

                # Derivadas (simplificadas)
                dP_dtap = -2 / tap**3 * v[i]**2 * g[i,i] + v[i]*v[j]/tap**2 * (g_ij*np.cos(teta[i]-teta[j]) + b_ij*np.sin(teta[i]-teta[j]))
                dQ_dtap = -2 / tap**3 * v[i]**2 * b[i,i] + v[i]*v[j]/tap**2 * (g_ij*np.sin(teta[i]-teta[j]) - b_ij*np.cos(teta[i]-teta[j]))

                # Preenche as colunas extra
                H_tap[i, idx] = dP_dtap
                M_tap[i, idx] = dQ_dtap

                break   # achou o trafo, não precisa buscar mais

    # Expande Jacobiana
    H_exp = np.hstack((H, H_tap))
    N_exp = np.hstack((N, np.zeros((nbar, n_ctrl))))
    M_exp = np.hstack((M, M_tap))
    L_exp = np.hstack((L, np.zeros((nbar, n_ctrl))))

    # Monta Jacobiana final
    Jac = np.vstack((
        np.hstack((H_exp, N_exp)),
        np.hstack((M_exp, L_exp))
    ))


    return Jac
