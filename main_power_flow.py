from pathlib import Path
import numpy as np
import pandas as pd
import cmath
from leitura_pwf import *



# === Funções para transformar os dados de Entrada em pu === #
def inicializa_barras(dbar, pbase=100):
    dbar = dbar.copy()
    dbar[["pg", "qg", "qn", "qm", "pd", "qd", "shunt"]] /= pbase
    dbar["teta"] = np.deg2rad(dbar["teta"])
    return dbar

def inicializa_linhas(dlin, pbase=100):
    dlin = dlin.copy()
    dlin["r"] /= 100
    dlin["x"] /= 100
    dlin["bsh"] = (dlin["bsh"] / 2) / pbase
    return dlin

# === Função para cálculo da Ybarra === #
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

# === Funções auxiliares para cálculo de Fluxo de Potência === #
def calcula_residuos(v, teta, ybarra, tipo, pesp, qesp):

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
    return delta_y, pcalc, qcalc

def calcula_jacobiano(v, teta, pcalc, qcalc, ybarra, tipo):
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
    return Jac

# === Função para cálculo de Fluxo de Potência === #
def calcula_fluxo_de_potencia_newt(dbar, ybarra, tol=1e-4, iter_max=25):
    nbar = len(dbar)
    tipo = dbar.tipo.values
    v = dbar.v.values.copy()
    teta = dbar.teta.values.copy()
    pg, pd = dbar.pg.values, dbar.pd.values
    qg, qd = dbar.qg.values, dbar.qd.values
    pesp, qesp = pg - pd, qg - qd

    for it in range(iter_max):
        delta_y, pcalc, qcalc = calcula_residuos(v, teta, ybarra, tipo, pesp, qesp)
        erro_max = np.max(np.abs(delta_y))
        print(f"Iteração {it:>3}: Erro máximo = {erro_max:.4e}")

        if erro_max < tol:
            print(f"Convergência atingida em {it} iterações. Erro máximo: {erro_max:.4e}")
            return v, teta, pcalc, qcalc, True

        Jac = calcula_jacobiano(v, teta, pcalc, qcalc, ybarra, tipo)
        try:
            delta_sol = np.linalg.solve(Jac, delta_y)
        except np.linalg.LinAlgError:
            print("Erro: Matriz Jacobiana singular.")
            break

        teta += delta_sol[:nbar]
        v += delta_sol[nbar:]

    print("Fluxo de potência não convergiu!")
    return v, teta, pcalc, qcalc, False

# === Função para impressão de resultados === #
def imprime_balanco_potencia(dbar, pcalc, qcalc, pbase=100):
    pd_total = dbar.pd.sum() * pbase
    qd_total = dbar.qd.sum() * pbase
    pg_total = dbar.pg.sum() * pbase
    qg_total = dbar.qg.sum() * pbase

    # Perdas = geração - carga (baseado nas potências líquidas injetadas)
    perdas_ativas = pg_total - pd_total
    perdas_reativas = qg_total - qd_total

    print("\n=== Balanço de Potência ===")
    print(f"→ Carga Total      : {pd_total:,.6f} MW  | {qd_total:,.6f} MVAr")
    print(f"→ Geração Total    : {pg_total:,.6f} MW  | {qg_total:,.6f} MVAr")
    print(f"→ Perdas no Sistema: {perdas_ativas:,.6f} MW ")# | {perdas_reativas:,.6f} MVAr")

def imprime_resultados_barras(dbar, v, teta, pcalc, qcalc, pbase=100):
    tipo_map = {0: "PQ", 1: "PV", 2: "Slack"}
    dados = {
        "Barra": dbar.index + 1,
        "Tipo": [tipo_map[t] for t in dbar.tipo],
        "V (pu)": np.round(v, 6),
        "Teta (graus)": np.round(np.rad2deg(teta), 6),
        "P Inj (MW)": np.round(pcalc * pbase, 6),
        "Q Inj (MVAr)": np.round(qcalc * pbase, 6)
    }

    df = pd.DataFrame(dados)
    print("\n=== Resultados do Fluxo de Potência ===")
    print(df.to_string(index=False))

def imprime_resultados_circuitos(dlin, v, teta, ybarra, pbase=100):
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
            "Pkm (MW)": round(pkm, 6),
            "Qkm (MVAr)": round(qkm, 6),
            "R": round(r, 6),
            "X": round(x, 6),
            "Bsh": round(bsh, 6),
            "TAP": round(tap, 6),
            "Defasagem (rad)": round(defas, 6),
            "Pmax (MW)": round(row.pkmax, 6) if "pkmax" in row else None
        })

    df_fluxo = pd.DataFrame(linhas)
    print(df_fluxo.to_string(index=False))



# === Dados de Entrada === #
dbar, dlin = read_pwf(arquivo_pwf)

dbar_dataframe = converter_dbar(dbar)
dlin_dataframe = converter_dlin(dlin)

#print_dbar_info(DBAR)
#print_dlin_info(DLIN)

# === Configuração Inicial === #
arquivo_pwf = 'ieee14.pwf'
pbase = 100.
tol = 0.0000003
iter_max = 25

# === Execução Principal === #
dbar = inicializa_barras(dbar_dataframe, pbase)
dlin = inicializa_linhas(dlin_dataframe, pbase)

ybarra = calcula_ybarra(len(dbar), dlin, dbar.shunt.values)

v, teta, pcalc, qcalc, convergiu = calcula_fluxo_de_potencia_newt(dbar, ybarra, tol, iter_max)

# === Saídas para simples de verificação === #
imprime_balanco_potencia(dbar, pcalc, qcalc, pbase)
imprime_resultados_barras(dbar, v, teta, pcalc, qcalc, pbase)
imprime_resultados_circuitos(dlin, v, teta, ybarra, pbase)