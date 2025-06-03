from pathlib import Path
import numpy as np
import pandas as pd
import cmath
from leitura_pwf import *


# === Configuração Inicial === #
arquivo_pwf = 'ieee14.pwf'
Pbase = 100.
TOL = 0.003
iter_max = 25


# === Dados de Entrada === #
DBAR, DLIN = read_pwf(arquivo_pwf)

dbar_dataframe = converter_dbar(DBAR)
dlin_dataframe = converter_dlin(DLIN)

#print_dbar_info(DBAR)
#print_dlin_info(DLIN)


# === Funções para transformar os dados de Entrada em pu === #
def inicializa_barras(dbar, pbase):
    dbar = dbar.copy()
    dbar[["pg", "qg", "qn", "qm", "pd", "qd", "shunt"]] /= pbase
    dbar["teta"] = np.deg2rad(dbar["teta"])
    return dbar

def inicializa_linhas(dlin, pbase):
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


# === Funções para cálculo de Fluxo de Potência === #
def calcular_residuos(v, teta, ybarra, tipo, pesp, qesp):
    x, y = np.cos(teta) * v, np.sin(teta) * v
    vret = x + 1j * y
    i_inj = ybarra @ vret
    s = vret * np.conj(i_inj)
    pcalc = np.real(s)
    qcalc = np.imag(s)

    delta_p = pesp - pcalc
    delta_q = qesp - qcalc

    delta_p[tipo == 2] = 0  # barra swing
    delta_q[tipo == 2] = 0
    delta_q[tipo == 1] = 0  # barra PV

    delta_y = np.concatenate([delta_p, delta_q])
    return delta_y, pcalc, qcalc

def montar_jacobiano(v, teta, pcalc, qcalc, ybarra, tipo):
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

def fluxo_de_potencia(dbar, ybarra, tol=1e-4, iter_max=25):
    nbar = len(dbar)
    tipo = dbar.tipo.values
    v = dbar.v.values.copy()
    teta = dbar.teta.values.copy()
    pg, pd = dbar.pg.values, dbar.pd.values
    qg, qd = dbar.qg.values, dbar.qd.values
    pesp, qesp = pg - pd, qg - qd

    for it in range(iter_max):
        delta_y, pcalc, qcalc = calcular_residuos(v, teta, ybarra, tipo, pesp, qesp)
        erro_max = np.max(np.abs(delta_y))
        print(f"Iteração {it:>3}: Erro máximo = {erro_max:.4e}")

        if erro_max < tol:
            print(f"Convergência atingida em {it} iterações. Erro máximo: {erro_max:.4e}")
            return v, teta, pcalc, qcalc, True

        Jac = montar_jacobiano(v, teta, pcalc, qcalc, ybarra, tipo)
        try:
            delta_sol = np.linalg.solve(Jac, delta_y)
        except np.linalg.LinAlgError:
            print("Erro: Matriz Jacobiana singular.")
            break

        teta += delta_sol[:nbar]
        v += delta_sol[nbar:]

    print("Fluxo de potência não convergiu!")
    return v, teta, pcalc, qcalc, False



# === Execução Principal === #
dbar = inicializa_barras(dbar_dataframe, Pbase)
dlin = inicializa_linhas(dlin_dataframe, Pbase)

NBAR = len(dbar)
NLIN = len(dlin)
TIPO = dbar.tipo.values

Ybarra = calcula_ybarra(NBAR, dlin, dbar.shunt.values)
G = np.real(Ybarra)
B = np.imag(Ybarra)

v, teta, pcalc, qcalc, convergiu = fluxo_de_potencia(dbar, Ybarra, TOL, iter_max)

# Saída simples de verificação
print("Tensões finais (pu):", ["{:.6f}".format(val) for val in v])
print("Angulos finais (graus):", ["{:.6f}".format(val) for val in np.rad2deg(teta)])