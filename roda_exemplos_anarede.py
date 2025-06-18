from pathlib import Path
import numpy as np
import pandas as pd
import cmath
from leitura_pwf import *
from main_power_flow import *

# ================================================================================= #
# CASO 1:
# ================================================================================= #
# === Configuração Inicial === #
arquivo_pwf = 'arquivos_do_trabalho/IEEE14_Caso1.pwf'
# arquivo_pwf = 'arquivos_do_trabalho/teste.pwf'
pbase = 100.
tol = 1e-8	
iter_max = 25

# === Dados de Entrada === #
dbar, dlin = read_pwf(arquivo_pwf)
dbar, dlin = inicializa_dbar_dlin(dbar, dlin, pbase)

# === Execução Principal === #
v, teta, pcalc, qcalc, convergiu = calcula_fluxo_de_potencia_newt(
    dbar, dlin, tol, iter_max, flat_start=True, QLIM=False, CREM=True, CTAP=False
)

# === Entradas para simples de verificação === #
imprime_info_dbar(dbar)
imprime_info_dlin(dlin)
# === Saídas para simples de verificação === #
imprime_resultados_barras(dbar, v, teta, pcalc, qcalc, pbase, 3)
imprime_resultados_circuitos(dlin, v, teta, pbase, 3)
imprime_balanco_potencia(dbar, pcalc, pbase, 3)


# ================================================================================= #
# CASO 2:
# ================================================================================= #



# ================================================================================= #
# CASO 3:
# ================================================================================= #



# ================================================================================= #
# CASO 4:
# ================================================================================= #



# ================================================================================= #
# CASO 5:
# ================================================================================= #