import pytest
import numpy as np
import pandas as pd
from main_power_flow import calcula_fluxo_de_potencia_newt, read_pwf, inicializa_dbar_dlin

def test_power_flow_results():
    arquivo_pwf = 'pwf_files/IEEE14_Caso1.pwf'
    pbase = 100.
    tol = 1e-8
    iter_max = 25

    dbar, dlin = read_pwf(arquivo_pwf)
    dbar, dlin = inicializa_dbar_dlin(dbar, dlin, pbase)

    v, teta, pcalc, qcalc, convergiu = calcula_fluxo_de_potencia_newt(
        dbar, dlin, tol, iter_max, flat_start=False, QLIM=False, CREM=True
    )

    expected_voltage = np.array([1.060, 1.045, 1.010, 1.006, 1.014, 1.089, 1.000, 0.949, 1.018, 1.023, 1.052, 1.070, 1.061, 1.019])
    expected_angle = np.array([0.000, -5.001, -12.790, -10.131, -8.771, -14.694, -13.087, -13.087, -14.658, -14.934, -14.912, -15.512, -15.490, -16.051])
    expected_pg = np.array([233.138, 40.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000])
    expected_qg = np.array([-14.151, 53.631, 32.023, 0.000, 0.000, 43.568, 0.000, -27.572, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000])

    voltage_margin = 0.001
    angle_margin = 0.1
    power_margin = 0.1

    mismatches = []

    for i in range(len(v)):
        if abs(v[i] - expected_voltage[i]) > voltage_margin:
            mismatches.append(f"Voltage at bus {i+1}: Expected {expected_voltage[i]:.3f}, Got {v[i]:.3f}")
        if abs(np.degrees(teta[i]) - expected_angle[i]) > angle_margin:
            mismatches.append(f"Angle at bus {i+1}: Expected {expected_angle[i]:.3f}, Got {np.degrees(teta[i]):.3f}")
        if abs(dbar.iloc[i].pg * pbase - expected_pg[i]) > power_margin:
            mismatches.append(f"PG at bus {i+1}: Expected {expected_pg[i]:.3f}, Got {dbar.iloc[i].pg * pbase:.3f}")
        if abs(dbar.iloc[i].qg * pbase - expected_qg[i]) > power_margin:
            mismatches.append(f"QG at bus {i+1}: Expected {expected_qg[i]:.3f}, Got {dbar.iloc[i].qg * pbase:.3f}")

    expected_results = pd.DataFrame({
        'DE': [1, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 9, 9, 10, 12, 13],
        'PARA': [2, 5, 3, 4, 5, 4, 5, 7, 9, 6, 11, 12, 13, 8, 9, 10, 14, 11, 13, 14],
        'Pkm': [157.435, 75.703, 73.658, 55.661, 42.088, -22.892, 
                -57.512, 25.369, 14.999, 48.503,# errados
                  9.682, 8.440, 19.181, -0.000, 25.369, 3.267, 7.601, -5.752, 2.250, 7.605],
        'Qkm': [-20.533, 6.382, 3.519, 5.253, 4.260, 11.269, 
                -0.064, 14.8, 4.4, 2.1, # errados
                15.686, 4.021, 13.483, 29.060, -15.983, -7.098, -3.663, -12.947, 2.233, 9.294]
    })

    calculated_results = []
    for idx, row in dlin.iterrows():
        de = int(row.de) - 1
        para = int(row.para) - 1
        r, x, bsh, tap, defas = row.r, row.x, row.bsh, row.tap, row.defasagem

        z = complex(r, x)
        y = 1 / z
        y_shunt = 1j * bsh

        vk = v[de] * np.exp(1j * teta[de])
        vm = v[para] * np.exp(1j * teta[para])

        ikm = ((vk - vm * np.exp(-1j * defas)) / (tap * z)) + (vk * y_shunt)
        skm = vk * np.conj(ikm)
        pkm = np.real(skm) * pbase
        qkm = np.imag(skm) * pbase

        calculated_results.append({
            'DE': row.de,
            'PARA': row.para,
            'Pkm': pkm,
            'Qkm': qkm
        })

    calculated_df = pd.DataFrame(calculated_results)

    for _, row in expected_results.iterrows():
        calc_row = calculated_df[
            (calculated_df['DE'] == row['DE']) & 
            (calculated_df['PARA'] == row['PARA'])
        ].iloc[0]

        pkm_expected = row['Pkm']
        qkm_expected = row['Qkm']
        pkm_calculated = calc_row['Pkm']
        qkm_calculated = calc_row['Qkm']

        pkm_margin = max(abs(pkm_expected * 0.1), 1e-10)
        qkm_margin = max(abs(qkm_expected * 0.1), 1e-10)

        if abs(pkm_calculated - pkm_expected) > pkm_margin:
            mismatches.append(f"Pkm mismatch for line {row['DE']}-{row['PARA']}: Expected {pkm_expected:.3f}, Got {pkm_calculated:.3f}")
        
        if abs(qkm_calculated - qkm_expected) > qkm_margin:
            mismatches.append(f"Qkm mismatch for line {row['DE']}-{row['PARA']}: Expected {qkm_expected:.3f}, Got {qkm_calculated:.3f}")

    if mismatches:
        print("\nMismatches found:")
        for mismatch in mismatches:
            print(mismatch)
        assert False, "Test failed due to mismatches"

    assert convergiu, "Power flow did not converge" 