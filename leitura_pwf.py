import numpy as np
import pandas as pd
import math


def read_pwf(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    DBAR = []
    DLIN = []
    
    line_index = 0
    
    while line_index < len(lines):
        line = lines[line_index].strip()
        
        if line == 'DBAR':
            line_index += 1 
            
            while line_index < len(lines):
                line = lines[line_index]
                
                if line.strip().startswith('99999'):
                    break
                if not line.strip() or line.strip().startswith('('):
                    line_index += 1
                    continue
                
                # Initialize with None for all fields (except bus number and type)
                dbar_row = [0, 0, None, None, None, None, None, None, None, None, None, None]
                
                # Bus number (columns 1-5)
                barra_str = line[0:5].strip()
                if barra_str:
                    dbar_row[0] = int(barra_str)
                
                # Type (column 8)
                tipo_str = line[7:8].strip()
                if tipo_str:
                    dbar_row[1] = int(tipo_str)
                
                # Voltage (columns 25-28)
                voltage_str = line[24:28].strip()
                if voltage_str:
                    # For voltage only, convert 10600. to 1.06
                    dbar_row[2] = float(voltage_str) / 1000.0

                # Angle
                teta_str = line[28:32].strip()
                if teta_str:
                    dbar_row[3] = float(teta_str)   
                
                # Active generation PG (columns 33-37)
                pg_str = line[32:37].strip()
                if pg_str:
                    dbar_row[4] = float(pg_str)
                
                # Reactive generation QG (columns 38-42)
                qg_str = line[37:42].strip()
                if qg_str:
                    dbar_row[5] = float(qg_str)
                
                # Minimum reactive generation QN (columns 43-47)
                qn_str = line[42:47].strip()
                if qn_str:
                    dbar_row[6] = float(qn_str)
                elif dbar_row[1] == 2:
                    dbar_row[6] = -9999.0    
                
                # Maximum reactive generation QM (columns 48-52)
                qm_str = line[47:52].strip()
                if qm_str:
                    dbar_row[7] = float(qm_str)
                elif dbar_row[1] == 2:
                    dbar_row[7] = 99999.0

                # Control Bus
                bc_str = line[52:58].strip()
                if bc_str:
                    dbar_row[8] = int(bc_str)
                
                # Active load PL (columns 59-63)
                pl_str = line[58:63].strip()
                if pl_str:
                    dbar_row[9] = float(pl_str)
                
                # Reactive load QL (columns 64-68)
                ql_str = line[63:68].strip()
                if ql_str:
                    dbar_row[10] = float(ql_str)
                
                # Shunt component SH (columns 69-73)
                shunt_str = line[68:73].strip()
                if shunt_str:
                    dbar_row[11] = float(shunt_str)
                
                DBAR.append(dbar_row)
                
                line_index += 1
                if line_index < len(lines):
                    line = lines[line_index]
                else:
                    break
        
        elif line == 'DLIN':
            line_index += 1  
            
            while line_index < len(lines):
                line = lines[line_index]    

                if line.strip().startswith('99999'):
                    break
                if not line.strip() or line.strip().startswith('('):
                    line_index += 1
                    continue
                
                # Initialize with None for all fields
                dlin_row = [0, 0, None, None, None, None, None]
                
                # From bus (columns 1-5)
                from_str = line[0:5].strip()
                if from_str:
                    dlin_row[0] = int(from_str)
                
                # To bus (columns 11-15)
                to_str = line[10:15].strip()
                if to_str:
                    dlin_row[1] = int(to_str)
                
                # Circuit (columns 16-17)
                circuit_str = line[15:17].strip()
                if circuit_str:
                    circuit = int(circuit_str)
                else:
                    circuit = 1
                
                # Resistance (columns 21-26)
                r_str = line[20:26].strip()
                if r_str:
                    dlin_row[2] = float(r_str)
                
                # Reactance (columns 27-32)
                x_str = line[26:32].strip()
                if x_str:
                    dlin_row[3] = float(x_str)
                
                # Susceptance (columns 33-38)
                b_str = line[32:38].strip()
                if b_str:
                    dlin_row[4] = float(b_str)
                
                # Tap (columns 39-43)
                tap_str = line[38:43].strip()
                if tap_str:
                    dlin_row[5] = float(tap_str)
                else:
                    dlin_row[5] = 1.0
                
                # Phase shift (columns 54-58)
                if len(line) > 53:
                    phase_str = line[53:58].strip()
                    if phase_str:
                        dlin_row[6] = float(phase_str) / 180 * math.pi
                
                DLIN.append(dlin_row)
                
                line_index += 1
                if line_index < len(lines):
                    line = lines[line_index]
                else:
                    break
        
        line_index += 1

    #return np.array(DBAR, dtype=object), np.array(DLIN, dtype=object)

    DBAR, DLIN = np.array(DBAR, dtype=object), np.array(DLIN, dtype=object)
    
    # Conversão para DataFrame diretamente no final
    DBAR_sem_none = [[0 if val is None else val for val in row] for row in DBAR]
    dbar = pd.DataFrame(DBAR_sem_none, columns=["barra", "tipo", "v", "teta", "pg", "qg", "qn", "qm", "bc", "pl", "ql", "shunt"])
    dbar = dbar.fillna(0).astype({
        "barra": int, "tipo": int, "v": float, "teta": float,
        "pg": float, "qg": float, "qn": float, "qm": float,
        "bc": int, "pl": float, "ql": float, "shunt": float
    })

    DLIN_sem_none = [[0 if val is None else val for val in row] for row in DLIN]
    dlin = pd.DataFrame(DLIN_sem_none, columns=["de", "para", "r", "x", "bsh", "tap", "defasagem"])
    dlin = dlin.fillna(0).astype({
        "de": int, "para": int, "r": float, "x": float,
        "bsh": float, "tap": float, "defasagem": float
    })

    return dbar, dlin

def imprime_info_dbar(dbar):
    print(f"Number of buses: {len(dbar)}")
    print("\nBus Information:")
    print("=" * 120)
    print(f"{'Bus':^5} | {'Type':^5} | {'Voltage':^10} | {'Angle':^8} | {'PG':^8} | {'QG':^8} | {'QN':^8} | {'QM':^8} | {'PL':^8} | {'QL':^8} | {'SH':^8}")
    print("-" * 120)

    for _, row in dbar.iterrows():
        print(f"{row['barra']:^5} | {row['tipo']:^5} | {row['v']:^10.4f} | {row['teta']:^8.2f} | "
              f"{row['pg']:^8.2f} | {row['qg']:^8.2f} | {row['qn']:^8.2f} | {row['qm']:^8.2f} | "
              f"{row['pl']:^8.2f} | {row['ql']:^8.2f} | {row['shunt']:^8.2f}")

def imprime_info_dlin(dlin):
    print(f"Number of branches: {len(dlin)}")
    print("\nBranch Information:")
    print("=" * 80)
    print(f"{'From':^5} | {'To':^5} | {'R (%)':^8} | {'X (%)':^8} | {'B (Mvar)':^10} | {'Tap':^8}")
    print("-" * 80)

    for _, row in dlin.iterrows():
        print(f"{row['de']:^5} | {row['para']:^5} | {row['r']:^8.4f} | {row['x']:^8.4f} | {row['bsh']:^10.4f} | {row['tap']:^8.4f}")


# Teste na leitura
dbar, dlin = read_pwf('p2_ex1.pwf')
imprime_info_dbar(dbar)
imprime_info_dlin(dlin)