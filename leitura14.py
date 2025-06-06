import numpy as np
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
            line_index += 2  # Skip header line
            line = lines[line_index]  # Don't strip, we need exact character positions
            
            while line_index < len(lines) and '99999' not in line:
                if not line.strip():
                    line_index += 1
                    if line_index < len(lines):
                        line = lines[line_index]
                    continue
                
                # Initialize with None for all fields (except bus number and type)
                dbar_row = [0, 0, None, None, None, None, None, None, None, None, None]
                
                # Bus number (columns 1-5)
                bus_str = line[0:5].strip()
                if bus_str:
                    dbar_row[0] = int(bus_str)
                
                # Type (column 8)
                type_str = line[7:8].strip()
                if type_str:
                    dbar_row[1] = int(type_str)
                
                # Voltage (columns 25-29)
                voltage_str = line[24:29].strip()
                if voltage_str:
                    # For voltage only, convert 10600. to 1.06
                    dbar_row[2] = float(voltage_str) / 10000.0
                
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
                
                # Maximum reactive generation QM (columns 48-52)
                qm_str = line[47:52].strip()
                if qm_str:
                    dbar_row[7] = float(qm_str)
                
                # Active load PL (columns 59-63)
                pl_str = line[58:63].strip()
                if pl_str:
                    dbar_row[8] = float(pl_str)
                
                # Reactive load QL (columns 64-68)
                ql_str = line[63:68].strip()
                if ql_str:
                    dbar_row[9] = float(ql_str)
                
                # Shunt component SH (columns 69-73)
                sh_str = line[68:73].strip()
                if sh_str:
                    dbar_row[10] = float(sh_str)
                
                DBAR.append(dbar_row)
                
                line_index += 1
                if line_index < len(lines):
                    line = lines[line_index]
                else:
                    break
        
        elif line == 'DLIN':
            line_index += 2  # Skip header line
            line = lines[line_index]  # Don't strip, we need exact character positions
            
            while line_index < len(lines) and '99999' not in line:
                if not line.strip():
                    line_index += 1
                    if line_index < len(lines):
                        line = lines[line_index]
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
    
    return np.array(DBAR, dtype=object), np.array(DLIN, dtype=object)

def print_dbar_info(dbar):
    print("\nBus Information:")
    print("=" * 100)
    print(f"{'Bus':^5} | {'Type':^5} | {'Voltage':^10} | {'PG':^8} | {'QG':^8} | {'QN':^8} | {'QM':^8} | {'PL':^8} | {'QL':^8} | {'SH':^8}")
    print("-" * 100)
    
    for row in dbar:
        bus_num = int(row[0])
        bus_type = int(row[1])
        voltage = row[2] if row[2] is not None else "None"
        pg = row[4] if row[4] is not None else "None"
        qg = row[5] if row[5] is not None else "None"
        qn = row[6] if row[6] is not None else "None"
        qm = row[7] if row[7] is not None else "None"
        pl = row[8] if row[8] is not None else "None"
        ql = row[9] if row[9] is not None else "None"
        sh = row[10] if row[10] is not None else "None"
        
        print(f"{bus_num:^5} | {bus_type:^5} | {voltage:^10} | {pg:^8} | {qg:^8} | {qn:^8} | {qm:^8} | {pl:^8} | {ql:^8} | {sh:^8}")

def print_dlin_info(dlin):
    print("\nBranch Information:")
    print("=" * 80)
    print(f"{'From':^5} | {'To':^5} | {'R (%)':^8} | {'X (%)':^8} | {'B (Mvar)':^10} | {'Tap':^8}")
    print("-" * 80)
    
    for row in dlin:
        from_bus = int(row[0])
        to_bus = int(row[1])
        r = row[2] if row[2] is not None else 0.0
        x = row[3] if row[3] is not None else 0.0
        b = row[4] if row[4] is not None else 0.0
        tap = row[5] if row[5] is not None else 1.0
        
        print(f"{from_bus:^5} | {to_bus:^5} | {r:^8.4f} | {x:^8.4f} | {b:^10.4f} | {tap:^8.4f}")

if __name__ == "__main__":
    DBAR, DLIN = read_pwf('ieee14.pwf')
    print(f"Number of buses: {len(DBAR)}")
    print(f"Number of lines: {len(DLIN)}")
    
    print_dbar_info(DBAR)
    print_dlin_info(DLIN) 