# Power Flow Analysis

A Python-based power flow analysis tool that can parse and analyze IEEE standard power flow files.

## Features

- Parse IEEE standard PWF (Power Flow) files
- Extract bus data (DBAR) and branch data (DLIN)
- Display power system information in a formatted table
- Support for voltage, power generation, and load data

## File Format

The project works with IEEE format PWF files that contain:

- DBAR: Bus data section with bus number, type, voltage, generation, and load information
- DLIN: Branch data section with line/transformer parameters
- DGLT: System parameters

### DBAR Format

Each bus record contains:

- Bus number (columns 1-5)
- Type (column 8): 0=PQ, 1=PV, 2=slack
- Voltage in p.u. (columns 25-29): Stored as 5 digits (10600. means 1.06 p.u.)
- Active generation PG (columns 33-37)
- Reactive generation QG (columns 38-42)
- Minimum reactive generation QN (columns 43-47)
- Maximum reactive generation QM (columns 48-52)
- Active load PL (columns 59-63)
- Reactive load QL (columns 64-68)
- Shunt component SH (columns 69-73)

### DLIN Format

Each branch record contains:

- From bus (columns 1-5)
- To bus (columns 11-15)
- Circuit ID (columns 16-17)
- Resistance in % (columns 21-26)
- Reactance in % (columns 27-32)
- Susceptance in Mvar (columns 33-38)
- Tap ratio (columns 39-43)
- Phase shift (columns 54-58)

## Usage

```python
# Parse a PWF file
python leitura14.py

# Output example:
# Bus Information:
# ====================================================================================================
# Bus  | Type  |  Voltage   |    PG    |    QG    |    QN    |    QM    |    PL    |    QL    |    SH
# ----------------------------------------------------------------------------------------------------
#  1   |   2   |    1.06    |  232.4   |  -16.9   |   None   |   None   |   None   |   None   |   None
#  2   |   1   |   1.045    |   40.0   |   42.4   |  -40.0   |   50.0   |   21.7   |   12.7   |   None
# ...
```

## Files

- `leitura14.py`: Python parser for PWF files
- `ieee14.pwf`: IEEE 14 bus test case
- `Leitura14.m`: MATLAB version of the parser (reference implementation)

## Requirements

- Python 3.6+
- NumPy

## License

This project is available as open source under the terms of the MIT License.
