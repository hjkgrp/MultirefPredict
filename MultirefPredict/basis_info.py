import basis_set_exchange
import re
import qcelemental as qcel

def element_to_num_AO(symbol, basis_set):
    parsed_basis_set = basis_set
    if basis_set.lower() == "lacvps_ecp" or basis_set.lower() == "lacvp*_ecp":
        atomic_num = qcel.periodictable.to_atomic_number(symbol)
        if atomic_num <= 18:
            parsed_basis_set = "6-31g*"
        else:
            parsed_basis_set = "lanl2dz"
    basis_str = basis_set_exchange.get_basis(parsed_basis_set 
                ,elements = [symbol],fmt = 'turbomole', header = False) 
    basis_lines = basis_str.split('\n')
    AO = 0
    for line in basis_lines:
        mobj = re.search(r'\s+\d+\s+[a-z]+', line )
        if mobj:
            if line.split()[1] == "s":
                AO += 1
            if line.split()[1] == "p":
                AO += 3
            if line.split()[1] == "d":
                AO += 6
    return AO

def molecule_to_num_AO(molecule, basis_set):
    elements = {}
    # Build a dictionary of "element":count
    for symbol in molecule.symbols:
         if symbol not in elements:
             elements[symbol] = 1
         else:
             elements[symbol] += 1
    # Sum up the AO for all atoms
    AO = 0
    for symbol in elements:
        AO += element_to_num_AO(symbol, basis_set)*elements[symbol]

    core_AO = 0
    for symbol in elements:
        core_AO += element_to_core_AO(symbol, basis_set)*elements[symbol]
  
    return AO, core_AO

def element_to_core_AO(symbol, basis_set):
    core_AO = 0
    using_ecp = False
    # May add more cases about ecp in the future
    if basis_set.lower() in [ "lacvp*_ecp", "lacvps_ecp", "lanl2dz" ]:
        using_ecp = True
    
    atomic_num = qcel.periodictable.to_atomic_number(symbol)
    if atomic_num <= 2:
        core_AO = 0
    elif atomic_num in range(3,11):
        core_AO = 1
    elif atomic_num in range(11,19):
        core_AO = 5
    elif atomic_num > 18:
        if not using_ecp:
            if atomic_num in range(19,37):
                core_AO = 5
            elif atomic_num in range(37,55):
                core_AO = 9  #1s2s2p3s3p
            elif atomic_num in range(55,87):
                core_AO = 18
            else:
                core_AO = 27
        else:
            core_AO = 0

    return core_AO
