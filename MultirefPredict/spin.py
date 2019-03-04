import qcelemental

class atomSpinMultDict:
    def __init__(self):
        self.mult_dict = {
            1: 2,
            2: 1,
            3: 2,
            4: 1,
            5: 2,
            6: 3,
            7: 4,
            8: 3,
            9: 2,
            10: 1,
            11: 2,
            12: 1,
            13: 2,
            14: 3,
            15: 4,
            16: 3,
            17: 2,
            18: 1,
            19: 2,
            20: 1,
            21: 2,
            22: 3,
            23: 4,
            24: 7, #Cr has 4s(1)3d(5) configuration
            25: 6,
            26: 5,
            27: 4,
            28: 3,
            29: 2,
            30: 1,
            31: 2,
            32: 3,
            33: 4,
            34: 3,
            35: 2,
            36: 1,
            37: 2,
            38: 1,
        }

    """
    determine the spin of the ground state of an atom or ion
    """
    def get_spinmult(self, symbol, charge=0):
        anumber = qcelemental.periodictable.to_atomic_number(symbol)  
        nele = anumber - charge
        spinmult = 0
        #Neutral atoms
        if charge == 0:
            if anumber <= len(self.mult_dict):
                spinmult = self.mult_dict[nele]
            else:
                raise ValueError("Determination of ion spinmult hasn't been\
                        implemented for elements with electrons occupying\
                        4d orbitals")
        else:
            if anumber <= 20: # H to Ca
                spinmult = self.mult_dict[nele]
            elif anumber <= 30: #Sc - Zn
                if charge == 1:#lose one 4s electron
                    # Assume that electrons are always removed from 4s
                    #TODO: check this against literature
                    ms_4s = 1
                    nele_3d = nele-18-1
                    ms_3d = nele_3d
                    if nele_3d > 5:
                        ms_3d = 10 - nele_3d
                    spinmult = ms_4s + ms_3d + 1
                elif charge >= 2:#lose both 4s electrons
                    n_3d = anumber - 18 - charge
                    if n_3d <= 5:
                        spinmult = n_3d + 1
                    else:
                        spinmult = 10 - n_3d + 1
            elif anumber <=36: # Ga - Kr
                nele_4p = nele-30
                #no more than 4p^25s^2 ocupation
                if charge >= nele_4p-8:
                    spinmult = self.mult_dict[nele]
                else:
                    raise ValueError("Determination of ion spinmult hasn't been\
                        implemented for elements with electrons occupying\
                        4d orbitals")
            else:
                raise ValueError("Determination of ion spinmult hasn't\
                        been implemented for ions with electrons occupying\
                        4d orbitals")
        return spinmult


'''class atomSpinMultDictGen:
    def __init__(self):
        #1s, 1s2s
        mult_dict = {1: 2, 2: 1, 3: 2, 4: 1}
        # 1s2s2p, n(2p) <=3
        for i in range(5, 8):
            mult_dict[i] = i - 4 + 1
        # 1s2s2p, n(2p) > 3
        for i in range(8, 11):
            mult_dict[i] = 10 - i + 1
        # 1s2s2p3s
        mult_dict[11] = 2
        mult_dict[12] = 1
        # 1s2s2p3s3p, n(3p) <=3
        for i in range(13, 16):
            mult_dict[i] = i - 12 + 1
        # 1s2s2p3s3p, n(3p) >3
        for i in range(16, 19):
            mult_dict[i] = 18 - i + 1
        # 1s2s2p3s3p4s
        mult_dict[19] = 2
        mult_dict[20] = 1
        # 1s2s2p3s3p4s3d, n(3d) <= 5
        for i in range(21, 26):
            mult_dict[i] = i - 20 + 1
        # 1s2s2p3s3p4s3d, n(3d) >5
        for i in range(26, 31):
            mult_dict[i] = 30 - i + 1
        # 1s2s2p3s3p4s3d4p, n(4p) <=3
        for i in range(31, 34):
            mult_dict[i] = i - 30 + 1
        # 1s2s2p3s3p4s3d4p, n(4p) >3
        for i in range(34, 37):
            mult_dict[i] = 36 - i + 1
        self.mult_dict = mult_dict

    def getDict(self, nele):
        return self.mult_dict

'''
