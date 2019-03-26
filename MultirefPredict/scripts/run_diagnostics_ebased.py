import MultirefPredict

diagnostics_list = ["B1", "A25PBE"]

input_dict = {"xyzfile": "./geo.xyz",
              "molname": "geo",
              "charge": 0,
              "spinmult": 1,
              "rundir": './',
              "program": 'psi4',
              "record": True}
for diag in diagnostics_list:
    MultirefPredict.diagnostic_factory(diag, **input_dict).computeDiagnostic()
