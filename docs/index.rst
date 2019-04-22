.. MultirefPredict documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MultirefPredict's documentation!
=========================================================
*Automated workflow to predict multireference character of molecules in quantum chemistry calculation*

Program Execution
-----------------

A simple example of MultirefPredict's capabilities of as follows:

.. code:: python

   >>> import MultirefPredict
   >>> import qcelemental
   >>> mol = qcelemental.models.Molecule.from_data("""
   O                 0.000000000000     0.000000000000    -0.068516245955
   H                 0.000000000000    -0.790689888800     0.543701278274
   H                 0.000000000000     0.790689888800     0.543701278274
   """)

Here we just built a molecule object in qcelemental format. Various multireference diagnostics can be calculated for this molecule with the ``computeDiagnostic`` syntax of the diagnostic object, which is initialized with the ``diagnostic_factory`` function.

.. code:: python

   >>> b1 = MultirefPredict.diagnostic_factory("B1",molecule=mol, molname="water", record=False).computeDiagnostic()

``b1`` is the returned diagnostic value

.. code:: python

   >>> b1
   0.006334860365228678

The I/O of the backend quantum chemistry package is totally hidden from the user. If one would like to get detailed information about the quantum chemistry calculation, it can be easily achieved by setting the ``record`` keyword to ``True``.

.. code:: python

   >>> b1 = MultirefPredict.diagnostic_factory("B1",molecule=mol, molname="water", record=True).computeDiagnostic()

Then the jason files recording the quantum chemistry calculations associated with this diagnostic calculation are automatically dumped in the working directory

.. code:: bash
  
   >>> ls 
   B1_water_b1lyp_H.json
   B1_water_b1lyp_O.json
   B1_water_b1lyp_whole.json
   B1_water_blyp_H.json
   B1_water_blyp_O.json
   B1_water_blyp_whole.json

These files are in the qcelemental result format. For example:

.. code:: python
  
   >>> with open('B1_water_b1lyp_H.json') as f:
            data = json.load(f)
                print(json.dumps(data,indent=4))
   {
       "molecule": {
               "symbols": [
                   "H"
               ],
               "geometry": [
                    0.0,
                    0.0,
                    0.0
                ],
                "molecular_charge": 0.0,
                "molecular_multiplicity": 2
        },
        "driver": "energy",
        "model": {
             "method": "b1lyp",
              "basis": "6-31g"
         },
         ...
   }

Backends
--------

Currently available compute backends for single results are as follow:

- `Psi4 <http://www.psicode.org/>`_
- `TeraChem <http://www.petachem.com/products.html>`_

All backends are driven by `QCEngine <https://github.com/MolSSI/QCEngine>`_

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   examples
