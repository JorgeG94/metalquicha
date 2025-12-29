.. _input_files:

Input file formats
======================

.. contents::


This is an example of a JSON input file that uses fragmentation with the MBE method:

.. code-block:: json

    {
        "schema": { "name": "mqc-frag", "version": "1.0" },
        "molecules": [{
            "xyz": "sample_inputs/prism.xyz",
            "fragments": [
              [0,1,2],
              [3,4,5],
              [6,7,8],
              [9,10,11],
              [12,13,14],
              [15,16,17]
            ],
            "fragment_charges": [
                0,
                0,
                0,
                0,
                0,
                0
            ],
            "fragment_multiplicities": [1, 1, 1, 1, 1, 1],
            "molecular_charge": 0,
            "molecular_multiplicity": 1
        }],
      "model": {
        "method": "XTB-GFN1"
      },
      "keywords":{
        "scf": {
          "maxiter": 300,
          "tolerance": 1e-6
        },
        "fragmentation":{
            "method": "MBE",
            "level": 2
        }
      },
      "system": {
        "logger":{
          "level": "Verbose"
        }
      },
      "driver": "Energy"
    }

For simplicity, you can use the "xyz" to provide the path to the geometry file in XYZ format. However,
you can also provide the geometry direclty in the json, for example:

.. code-block:: json

    {
        "schema": { "name": "mqc-frag", "version": "1.0" },
        "molecules": [{
            "geometry": {
                "symbols": ["H", "O", "H"],
                "coordinates": [
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.96],
                    [0.76, 0.0, -0.24]
                ]
            },
            ...
        }],
        ...
    }

Fragments group
================

In the "fragments" field, you can provide a list of fragments, where each fragment is a list of
atom indices (0-based). The "fragment_charges" and "fragment_multiplicities" fields provide the
charge and multiplicity of each fragment, respectively. The "molecular_charge" and "molecular_multiplicity"
fields provide the overall charge and multiplicity of the molecule.

Metalquicha will check that the sum of the fragment charges and multiplicities match the overall
molecular charge and multiplicity.


Model group
============

The "model" field provides the method and basis set to be used in the calculation. The "method"
field can be any method supported by Metalquicha, such as "HF", "DFT", "XTB-GFN1", etc.

Right now, only the XTB methods are supported. So anything that's not XTB-GFN1 or XTB-GFN2 will result
in an error. Support for more methods will be added in the future.

Driver group
============

The "driver" field provides the type of calculation to be performed. The "driver" field can be
"Energy", "Gradient", "Hessian", etc. Right now, only "Energy" and "Gradient" are supported.
