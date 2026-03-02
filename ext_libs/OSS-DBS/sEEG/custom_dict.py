#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 19:10:29 2025

@author: forel
"""

oss_params = {
  "BrainRegion": {
    "Center": {
      "x[mm]": 8.04146971649281,
      "y[mm]": -14.364058770144673,
      "z[mm]": 0.9143724523894229
    },
    "Dimension": {
      "x[mm]": 54.383297871632166,
      "y[mm]": 53.38734295967949,
      "z[mm]": 60.67471786493622
    },
    "Shape": "Ellipsoid"
  },
  "Electrodes": [
    {
      "Name": "BostonScientificVerciseDirectedCustom",
      "CustomParameters": {
        "tip_length": 1.5,
        "contact_length": 1.5,
        "contact_spacing": 0.5,
        "lead_diameter": 1.3,
        "total_length": 450.0
      },
      "Rotation[Degrees]": -30.333156848528787,
      "Direction": {
        "x[mm]": 0.36446958383539824,
        "y[mm]": 0.2816563042206489,
        "z[mm]": 0.8875988106975004
      },
      "TipPosition": {
        "x[mm]": 7.76811752861626,
        "y[mm]": -14.57530099831016,
        "z[mm]": 0.24867334436629762
      },
      "EncapsulationLayer": {
        "Thickness[mm]": 0.0,
        "Material": "Gray matter",
        "DielectricModel": "ColeCole4",
        "DielectricParameters": None,
        "MaxMeshSize": 0.1
      },
      "Contacts": [
        {
          "Contact_ID": 1,
          "Active": False,
          "Current[A]": 0.0,
          "Voltage[V]": 0.0,
          "Floating": True,
          "SurfaceImpedance[Ohmm]": {
            "real": 0.0,
            "imag": 0.0
          },
          "MaxMeshSize": 1000000.0,
          "MaxMeshSizeEdge": 1000000.0
        },
        {
          "Contact_ID": 2,
          "Active": False,
          "Current[A]": 0.0,
          "Voltage[V]": 0.0,
          "Floating": True,
          "SurfaceImpedance[Ohmm]": {
            "real": 0.0,
            "imag": 0.0
          },
          "MaxMeshSize": 1000000.0,
          "MaxMeshSizeEdge": 1000000.0
        },
        {
          "Contact_ID": 3,
          "Active": False,
          "Current[A]": 0.0,
          "Voltage[V]": 0.0,
          "Floating": True,
          "SurfaceImpedance[Ohmm]": {
            "real": 0.0,
            "imag": 0.0
          },
          "MaxMeshSize": 1000000.0,
          "MaxMeshSizeEdge": 1000000.0
        },
        {
          "Contact_ID": 4,
          "Active": False,
          "Current[A]": 0.0,
          "Voltage[V]": 0.0,
          "Floating": True,
          "SurfaceImpedance[Ohmm]": {
            "real": 0.0,
            "imag": 0.0
          },
          "MaxMeshSize": 1000000.0,
          "MaxMeshSizeEdge": 1000000.0
        },
        {
          "Contact_ID": 5,
          "Active": False,
          "Current[A]": 0.0,
          "Voltage[V]": 0.0,
          "Floating": True,
          "SurfaceImpedance[Ohmm]": {
            "real": 0.0,
            "imag": 0.0
          },
          "MaxMeshSize": 1000000.0,
          "MaxMeshSizeEdge": 1000000.0
        },
        {
          "Contact_ID": 6,
          "Active": False,
          "Current[A]": 0.0,
          "Voltage[V]": 0.0,
          "Floating": True,
          "SurfaceImpedance[Ohmm]": {
            "real": 0.0,
            "imag": 0.0
          },
          "MaxMeshSize": 1000000.0,
          "MaxMeshSizeEdge": 1000000.0
        },
        {
          "Contact_ID": 7,
          "Active": False,
          "Current[A]": 0.0,
          "Voltage[V]": 0.0,
          "Floating": True,
          "SurfaceImpedance[Ohmm]": {
            "real": 0.0,
            "imag": 0.0
          },
          "MaxMeshSize": 1000000.0,
          "MaxMeshSizeEdge": 1000000.0
        },
        {
          "Contact_ID": 8,
          "Active": False,
          "Current[A]": -0.005,
          "Voltage[V]": 0.0,
          "Floating": True,
          "MaxMeshSizeEdge": 0.08168140899333462,
          "SurfaceImpedance[Ohmm]": {
            "real": 0.0,
            "imag": 0.0
          },
          "MaxMeshSize": 1000000.0
        }
      ]
    }
  ],
  "Surfaces": [
    {
      "Name": "BrainSurface",
      "Active": True,
      "Current[A]": 0.005,
      "Voltage[V]": 0.0
    }
  ],
  "MaterialDistribution": {
    "MRIPath": "/home/interscan/Documents/data/RajamaniTWEEED/derivatives/leaddbs/sub-TWEED07/stimulations/native/StimFitOSS/segmask.nii",
    "MRIMapping": {
      "Blood": 4,
      "Gray matter": 1,
      "White matter": 2,
      "CSF": 3,
      "Unknown": 0
    },
    "DiffusionTensorActive": False,
    "DTIPath": ""
  },
  "DielectricModel": {
    "Type": "ColeCole4",
    "CustomParameters": None
  },
  "Mesh": {
    "LoadMesh": False,
    "LoadPath": "",
    "MeshingHypothesis": {
      "Type": "Default",
      "MaxMeshSize": 1000000.0,
      "MeshSizeFilename": ""
    },
    "HPRefinement": {
      "Active": False,
      "Levels": 2,
      "Factor": 0.125
    },
    "AdaptiveMeshRefinement": {
      "Active": False,
      "MaxIterations": 10,
      "ErrorTolerance": 0.1
    },
    "MaterialRefinementSteps": 1,
    "MeshSize": {
      "Edges": {},
      "Faces": {},
      "Volumes": {}
    },
    "SaveMesh": False,
    "SavePath": "mesh"
  },
  "EQSMode": False,
  "FEMOrder": 2,
  "StimulationSignal": {
    "Type": "Multisine",
    "ListOfFrequencies": [
      10000
    ],
    "Frequency[Hz]": 130.0,
    "PulseWidth[us]": 60.0,
    "PulseTopWidth[us]": 0.0,
    "CounterPulseWidth[us]": 0.0,
    "InterPulseWidth[us]": 0.0,
    "SpectrumMode": "FullSpectrum",
    "CounterAmplitude": 1.0,
    "CutoffFrequency": 1000000.0,
    "CurrentControlled": True
  },
  "Solver": {
    "Type": "CG",
    "Preconditioner": "local",
    "PreconditionerKwargs": {},
    "MaximumSteps": 500,
    "Precision": 1e-10
  },
  "PointModel": {
    "Pathway": {
      "Active": False,
      "FileName": "/home/interscan/Documents/data/RajamaniTWEEED/derivatives/leaddbs/sub-TWEED07/stimulations/native/StimFitOSS/OSS_sim_files_rh/Allocated_axons.h5",
      "ExportField": False
    },
    "Lattice": {
      "Active": True,
      "Center": {
        "x[mm]": 10.233118652308892,
        "y[mm]": -12.670387290304925,
        "z[mm]": 6.251731384857534
      },
      "Shape": {
        "x": 71,
        "y": 71,
        "z": 71
      },
      "Direction": {
        "x[mm]": 0,
        "y[mm]": 0,
        "z[mm]": 1
      },
      "PointDistance[mm]": 0.33,
      "CollapseVTA": False,
      "ExportField": False
    },
    "VoxelLattice": {
      "Active": False,
      "Shape": {
        "x": 31,
        "y": 31,
        "z": 41
      },
      "TimeDomain": False,
      "ExportField": False
    }
  },
  "OutputPath": "/home/interscan/Documents/data/RajamaniTWEEED/derivatives/leaddbs/sub-TWEED07/stimulations/native/StimFitOSS/OSS_sim_files_rh/Results",
  "ComputeImpedance": False,
  "ExportVTK": True,
  "ExportFrequency": None,
  "ExportElectrode": False,
  "ModelSide": 0,
  "CalcAxonActivation": False,
  "ActivationThresholdVTA[V-per-m]": 200.0,
  "FailFlag": "rh",
  "OutOfCore": False,
  "PathwayFile": None,
  "StimSets": {
    "Active": False,
    "StimSetsFile": "/home/interscan/Documents/data/RajamaniTWEEED/derivatives/leaddbs/sub-TWEED07/stimulations/native/StimFitOSS/OSS_sim_files_rh/Current_protocols_0.csv"
  },
  "StimulationFolder": "/home/interscan",
  "TruncateAfterActivePartRatio": None
}