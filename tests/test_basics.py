#-----------------------------------------------------------------------------
# Copyright (c) 2020, Greg Landrum
# All rights reserved.
#
# The full license is in the LICENSE file, distributed with this software.
#-----------------------------------------------------------------------------

import unittest
import intake
from rdkit import Chem
from pathlib import Path


class TestCase(unittest.TestCase):
    def test_smiles(self):
        for fn in ("data/first_200.tpsa.csv", "data/first_200.tpsa.csv.gz"):
            fn = Path(__file__).parent / fn
            ds = intake.open_smiles(str(fn), header=None, smilesColumn=0)
            df = ds.read()
            self.assertEqual(len(df), 200)
            self.assertEqual(Chem.MolToSmiles(df.iloc[0]['mol']),
                             'CC1=CC(=O)C=CC1=O')
            self.assertEqual((df.iloc[0][1]), 34.14)

    def test_sdf(self):
        for fn in ("data/first_200.props.sdf", "data/first_200.props.sdf.gz"):
            fn = Path(__file__).parent / fn
            ds = intake.open_sdf(str(fn))
            df = ds.read()
            self.assertEqual(len(df), 200)
            self.assertEqual(Chem.MolToSmiles(df.iloc[0]['mol']),
                             'CC1=CC(=O)C=CC1=O')
