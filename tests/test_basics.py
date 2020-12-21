#-----------------------------------------------------------------------------
# Copyright (c) 2020, Greg Landrum
# All rights reserved.
#
# The full license is in the LICENSE file, distributed with this software.
#-----------------------------------------------------------------------------

import unittest
import intake
from rdkit import Chem


class TestCase(unittest.TestCase):
    def test_smiles(self):
        ds = intake.open_smiles('./data/first_200.tpsa.csv',
                                header=None,
                                smilesColumn=0)
        df = ds.read()
        self.assertEqual(len(df), 200)
        self.assertEqual(Chem.MolToSmiles(df.iloc[0]['mol']),
                         'CC1=CC(=O)C=CC1=O')
        self.assertEqual((df.iloc[0][1]), 34.14)

        ds = intake.open_smiles('./data/first_200.tpsa.csv.gz',
                                header=None,
                                smilesColumn=0)
        df = ds.read()
        self.assertEqual(len(df), 200)
        self.assertEqual(Chem.MolToSmiles(df.iloc[0]['mol']),
                         'CC1=CC(=O)C=CC1=O')
        self.assertEqual((df.iloc[0][1]), 34.14)

    def test_sdf(self):
        ds = intake.open_sdf('./data/first_200.props.sdf')
        df = ds.read()
        self.assertEqual(len(df), 200)
        self.assertEqual(Chem.MolToSmiles(df.iloc[0]['mol']),
                         'CC1=CC(=O)C=CC1=O')

        ds = intake.open_sdf('./data/first_200.props.sdf.gz')
        df = ds.read()
        self.assertEqual(len(df), 200)
        self.assertEqual(Chem.MolToSmiles(df.iloc[0]['mol']),
                         'CC1=CC(=O)C=CC1=O')
