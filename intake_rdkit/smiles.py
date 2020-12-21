#-----------------------------------------------------------------------------
# Copyright (c) 2020, Greg Landrum
# All rights reserved.
#
# The full license is in the LICENSE file, distributed with this software.
#-----------------------------------------------------------------------------

from rdkit import Chem
import pandas as pd
import intake
import numpy as np
import gzip

class SmilesSource(intake.source.base.DataSource):
    " Inspired by example here: https://github.com/intake/intake-examples/blob/master/tutorial/data_engineer.ipynb "
    container = 'dataframe'
    name = 'smiles'
    version = '0.1.0'
    partition_access = True
    suppl_factory = Chem.SmilesMolSupplier
    
    def __init__(self, filename, metadata=None, nPartitions=1, smilesColumn='smiles', 
        useSmarts=False, **kwargs):
        super(SmilesSource, self).__init__(
            metadata=metadata
        )
        self._fname = filename
        self._smilesColumn = smilesColumn
        self._useSmarts = useSmarts
        self._supplkwargs = kwargs
        # this is, hopefully, for future use:
        self._nPartitions = nPartitions

    def _get_schema(self):
        " reads property names from the first molecule and uses those "
        # FIX handle error if the first molecule doesn't exist
        df = pd.read_csv(self._fname,nrows=1,header=0,**self._supplkwargs)
        dt = {'mol':'O'}
        for pn,pv in zip(df.columns,df.iloc[0]):
            dt[pn] = np.dtype(type(pv)).name
        self._colNames = list(df.columns)
        return intake.source.base.Schema(
            datashape=None,
            dtype=dt,
            shape=(None, len(dt)),
            npartitions=self._nPartitions,
            extra_metadata=self._supplkwargs
        )

    def _get_partition(self, i):
        df = pd.read_csv(self._fname,**self._supplkwargs)
        partSize = len(df)//self._nPartitions
        start = i*partSize
        if i==self._nPartitions-1:
            end = len(df)
        else:
            end = start + partSize
        res = [] 
        for idx in range(start,end):
            row = df.iloc[idx]
            if not self._useSmarts:
                m = Chem.MolFromSmiles(row[self._smilesColumn])
            else:
                m = Chem.MolFromSmarts(row[self._smilesColumn])
            if m is not None:
                nrow = {'mol':m}
                for k in self._colNames:
                    nrow[k] = row[k]
                res.append(nrow)
        return pd.DataFrame(res)

    def read(self):
        self._load_metadata()
        return pd.concat([self.read_partition(i) for i in range(self.npartitions)])
        
    def _close(self):
        # close any files, sockets, etc
        self.suppl=None
