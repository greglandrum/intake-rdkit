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
class SDFSource(intake.source.base.DataSource):
    " Inspired by example here: https://github.com/intake/intake-examples/blob/master/tutorial/data_engineer.ipynb "
    container = 'dataframe'
    name = 'sdf'
    version = '0.1.0'
    partition_access = True
    suppl_factory = Chem.SDMolSupplier
    
    def __init__(self, filename, metadata=None, nPartitions=1, **kwargs):
        super(SDFSource, self).__init__(
            metadata=metadata
        )
        self._fname = filename
        if filename.endswith('.gz'):
            self.suppl_factory = lambda fname,**kwargs:Chem.ForwardSDMolSupplier(gzip.open(fname),**kwargs)
            if nPartitions != 1:
                raise ValueError('can only read compressed files with a single partition')
        self._supplkwargs = kwargs
        # this is, hopefully, for future use:
        self._nPartitions = nPartitions

    def _get_schema(self):
        " reads property names from the first molecule and uses those "
        # FIX handle error if the first molecule doesn't exist
        suppl = self.suppl_factory(self._fname,**self._supplkwargs)
        m0 = next(suppl)
        pd = m0.GetPropsAsDict()
        dt = {'mol':'O'}
        for pn,pv in pd.items():
            dt[pn] = np.dtype(type(pv)).name
        return intake.source.base.Schema(
            datashape=None,
            dtype=dt,
            shape=(None, len(dt)),
            npartitions=self._nPartitions,
            extra_metadata=self._supplkwargs
        )

    def _get_partition(self, i):
        suppl = self.suppl_factory(self._fname,**self._supplkwargs)
        res = []
        if self._nPartitions>1:
            partSize = len(suppl)//self._nPartitions
            start = i*partSize
            if i==self._nPartitions-1:
                end = len(suppl)
            else:
                end = start + partSize
                
            for idx in range(start,end):
                m = suppl[idx]
                if m is not None:
                    row = {'mol':m}
                    row.update(m.GetPropsAsDict())
                    res.append(row)
        else:
            for m in suppl:
                if m is not None:
                    row = {'mol':m}
                    row.update(m.GetPropsAsDict())
                    res.append(row)
        return pd.DataFrame(res)

    def read(self):
        self._load_metadata()
        return pd.concat([self.read_partition(i) for i in range(self.npartitions)])
        
    def _close(self):
        # close any files, sockets, etc
        self.suppl=None
