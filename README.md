# RDKit-based chemistry support for intake

This adds support for a couple of chemistry datatypes to [intake](https://intake.readthedocs.io/en/latest/index.html)

## To install

You'll need to have both intake and the rdkit installed in your environment:
```
conda install -c conda-forge intake rdkit`
```

And then you can `pip install` intake-rdkit directly from this repo:
```
python -m pip install git+git://github.com/greglandrum/intake-rdkit.git
```

## Usage

### Directly loading SDF or SMILES files

Here's reading a compressed CSV file (of course the file doesn't have to be compressed):
```
>>> ds = intake.open_smiles('./files/CHEMBL1821_Ki_set.csv.gz',smilesColumn='canonical_smiles')
>>> df = ds.read()
>>> df.head()
                                                mol      chembl_id  ... pchembl_value doc_id
0  <rdkit.Chem.rdchem.Mol object at 0x7fa5df5e6490>  CHEMBL1794855  ...           NaN   6491
1  <rdkit.Chem.rdchem.Mol object at 0x7fa5df6bc6c0>  CHEMBL2112955  ...          7.85  16845
2  <rdkit.Chem.rdchem.Mol object at 0x7fa5df6bc760>  CHEMBL2112957  ...          8.64  16845
3  <rdkit.Chem.rdchem.Mol object at 0x7fa5df5e6620>   CHEMBL369062  ...          7.22  16845
4  <rdkit.Chem.rdchem.Mol object at 0x7fa5df5e6580>  CHEMBL2112961  ...          7.37  16845

[5 rows x 7 columns]
```

And here's an SDF (you can read `.sdf.gz` too):
```
>>> ds = intake.open_sdf('./files/jm200186n.sdf')
>>> df = ds.read()
>>> df.head()
                                                mol  scaffold
0  <rdkit.Chem.rdchem.Mol object at 0x7fa5df5763f0>       1.0
1  <rdkit.Chem.rdchem.Mol object at 0x7fa5df576080>       NaN
2  <rdkit.Chem.rdchem.Mol object at 0x7fa5df5765d0>       NaN
3  <rdkit.Chem.rdchem.Mol object at 0x7fa5df576760>       NaN
4  <rdkit.Chem.rdchem.Mol object at 0x7fa5df5760d0>       NaN
```

Note that calling `ds.read()` parses all the molecules in the dataset and reads
them into a pandas DataFrame, so be careful with big data files.

### Working with a data catalog

This is more interesting.

Here's the catalog I'm working with, which I have saved in a file called `literature.yaml`:
```
metadata:
  version: 1
  creator: 
    name: greg landrum
    email: greg.landrum@t5informatics.com

  summary: |
    Collection of datasets pulled from the literature

sources:
  cdk2_project:
    description: screening results and synthesized compounds for a CDK2 project.
    args:
      filename: '{{ CATALOG_DIR }}/files/jm020472j_2.csv.gz'
      smilesColumn: Smiles
    metadata:
      journal_url: https://pubs.acs.org/doi/10.1021/jm020472j
      additional_notes: |
        The scaffold column is a manually assignment to chemical series.
        
        The sourcepool column indicates whether the compound comes from the
        screening deck (divscreen) or was synthesized for the project
        (synscreen)
    driver: intake_rdkit.smiles.SmilesSource

  platinum_2017:
    description: Platinum 2017 set for testing conformation generators
    args:
      filename: '{{ CATALOG_DIR }}/files/platinum_dataset_2017_01.sdf.gz'
    metadata:
      journal_url: https://pubs.acs.org/doi/10.1021/acs.jcim.6b00613
      additional_notes: |
        This is the 2017 update of the platinum set as described here: 
        https://pubs.acs.org/doi/10.1021/acs.jcim.7b00505
    driver: intake_rdkit.sdf.SDFSource
```

And here's how I work with it:
```
>>> cat = intake.open_catalog('./literature.yaml')
>>> cat.metadata
{'version': 1, 'creator': {'name': 'greg landrum', 'email': 'greg.landrum@t5informatics.com'}, 'summary': 'Collection of datasets pulled from the literature\n'}
>>> for entry in cat:
...     print(entry)
... 
cdk2_project
platinum_2017
>>> cat.platinum_2017
sources:
  platinum_2017:
    args:
      filename: /scratch/cheminformatics_datasets/.//files/platinum_dataset_2017_01.sdf.gz
    description: Platinum 2017 set for testing conformation generators
    driver: intake_rdkit.sdf.SDFSource
    metadata:
      additional_notes: 'This is the 2017 update of the platinum set as described
        here: https://pubs.acs.org/doi/10.1021/acs.jcim.7b00505'
      catalog_dir: /scratch/cheminformatics_datasets/./
      journal_url: https://pubs.acs.org/doi/10.1021/acs.jcim.6b00613
>>> df = cat.platinum_2017.read()
>>> len(df)
4548
>>> df.head()
                                                mol
0  <rdkit.Chem.rdchem.Mol object at 0x7fa5df6205d0>
1  <rdkit.Chem.rdchem.Mol object at 0x7fa5df620a80>
2  <rdkit.Chem.rdchem.Mol object at 0x7fa5df620670>
3  <rdkit.Chem.rdchem.Mol object at 0x7fa5df620210>
4  <rdkit.Chem.rdchem.Mol object at 0x7fa5f8c69e40>
```