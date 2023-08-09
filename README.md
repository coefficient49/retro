# retro
reverse translation and optimization

installation:

```bash
pip install git+https://github.com/coefficient49/retro.git
```
running basic reverse translation and optimization
```bash
retro -f test.csv
```

using it in a script

```python
from retro.main import *

peptides = ["ACTYACTYACTAYCTYA","ACTYACTYXCTXYCTYA"]
enzymes = ["BsaI","EcoRV"]
enzymes = enzymes_cls = RestrictionBatch([get_enzyme_class_from_str(x) for x in enzymes])

outputs = run_all(peptides,enzyme_filter=enzymes_cls,hamming_check=False)

outputs["final_df_output"]

```
