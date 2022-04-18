# Primer-Designer

This script will calculate the optimal PCR primers for a given gene.

The script finds the DNA sequence of the target gene using the [Ensembl REST API](https://rest.ensembl.org), finds overlapping sections between any splice variants, then initialises [`primer3`](https://github.com/primer3-org/primer3) to find the best primer pairs according to a set of parameters.

## Prerequisites

You must have `python3` installed. You will need to install any other outstanding requirements:

| Command | Installation |
| --- | --- |
| `primer3` | Install via [brew](https://brew.sh) by running `brew install primer3` |

## Usage

Download with:
```
git clone https://github.com/phenotypic/Primer-Designer.git
pip3 install -r requirements.txt
```

Run from the same directory with:
```
python3 designer.py
```

Here are some flags you can add:

| Flag | Description |
| --- | --- |
| `-s <species>` | Species: Define a species (script will prompt you otherwise) |
| `-g <gene>` | Gene: Define a target gene (script will prompt you otherwise) |
| `-c <method>` | Interface: Set a thermodynamic parameters path (script will use default otherwise) |
| `-v` | Verbose: Activate verbose logging |

After running the script and selecting a target species and gene, the overlapping regions of any splice variants will be found.

Next, you will be asked if you want to adjust the default primer search parameters:

```
+-------+--------------------+---------+
| Index |     Parameter      |  Value  |
+-------+--------------------+---------+
|   1   |  Primer min size   |    20   |
|   2   |  Primer opt size   |    23   |
|   3   |  Primer max size   |    26   |
|   4   |   Primer min Tm    |    56   |
|   5   |   Primer opt Tm    |    60   |
|   6   |   Primer max Tm    |    64   |
|   7   | Product size range | 300-500 |
+-------+--------------------+---------+
```

Finally, the script will generate pairs of primers. Primer pair 1 is usually the most optimal:

```
Primer pair 1
---------------------
Gene name: wnt10a
Species: zebrafish
Product size: 458
Any compl   : 0.00
End compl   : 0.00
+-----------+--------+--------+--------+--------------+---------+---------------+-------------------------+
| Direction | Length |   Tm   |  GC%   | Self-binding | Hairpin | End Stability |         Sequence        |
+-----------+--------+--------+--------+--------------+---------+---------------+-------------------------+
|    Left   |   23   | 59.995 | 47.826 |     0.00     |   0.00  |     3.2800    | CCTCTTCTGTTCTTTGGCTTTGG |
|   Right   |   23   | 60.121 | 47.826 |     0.00     |   0.00  |     3.1800    | GACGGTTGAGCTTGATTCTGAAC |
+-----------+--------+--------+--------+--------------+---------+---------------+-------------------------+
```

## Notes

- If your target gene has lots of splice varaints (e.g. [DARS2](https://ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000117593;r=1:173824653-173858808)), it is possible that there will not be any overlapping sequences across all varaints. In this case, you should manually select the splice varaints which are most prevalent.
- The script has been thoroughly tested with these `zebrafish` genes: 
  - `wnt10a` (direction: 1, variants: multiple)
  - `lamb1a` (direction: -1, variants: multiple)
  - `slc7a5` (direction: -1, variants: single)

## To-do

- Use parameters file as default values (not coded into script)
- BLAST primer check
- Allow command arguments for primer parameters
- Allow input for left or right primer
- Save info and sequence files locally for offline checks
