### Metagen pipeline v2

#### HOW-TO

Test dataset

```
snakemake -s main.snk --config-file config.json --cores 8
```

#### Test configuration file

`test/config_sim`:

**Simple**
Ref: GCF_011881725.1: Escherichia coli_E 
* GCA_000012825.1 (10X) : Phocaeicola vulgatus (rep: itself)
* GCA_000435415.1 (15X) : Phocaeicola merdigallinarum (rep: GCA_902362595)
* GCF_011881725.1 (20X) : Escherichia coli_E (same than ref, rep: itself)
* GCF_000196555.1 (30X) : Bifidobacterium longum (rep: itself)
* GCA_000223435.2 (15X) : Klebsiella pneumonia (rep: GCF_000742135)

**Fergu**
Ref: GCF_011881725.1: Escherichia coli_E 
* GCA_000012825.1 (10X) : Phocaeicola vulgatus (rep: itself)
* GCA_000435415.1 (15X) : Phocaeicola merdigallinarum (rep: GCA_902362595)
* GCF_011881725.1 (20X) : Escherichia coli_E (same than ref, rep: itself)
* GCF_000196555.1 (30X) : Bifidobacterium longum (rep: itself)
* GCA_000223435.2 (15X) : Klebsiella pneumonia (rep: GCF_000742135)
* GCF_000026225.1 (20X) : Escherichia fergusonii (rep: itself)

**Fergu Alberti**
Ref: GCF_000026225.1: Escherichia fergusonii 
* GCA_000012825.1 (10X) : Phocaeicola vulgatus (rep: itself)
* GCA_000435415.1 (15X) : Phocaeicola merdigallinarum (rep: GCA_902362595)
* GCF_000196555.1 (30X) : Bifidobacterium longum (rep: itself)
* GCF_000026225.1 (20X) : Escherichia fergusonii (rep: itself)
* GCF_000759775.1 (30X) : GCF_000759775.1 (rep: itself)

**EStrain_simple** (to check result with next one, E. coli reference change)
Ref: GCF_003697165.2: Escherichia coli
* GCA_000012825.1 (10X) : Phocaeicola vulgatus (rep: itself)
* GCA_000435415.1 (15X) : Phocaeicola merdigallinarum (rep: GCA_902362595)
* GCF_900635635.1 (20X) : Escherichia coli (rep: RGCF_003697165.2)
* GCF_000196555.1 (30X) : Bifidobacterium longum (rep: itself)
* GCA_000223435.2 (15X) : Klebsiella pneumonia (rep: GCF_000742135)
* GCF_000026225.1 (20X) : Escherichia fergusonii (rep: itself)

**EStrain**
Ref: GCF_003697165.2: Escherichia coli
* GCA_000012825.1 (10X) : Phocaeicola vulgatus (rep: itself)
* GCA_000435415.1 (15X) : Phocaeicola merdigallinarum (rep: GCA_902362595)
* GCF_000196555.1 (30X) : Bifidobacterium longum (rep: itself)
* GCA_000223435.2 (15X) : Klebsiella pneumonia (rep: GCF_000742135)
* GCA_003864075.1 (10X) : Escherichia albertii (rep: GCF_000759775.1)
* GCA_902385815.1 (10X) : Escherichia marmotae (rep: GCF_002900365.1)
* GCA_902385855.1 (10X) : Escherichia ruysiae (rep: GCF_902498915.1)
* GCF_004323655.1 (10X) : Escherichia sp004211955 (rep: GCF_004211955.1)
* GCF_004745115.1 (10X) : Escherichia sp005843885 (rep: GCF_005843885.1)
* GCF_013171325.1 (10X) : Escherichia fergusonii (rep: GCF_000026225.1)
* GCF_014836715.1 (10X) : Escherichia sp001660175 (rep: GCF_001660175.1)
* GCF_900635635.1 (10X) : Escherichia coli (rep: GCF_003697165.2)


### TODO

* See what happen if we use two strains from the same species or two similar species
	* Check the sourmash results first
* Should we stick to only one similar genome based on ANI search (ref genome)?
* Use a zymo to test