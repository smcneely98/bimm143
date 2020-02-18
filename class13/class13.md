Class 13: Bioinformatics in Drug Discovery and Design
================

## *In silico* docking of drugs to HIV-1 protease

### 1.1 Obtaining and inspecting our input structure

``` r
library(bio3d)
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
hiv <- read.pdb(file.name)
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

> Q1. What is the name of the two non protein resid values in this
> structure? What does resid correspond to and how would you get a
> listing of all reside values in this structure?

### 1.2 Prepare initial protein and ligand input files

``` r
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```

### 1.3 Using tools to setup protein docking input

> Q2. Can you locate the binding site visually? Note that crystal
> structures normally lack hydrogen atoms, why?

> Q3. Look at the charges. Does it make sense (e.g. based on your
> knowledge of the physiochemical properties of amino acids)?

### 1.4 Prepare the ligand

### 1.5 Prepare a docking configuration

## Docking Ligands into HIV-1 Protease

### 2.1 Docking indinavir into HIV-1 protease

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

> Q4. Qualitatively, how good are the docks? Is the crystal binding mode
> reproduced? Is it the best conformation according to AutoDock Vina?

``` r
ori <- read.pdb("ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.649  4.206 11.110 10.529  4.840 10.932 10.993  3.655 10.996 11.222
    ## [11] 10.567 10.372 11.019 11.338  8.390  9.063  8.254  8.978

> Q5. Quantitatively how good are the docks? Is the crystal binding mode
> reproduced within 1Å RMSD for all atoms?

> Q6. How would you determine the RMSD for heavy atoms only (i.e. non
> hydrogen atoms)?

## Exploring the Conformational Dynamics of Proteins

### 3.1 Normal Mode Analysis (NMA)

``` r
library(bio3d)
pdb <- read.pdb("1HEL")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.01 seconds.
    ##  Diagonalizing Hessian...    Done in 0.08 seconds.

``` r
plot(modes, sse=pdb)
```

![](class13_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
mktrj(modes, mode=7, file="nma_7.pdb")
```

> Q7. What are the most flexible portions of HIV-1 protease? Would this
> flexibility likely effect docking calculations to this protein?
