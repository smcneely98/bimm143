---
title: "Class 13: Bioinformatics in Drug Discovery and Design"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## *In silico* docking of drugs to HIV-1 protease

### 1.1 Obtaining and inspecting our input structure

```{r}
library(bio3d)
file.name <- get.pdb("1hsg")
```

```{r}
hiv <- read.pdb(file.name)
hiv
```

> Q1. What is the name of the two non protein resid values in this structure? What does resid
correspond to and how would you get a listing of all reside values in this structure?



### 1.2 Prepare initial protein and ligand input files

```{r}
prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```

### 1.3 Using tools to setup protein docking input

> Q2. Can you locate the binding site visually? Note that crystal structures normally lack hydrogen atoms, why?



> Q3.  Look at the charges. Does it make sense (e.g. based on your knowledge of the physiochemical properties of amino acids)?



### 1.4 Prepare the ligand



### 1.5 Prepare a docking configuration



## Docking Ligands into HIV-1 Protease


## 