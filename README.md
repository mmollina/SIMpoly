# Introduction

This vignette demonstrates how to use **SIMpoly** to simulate multiple
biparental populations and how you can use **mappoly2** to construct a
consensus genetic map. The workflow includes:

1.  Defining simulation parameters and pedigree.
2.  Simulating multiple chromosomes using **SIMpoly**.
3.  Constructing linkage maps using **mappoly2**.
4.  Integrating maps to generate a consensus map.

> **Note**: The consensus map generation requires **mappoly2**. If you
> only have **SIMpoly**, the first simulation steps will work, but the
> final consensus map will not be generated.

------------------------------------------------------------------------

# Step 1: Installing and Loading Packages

    library(SIMpoly)
    library(mappoly2)  # for consensus mapping

------------------------------------------------------------------------

# Step 2: Define Simulation Parameters and Pedigree

    # Define ploidy for each parent.
    ploidy.vec <- c(4, 2, 4, 2, 4, 4)  # for six parents
    names(ploidy.vec) <- c("P1", "P2", "P3", "P4", "P5", "P6")

    # Define parental crosses.
    parents.mat <- matrix(c("P1", "P2",
                            "P1", "P3",
                            "P2", "P2",
                            "P3", "P4",
                            "P4", "P5",
                            "P5", "P6",
                            "P5", "P2"),
                          ncol = 2, byrow = TRUE)

    # Define number of offspring per cross.
    n_offspring <- c(200, 100, 200, 100, 300, 200, 100)

    # Define allele set for each founder.
    alleles <- list(P1 = 0:1, P2 = 0:1, P3 = 0:1, P4 = 0:1, P5 = 0:1, P6 = 0:1)

    # Generate pedigree.
    pedigree <- simulate_pedigree(ploidy.vec, parents.mat, n_offspring)
    head(pedigree, n = 10)
    #>      parent_1 parent_2 individual cross
    #> 1        <NA>     <NA>         P1 NAxNA
    #> 2        <NA>     <NA>         P2 NAxNA
    #> 3        <NA>     <NA>         P3 NAxNA
    #> 4        <NA>     <NA>         P4 NAxNA
    #> 5        <NA>     <NA>         P5 NAxNA
    #> 6        <NA>     <NA>         P6 NAxNA
    #> P1.1       P1       P2    P1xP2_1 P1xP2
    #> P1.2       P1       P2    P1xP2_2 P1xP2
    #> P1.3       P1       P2    P1xP2_3 P1xP2
    #> P1.4       P1       P2    P1xP2_4 P1xP2

------------------------------------------------------------------------

# Step 3: Simulate Multiparental Data

    n.chr <- 3  # Number of chromosomes
    map.len <- c(100, 120, 90)  # Chromosome lengths
    n_mrk <- c(2000, 1500, 2000)  # Number of markers per chromosome

    # Simulate genetic data
    result <- simulate_multiparental_data(n.chr, map.len, pedigree, ploidy.vec, n_mrk, alleles,  missing = 0.1, p = .3 , rho = .7, seed = 1)

    # Extract results
    wide_df <- result$wide_df
    all_dat <- result$dat
    parent_homologs <- result$parent_homologs[[1]]
    plot_parent_phase_gg(parent_homologs)

![](vignette_simulation_files/figure-markdown_strict/simulate-data-1.png)

    plot_correlated_sets(result$counts)

![](vignette_simulation_files/figure-markdown_strict/simulate-data-2.png)

------------------------------------------------------------------------

# Step 4: Construct Linkage Maps with mappoly2

    MAPs <- vector("list", nrow(parents.mat))

    for(i in 1:nrow(parents.mat)){
      p1 <- parents.mat[i,1]
      p2 <- parents.mat[i,2]
      dat <- mappoly2:::table_to_mappoly(all_dat[[i]], ploidy.vec[p1], ploidy.vec[p2], p1, p2)
      
      # Filtering steps
      dat <- filter_markers_by_chisq_pval(dat, plot = FALSE)
      dat <- filter_markers_by_missing_rate(dat, mrk.thresh = .1, plot = FALSE)
      dat <- filter_individuals_by_missing_rate(dat, ind.thresh = .1, plot = FALSE)
      dat <- pairwise_rf(dat, mrk.scope = "all", ncpus = 8)
      plot(dat)
      dat <- rf_filter(dat, probs = c(0, 1), diagnostic.plot = FALSE)
      
      # Group and order markers
      g <- group(dat, expected.groups = n.chr, comp.mat = TRUE, inter = FALSE)
      plot(g)
      s <- make_sequence(g, ch = list(1,2,3))
      s <- order_sequence(s, type = "genome")
      s <- pairwise_phasing(s, type = "genome", parent = "p1p2")
      s <- mapping(s, type = "genome", parent = "p1p2", ncpus = n.chr)
      s <- calc_haplotypes(s, type = "genome", ncpus = n.chr)
      
      MAPs[[i]] <- s
    }
    plot_map(MAPs[[1]], lg = 1, type = "genome")

------------------------------------------------------------------------

# Step 5: Construct a Consensus Map

    # Visualize multiple maps
    plot_multi_map(MAPs)

    # Integrate maps
    w <- prepare_to_integrate(MAPs, type = "genome")
    plot(w)

    # Estimate consensus map
    x <- estimate_consensus_map(w, ncpus = 8)
    plot(x)
    plot(x, only.consensus = TRUE, col = mappoly::mp_pallet3(3))

    # Compute haplotypes
    system.time(x <- calc_consensus_haplo(x, ncpus = 8))

    # Visualize consensus haplotypes
    plot_consensus_haplo(x, lg = 1, ind = "P1xP2_1")
    plot_consensus_haplo(x, lg = 1, ind = "P5xP6_1")

------------------------------------------------------------------------

# Summary

This vignette demonstrated how to: 1. Define and simulate multiparental
populations using **SIMpoly**. 2. Construct linkage maps using
**mappoly2**. 3. Generate and visualize a consensus map.
=======
&#x20;&#x20;

# **SIMpoly - Simulation of Genetic Marker Data in Pedigreed Populations**

### *Marcelo Mollinari*

## 🚀 **Introduction**

This vignette demonstrates how to use **SIMpoly** to simulate multiple biparental populations and how to leverage **mappoly2** to construct a consensus genetic map. The workflow includes:

1. **Installing from GitHub**
2. **Defining simulation parameters and pedigree.**
3. **Simulating multiple chromosomes** using **SIMpoly**.
4. **Constructing linkage maps** using **mappoly2**.
5. **Integrating maps** to generate a consensus map.

> ⚠️ **Note:** The consensus map generation requires **mappoly2**. If you only have **SIMpoly**, the first simulation steps will work, but the final consensus map will not be generated.

## 📦 **Installation from GitHub**

You can install the development version from GitHub. Within R, install `devtools`:

```r
install.packages("devtools")
```

If you are using Windows, please install the latest recommended version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

To install SIMpoly from GitHub, use:

```r
devtools::install_github("mmollina/SIMpoly", dependencies=TRUE)
```

---

## 🛠️ **Step 1: Installing and Loading Packages**

```r
library(SIMpoly)
library(mappoly2)  # for consensus mapping
```

---

## 🔬 **Step 2: Define Simulation Parameters and Pedigree**

```r
# Define ploidy for each parent.
ploidy.vec <- c(4, 2, 4, 2, 4, 4)  # for six parents
names(ploidy.vec) <- c("P1", "P2", "P3", "P4", "P5", "P6")

# Define parental crosses.
parents.mat <- matrix(c("P1", "P2",
                        "P1", "P3",
                        "P2", "P2",
                        "P3", "P4",
                        "P4", "P5",
                        "P5", "P6",
                        "P5", "P2"),
                      ncol = 2, byrow = TRUE)

# Define number of offspring per cross.
n_offspring <- c(200, 100, 200, 100, 300, 200, 100)

# Define allele set for each founder.
alleles <- list(P1 = 0:1, P2 = 0:1, P3 = 0:1, P4 = 0:1, P5 = 0:1, P6 = 0:1)

# Generate pedigree.
pedigree <- simulate_pedigree(ploidy.vec, parents.mat, n_offspring)
head(pedigree, n = 10)
```

---

## 📊 **Step 3: Simulate Multiparental Data**

```r
n.chr <- 3  # Number of chromosomes
map.len <- c(100, 120, 90)  # Chromosome lengths
n_mrk <- c(2000, 1500, 2000)  # Number of markers per chromosome

# Simulate genetic data
result <- simulate_multiparental_data(n.chr, map.len, pedigree, ploidy.vec, n_mrk, alleles, missing = 0.1, p = .3, rho = .7)

# Extract results
wide_df <- result$wide_df
all_dat <- result$dat
parent_homologs <- result$parent_homologs[[1]]
plot_parent_phase_gg(parent_homologs)
```

---

## 🔗 **Step 4: Construct Linkage Maps with mappoly2**

```r
MAPs <- vector("list", nrow(parents.mat))

for(i in 1:nrow(parents.mat)){
  p1 <- parents.mat[i,1]
  p2 <- parents.mat[i,2]
  dat <- mappoly2:::table_to_mappoly(all_dat[[i]], ploidy.vec[p1], ploidy.vec[p2], p1, p2)
  
  # Filtering steps
  dat <- filter_markers_by_chisq_pval(dat, plot = FALSE)
  dat <- filter_markers_by_missing_rate(dat, mrk.thresh = .1, plot = FALSE)
  dat <- filter_individuals_by_missing_rate(dat, ind.thresh = .1, plot = FALSE)
  dat <- pairwise_rf(dat, mrk.scope = "all", ncpus = 8)
  dat <- rf_filter(dat, probs = c(0, 1), diagnostic.plot = FALSE)
  
  # Group and order markers
  g <- group(dat, expected.groups = n.chr, comp.mat = TRUE, inter = FALSE)
  s <- make_sequence(g, ch = list(1,2,3))
  s <- order_sequence(s, type = "genome")
  s <- mapping(s, type = "genome", parent = "p1p2", ncpus = n.chr)
  
  MAPs[[i]] <- s
}
plot_map(MAPs[[1]], lg = 1, type = "genome")
```

---

## 🏆 **Acknowledgment**

This package has been developed as part of the project [AFRI-Grant: A Genetics-Based Data Analysis System for Breeders in Polyploid Breeding Programs](https://app.dimensions.ai/details/grant/grant.9958948) and  [SCRI-Grant: Tools for polyploids](https://www.polyploids.org/), funded by USDA NIFA.


<div class="horizontalgap" style="width:5px">
     <a id="USDA-NIFA" href="https://portal.nifa.usda.gov/web/crisprojectpages/1027948-a-genetics-based-data-analysis-system-for-breeders-in-polyploid-breeding-programs.html"><img src="nifa-color-lockup.png" width="650" alt=""/></a> 
      <a id="NCSU" href="https://www.ncsu.edu/"><img src="https://brand.ncsu.edu/assets/logos/ncstate-brick-2x2-red.png" width="150" alt=""/></a>
    <span class="stretch"></span>
</div>


---
<sub>NC State University promotes equal opportunity and prohibits discrimination and harassment based upon one’s age, color, disability, gender identity, genetic information, national origin, race, religion, sex (including pregnancy), sexual orientation and veteran status.</sub>

