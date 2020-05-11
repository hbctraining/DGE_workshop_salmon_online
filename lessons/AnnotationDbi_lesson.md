## AnnotationDbi

AnnotationDbi is an R package that provides an interface for connecting and querying various annotation databases using SQLite data storage. The AnnotationDbi packages can query the *OrgDb*, *TxDb*, *EnsDb*, *Go.db*, and *BioMart* annotations. There is helpful [documentation](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf) available to reference when extracting data from any of these databases.

### org.Hs.eg.db

There are a plethora of organism-specific *orgDb* packages, such as `org.Hs.eg.db` for human and `org.Mm.eg.db` for mouse, and a list of organism databases can be found [here](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb). These databases are best for converting gene IDs or obtaining GO information for current genome builds, but not for older genome builds. These packages provide the current builds corresponding to the release date of the package, and update every 6 months. If a package is not available for your organism of interest, you can create your own using *AnnotationHub*.

```r
# Load libraries
library(org.Hs.eg.db)
library(AnnotationDbi)

# Check object metadata
org.Hs.eg.db
```

We can see the metadata for the database by just typing the name of the database, including the species, last updates for the different source information, and the source urls. Note the KEGG data from this database was last updated in 2011, so may not be the best site for KEGG pathway information.

```r
OrgDb object:
| DBSCHEMAVERSION: 2.1
| Db type: OrgDb
| Supporting package: AnnotationDbi
| DBSCHEMA: HUMAN_DB
| ORGANISM: Homo sapiens
| SPECIES: Human
| EGSOURCEDATE: 2018-Oct11
| EGSOURCENAME: Entrez Gene
| EGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
| CENTRALID: EG
| TAXID: 9606
| GOSOURCENAME: Gene Ontology
| GOSOURCEURL: ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/
| GOSOURCEDATE: 2018-Oct10
| GOEGSOURCEDATE: 2018-Oct11
| GOEGSOURCENAME: Entrez Gene
| GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
| KEGGSOURCENAME: KEGG GENOME
| KEGGSOURCEURL: ftp://ftp.genome.jp/pub/kegg/genomes
| KEGGSOURCEDATE: 2011-Mar15
| GPSOURCENAME: UCSC Genome Bioinformatics (Homo sapiens)
| GPSOURCEURL: 
| GPSOURCEDATE: 2018-Oct2
| ENSOURCEDATE: 2018-Oct05
| ENSOURCENAME: Ensembl
| ENSOURCEURL: ftp://ftp.ensembl.org/pub/current_fasta
| UPSOURCENAME: Uniprot
| UPSOURCEURL: http://www.UniProt.org/
| UPSOURCEDATE: Thu Oct 18 05:22:10 2018
```

We can easily extract information from this database using *AnnotationDbi* with the methods: `columns`, `keys`, `keytypes`, and `select`. For example, we will use our `org.Hs.eg.db` database to acquire information, but know that the same methods work for the *TxDb*, *Go.db*, *EnsDb*, and *BioMart* annotations.

```r
# Return the Ensembl IDs for a set of genes
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                     keys = res_tableOE_tb$gene,  # data to use for retrieval
                                     columns = c("SYMBOL", "ENTREZID","GENENAME"), # information to retreive for given data
                                     keytype = "ENSEMBL") # type of data given in 'keys' argument
```

We started from at about 57K genes in our results table, and the dimensions of our resulting annotation data frame also look quite similar. Let's take a peek to see if we actually returned annotations for each individual Ensembl gene ID that went in to the query:

```r
length(which(is.na(annotations_orgDb$SYMBOL)))
```

Looks like more than half of the input genes did not return any annotations. This is because the OrgDb family of database are primarily based on mapping using Entrez Gene identifiers. If you look at some of the Ensembl IDs from our query that returned NA, these map to pseudogenes (i.e [ENSG00000265439](https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000265439;r=6:44209766-44210063;t=ENST00000580735)) or non-coding RNAs (i.e. [ENSG00000265425](http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000265425;r=18:68427030-68436918;t=ENST00000577835)). The difference is due to the fact that each database implements different computational approaches for generating the gene builds. Let's get rid of those NA entries:

```r
# Determine the indices for the non-NA genes
non_na_idx <- which(is.na(annotations_orgDb$SYMBOL) == FALSE)

# Return only the genes with annotations using indices
annotations_orgDb <- annotations_orgDb[non_na_idx, ]
```

You may have also noted the *warning* returned: *'select()' returned 1:many mapping between keys and columns*. This is always going to happen with converting between different gene IDs (i.e. one geneID can map to more than one identifier in another databse) . Unless we would like to keep multiple mappings for a single gene, then we probably want to de-duplicate our data before using it.

```r
# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]
```

### EnsDb.Hsapiens.v86

To generate the Ensembl annotations, the *EnsDb* database can also be easily queried using AnnotationDbi. You will need to decide the release of Ensembl you would like to query. We know that our data is for GRCh38, and the most current *EnsDb* release for GRCh38 in Bioconductor is release 86, so we can install this database. All Ensembl releases are listed [here](http://useast.ensembl.org/info/website/archives/index.html). **NOTE: this is not the most current release of GRCh38 in the Ensembl database, but it's as current as we can obtain through AnnotationDbi.**

Since we are using *AnnotationDbi* to query the database, we can use the same functions that we used previously:

```r
# Load the library
library(EnsDb.Hsapiens.v86)

# Check object metadata
EnsDb.Hsapiens.v86

# Explore the fields that can be used as keys
keytypes(EnsDb.Hsapiens.v86)
```

Now we can return all gene IDs for our gene list:

```r
# Return the Ensembl IDs for a set of genes
annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                           keys = res_tableOE_tb$gene,
                                           columns = c("SYMBOL", "ENTREZID","GENEBIOTYPE"),
                                           keytype = "GENEID")
```

We can check for NA entries, and find that there are none:

```r
length(which(is.na(annotations_edb$SYMBOL) == FALSE))
```

Then we can again deduplicate, to remove the gene symbols which appear more than once:

```r
# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_edb$SYMBOL) == FALSE)

# Return only the non-duplicated genes using indices
annotations_edb <- annotations_edb[non_duplicates_idx, ]
```

> **NOTE:** In this case we used the same build but a slightly older release, and we found little discrepancy. If your analysis was conducted using an older genome build (i.e hg19), but used a newer build for annotation some genes may be found to be not annotated (NA). Some of the genes have changed names in between versions (due to updates and patches), so may not be present in the newer version of the database. 

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
