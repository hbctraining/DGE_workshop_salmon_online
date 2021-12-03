## Setting up

# 1. Let’s create a new project directory for this review:
#   
# - Create a new project called `R_refresher`

# File -> New Project -> New Directory - > New Project -> Name it "R_refresher" -> Create Project

# - Create a new R script called `reviewing_R.R`

# File -> New File -> R script. Then, File -> Save as... -> Name it "reviewing_R" -> Save

# - Create the following folders in the project directory - `data`, `figures`

# With the Files tab selected in the Files/Plots/Packages/Help window, click "New folder" and name it "data" or "figures" then click "OK".

# - Download a counts file to the `data` folder by [right-clicking here](https://github.com/hbctraining/DGE_workshop_salmon/blob/master/data/raw_counts_mouseKO.csv?raw=true)

# Right click the hyperlink and select "Download Linked File As" or "Save Link As". Navigate to your data folder within your R_refresher directory and click "Save".

# 2. Now that we have our directory structure setup, let's load our libraries and read in our data:
# 
#     - Load the `tidyverse` library
library(tidyverse)

#     - Use `read.csv()` to read in the downloaded file and save it in the object/variable `counts`
counts <- read.csv("data/raw_counts_mouseKO.csv")

#     - What is the syntax for a function?

# Object <- read.csv("path/to/file.txt")

#     - How do we get help for using a function?

# ?read.csv

#     - What is the data structure of `counts`?
class(counts)

# This will tell us that the counts object if a data frame.

#     - What main data structures are available in R?

# Data frames, vectors, lists, matrices and factors

#     - What are the data types of the columns?
str(counts)

# We can see that they are eight numeric columns.

#     - What data types are available in R?

# Numeric, character, integer, logical, complex and raw
       
# ## Creating vectors/factors and dataframes
# 
# 3. We are performing RNA-Seq on cancer samples with genotypes of p53 wildtype (WT) and knock-down (KO). You have 8 samples total, with 4 replicates per genotype. Write the R code you would use to construct your metadata table as described below.  
# 
#      - Create the vectors/factors for each column (Hint: you can type out each vector/factor, or if you want the process go faster try exploring the `rep()` function).

sex <- rep(c("M","F"), 4)
stage <- c(1,2,2,1,2,1,1,2)
genotype <- c(rep("KO", 4), rep("WT",4))
myc <- c(23,4,45,90,34,35,9, 10)

#      - Put them together into a dataframe called `meta`.

meta <- data.frame(sex, stage, genotype,myc)
View(meta)

#      - Use the `rownames()` function to assign row names to the dataframe (Hint: you can type out the row names as a vector, or if you want the process go faster try exploring the `paste0()` function).

rownames(meta) <- c(paste0(rep("KO", 4), 1:4), paste0(rep("WT",4), 1:4))

#     Your finished metadata table should have information for the variables `sex`, `stage`, `genotype`, and `myc` levels: 
# 
#     | |sex	| stage	| genotype	| myc |
#     |:--:|:--: | :--:	| :------:	| :--: |
#     |KO1 |	M	|1	|KO	|23|
#     |KO2|	F	|2	|KO	|4|
#     |KO3	|M	|2	|KO	|45|
#     |KO4	|F	|1	|KO	|90|
#     |WT1|	M	|2	|WT	|34|
#     |WT2|	F|	1|	WT|	35|
#     |WT3|	M|	1|	WT|	9|
#     |WT4|	F|	2|	WT|	10|



#### Exploring data
# 
# Now that we have created our metadata data frame, it's often a good idea to get some descriptive statistics about the data before performing any analyses. 
# 
#      - Summarize the contents of the `meta` object, how many data types are represented?
str(meta)

# We can see that two columns are numeric and two are characters.

#   - Check that the row names in the `meta` data frame are identical to the column names in `counts` (content and order).
all(rownames(meta) %in% colnames(counts))
all(rownames(meta) == colnames(counts))

# Either of these should return TRUE indicating that we have passed this check.

#      - Convert the existing `stage` column into a factor data type
meta$stage <- factor(meta$stage)

str(meta)

# This should show you that the "Stage" column is now a factor.

# 
# ## Extracting data
# 
# 4. Using the `meta` data frame created in the previous question, perform the following exercises (questions **DO NOT** build upon each other):
#   
#      - return only the `genotype` and `sex` columns using `[]`:
meta[,c(3,1)]
#Or
meta[,c("genotype","sex")]

#      - return the `genotype` values for samples 1, 7, and 8 using `[]`:
meta[c(1,7,8),3]
# Or
meta[c(1,7,8),"genotype"]

#      - use `filter()` to return all data for those samples with genotype `WT`:
filter(meta, genotype == "WT")
# OR
meta %>% filter(genotype == "WT")

#      - use `filter()`/`select()`to return only the `stage` and `genotype` columns for those samples with `myc` > 50:
meta %>% 
  filter(myc > 50) %>% 
  select(stage, genotype)

# Or 

select(filter(meta, myc > 50), stage, genotype)

#      - add a column called `pre_treatment` to the beginning of the dataframe with the values T, F, T, F, T, F, T, F 
pre_treatment <- c(T, F, T, F, T, F, T, F)

# Or

pre_treatment <- rep(c("T","F"),4)

meta <- cbind(pre_treatment, meta)

#      - why might this design be problematic?

# All of the Male samples are True for the pre-treatment and all of the females are False. Thus, we have confounded our design.

#      - Using `%>%` create a tibble of the `meta` object and call it `meta_tb` (make sure you don't lose the rownames!)
#      - change the names of the columns to: "A", "B", "C", "D", "E":
#      
meta_tb <- meta %>% 
              rownames_to_column(var="sampleIDs") %>%
              as.tibble()

colnames(meta_tb)[2:6] <- LETTERS[1:5]

# ## Visualizing data
# 
# 5. Often it is easier to see the patterns or nature of our data when we explore it visually with a variety of graphics. Let's use ggplot2 to explore differences in the expression of the Myc gene based on genotype.
#                                                                             
#      - Plot a boxplot of the expression of Myc for the KO and WT samples using `theme_minimal()` and give the plot new axes names and a centered title.
#                   
ggplot(meta) +
  geom_boxplot(aes(x = genotype, y = myc)) +
  theme_minimal() +
  ggtitle("Myc expression") +
  ylab("Myc level") +
  xlab("Genotype") +
  theme(plot.title = element_text(hjust=0.5, size = rel(2)))

### Preparing for downstream analysis tools
#                                                                             
# 6. Many different statistical tools or analytical packages expect all data needed as input to be in the structure of a list. Let's create a list of our count and metadata in preparation for a downstream analysis.
# 
#      - Create a list called `project1` with the `meta` and `counts` objects, as well as a new vector with all the sample names extracted from one of the 2 data frames.
#  
project1 <- list(meta, counts, rownames(meta))
project1
