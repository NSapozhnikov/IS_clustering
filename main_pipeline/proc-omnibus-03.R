## Process the results of omnibus tests for dataset A
## The parameters of clustering were optimized
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))
source("helpers.R")

# Initiate 
folder <- c(dbscan = "/mnt/wd/nsap/Clustering/main_pipeline/hap_assoc/csvs/dbscan",
            hdbscan = "/mnt/wd/nsap/Clustering/main_pipeline/hap_assoc/csvs/hdbscan")

# Get the files 
f1 <- lapply(folder, function(x) {
  list.files(x, pattern = "^A_chr", full.names = T) %>% 
  subset(x = (.), subset = grepl(pattern = "omnibus", x = (.)))
}) 
  
# Check file exists
sapply(c("dbscan", "hdbscan"), function(x) {
  sapply(f1[[x]], function(y) {
    if(!file.exists(y)) stop(y, " doesn't exist", call. = F)
  })
})

# Load data
l1 <- LoadOmnibus(f1["hdbscan"])

# Find best blocks in each cluster
b1 <- plyr::llply(l1, GetBlocks)

# Plot histogram of p-values
ggplot(b1[["hdbscan"]], aes(x = P)) + 
  geom_histogram(binwidth = 0.05, fill = "lightgrey", color = "darkgrey") +
  labs(x = "p-value", y = "Counts")


# Save best blocks
data.table::fwrite(b1[["hdbscan"]], "out/a-hap-assoc-omnibus-hdbscan.csv")


