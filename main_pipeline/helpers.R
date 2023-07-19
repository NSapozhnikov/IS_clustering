## Functions
suppressPackageStartupMessages(require(dplyr))

GetBlocks <- function(x) {
  # Get one best block for each cluster from omnibus test results.
  # Best block has minimum p-value and maximum size.
  # Adjust p-value for multiple testing. 
  # Args:
  #  x: a dataframe with the results of Omnibus test
  
  # For each cluster get blocks with minimum p-value and maximum size
  d1 <- x %>% group_by(CHR, START) %>%
    slice_min(P) %>%
    slice_max(SIZE) %>%
    ungroup()
 
  # Add adjusted p-values
  d1 %>% mutate(BONF = p.adjust(P, method = "bonferroni"),
                                HOLM = p.adjust(P, method = "holm"),
                                BH = p.adjust(P, method = "BH")) 
}

GetSNP <- function(x) {
  stringr::str_split(string = x, pattern = "\\|") %>% unlist()
}

GetHaps <- function(x) {
  ## Find haplotypes with minimum p-value, define the blocks, 
  ## subset haplotypes from these blocks and adjust p-value
  ## Args: 
  ##   x: a list with dataframes having the results of --hap-assoc Plink tests
  # Get distinct blocks with minimum p-value and maximum size in each cluster
  d1 <- x %>% group_by(CHR, START) %>% 
    slice_min(P) %>%
    slice_max(SIZE) %>%
    select(CHR, START, END) %>% 
    distinct() %>%
    ungroup()
  # Subset haplotypes in distinct blocks
  d2 <- left_join(d1, x, by = c("CHR", "START", "END"))
  # Add adjusted p-value
  d3 <- d2 %>% mutate(BONF = p.adjust(P, method = "bonferroni"),
                      HOLM = p.adjust(P, method = "holm"),
                      BH = p.adjust(P, method = "BH"))
}


GetSize <- function(x) {
  # Split the string by "|" and return the number of elements
  sapply(x, function(s) stringr::str_split(
    string = s, 
    pattern = "\\|") %>% unlist() %>% length(), USE.NAMES = F) 
}

LoadOmnibus <- function(l) {
  ## Subset the records having OMNIBUS from files created with --hap-assoc 
  ## argument of Plink 1.07. 
  ## Args:
  ##   l: a list of paths to csv files from different experiments
  plyr::llply(l, function(x) {
    plyr::ldply(x, function(y) {
      # Read header
      cn <- colnames(data.table::fread(y, nrows = 0))
      dt1 <- data.table::fread(cmd = paste("grep", "OMNIBUS", y))
      colnames(dt1) <- cn
      dt1 %>% 
        select(LOCUS, CHISQ, DF, P, SNPS) %>% 
        mutate(CHR = stringr::str_match(LOCUS, "([0-9]+):")[, 2] %>% as.integer(),
               START = stringr::str_match(LOCUS, ":([0-9]+)-")[, 2] %>% as.integer(),
               END = stringr::str_match(LOCUS, "-([0-9]+)")[, 2] %>% as.integer(),
               LENGTH = END - START,
               SIZE = GetSize(SNPS))
    }, .progress = "text")
  })
  
}

LoadNonOmnibus <- function(l) {
  ## Load the files created with --hap-assoc argument of Plink 1.07
  ## Args:
  ##   l: a list of paths to csv files from different experiments
  plyr::llply(l, function(x) {
    plyr::ldply(x, function(y) {
      data.table::fread(y) %>% 
        mutate(CHR = stringr::str_match(LOCUS, "([0-9]+):")[, 2] %>% as.integer(),
               START = stringr::str_match(LOCUS, ":([0-9]+)-")[, 2] %>% as.integer(),
               END = stringr::str_match(LOCUS, "-([0-9]+)")[, 2] %>% as.integer(),
               LENGTH = END - START,
               SIZE = GetSize(SNPS))
    }, .progress = "text")
  })
  
}



PlotLD <- function(r, gd, ld.folder) {
  
  snps1 <- r$SNPS %>% stringr::str_split(pattern = "\\|") %>% unlist()
  
  # Load map
  map <- read.table(sprintf("%s/chr%s.bim", gd, r$CHR))
  
  # Subset list of all SNPs in the region
  lim <- c(which(map$V2 == snps1[1]),
           which(map$V2 == snps1[length(snps1)]))
  
  snps2 <- map$V2[lim[1]:lim[2]]
  
  # Save
  snp.pref <- paste0(ld.folder, "/chr", paste0(r$CHR, "-", r$START, "-", r$END))
  snp.list <- paste0(snp.pref, ".txt") 
  write(snps2, snp.list)
  
  # Count LD
  system2("/home/gennady/tools/plink_5.2/plink", 
          c("--bfile", sprintf("%s/chr%s", gd, r$CHR), "--r2 square", "--make-founders", 
            "--extract", snp.list, "--out", snp.pref))
  
  # Load LD
  ld1 <- read.table(paste0(snp.pref, ".ld"))
  colnames(ld1) <- snps2
  rownames(ld1) <- snps2
  
  
  # Set pdf 
  pdf(paste0(snp.pref, ".pdf"), width = 11.7, height = 8.3)
  # Set graphical parameters
  old <- par(mfrow = c(1, 2))
  
  # Plot LD 
  snp.positions1 = map$V4[lim[1]:lim[2]]
  gaston::LD.plot(LD = ld1, snp.positions = snp.positions1, 
                  write.snp.id = T, draw.chr = T)
  
  # Leave only SNPs from the cluster
  ld2 <- ld1[rownames(ld1) %in% snps1, colnames(ld1) %in% snps1]
  snp.positions2 <- snp.positions1[rownames(ld1) %in% snps1]
  
  # Plot LD
  gaston::LD.plot(LD = ld2, snp.positions = snp.positions2, 
                  write.snp.id = T, draw.chr = T)
  invisible(dev.off())
  par(old)
  
}

ReduceHaps <- function(d4) {
  # Drop LOCUS column
  d4$LOCUS <- NULL
  # Create combinations of indexes
  n <- nrow(d4)
  ind1 <- combn(n, 2, simplify = F)
  # Iterate over the indexes and setup the rows to be removed
  ind2 <- sapply(ind1, function(x) {
    h1 <- d4 %>% slice(x[1]) %>% unlist()
    h2 <- d4 %>% slice(x[2]) %>% unlist()
    if(h1["CHR"] == h2["CHR"] & h1["START"] == h2["START"]) {
      ifelse(h1["SIZE"] < h2["SIZE"], return(x[1]), return(x[2]))
    } 
  }) %>% unlist() %>% unique()
  # Remove haplotypes of less size and order by genomic coordinates
  d5 <- d4 %>% slice(-ind2) %>% arrange(CHR, START)
  d5
}

GetNumClusters <- function(x) {
  # Group by CHR, START 
  dt <- x %>% count(CHR, START)
  # Get the number of clusters
  n <- nrow(dt)
  n
}

SubsetHaplotypesByPvalue <- function(x) {
  # Group by CHR, START 
  d1 <- x %>% count(CHR, START)
  # Get the number of clusters
  n <- nrow(d1)
  # Get threshold value of p-value
  p <- 0.05/n
  # Subset haplotypes with p-value less then threshold one
  x %>% filter(P < p)
}

Consolidate <- function(files) {
  ## Consolidate *sign.haps.csv files
  ## Args:
  ##  files: a vector with paths to *.csv files created with Plink
  ## Value:
  ##  A data frame with all haplotypes.
  # Load data
  d1 <- lapply(files, read.csv)
  d2 <- do.call("rbind", d1)
  # Extract the genomic coordinates
  d3 <- data.frame(
    CHR = stringr::str_match(d2$LOCUS, "([0-9]+):")[, 2] %>% as.integer(),
    START = stringr::str_match(d2$LOCUS, ":([0-9]+)-")[, 2] %>% as.integer(),
    END = stringr::str_match(d2$LOCUS, "-([0-9]+)")[, 2] %>% as.integer())
  # Create data frame with genomic coordinates
  d2$LOCUS <- NULL
  d4 <- cbind(d3, d2)
  # Create combinations of indexes
  n <- nrow(d4)
  ind1 <- combn(n, 2, simplify = F)
  # Iterate over the indexes and setup the rows to be removed
  ind2 <- sapply(ind1, function(x) {
    h1 <- d4 %>% slice(x[1]) %>% unlist()
    h2 <- d4 %>% slice(x[2]) %>% unlist()
    if(h1["CHR"] == h2["CHR"] & h1["START"] == h2["START"]) {
      ifelse(h1["SIZE"] < h2["SIZE"], return(x[1]), return(x[2]))
    } 
  }) %>% unlist() %>% unique()
  # Remove haplotypes of less size and order by genomic coordinates
  d5 <- d4 %>% slice(-ind2) %>% arrange(CHR, START)
  d5
}

