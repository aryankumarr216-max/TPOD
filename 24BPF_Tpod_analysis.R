#24BPF_Analysis
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("shiny", "shinyWidgets", "tidyverse", "future", "future.apply",
               "rhandsontable", "jsonlite", "data.table", "Rfast", "purrr", "dendsort", "shinyFiles", "viridis")
bioc_pkgs <- c("DESeq2", "edgeR", "BiocParallel")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE, ask = FALSE)
}


setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic")



df_24BPF <- read.delim(
  "24BPF_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_BMD_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 56
)

#============================================================================
# Function: LCRD.2
# ============================================================================
#bmc = same as BMD
#probe

LCRD.2 <- function(bmc, probe, cut = 1.66) {
  stopifnot(length(bmc) == length(probe))
  if (length(bmc) < 2) return(data.frame(Gene = NA, BMC = NA, N = NA))
  
  ord <- order(bmc)
  x <- bmc[ord]
  p <- probe[ord]
  r <- x[-1] / x[-length(x)]
  p <- p[-1]
  x <- x[-1]
  n <- length(x)
  
  g <- rep(1, n)
  for (j in 2:n) {
    if (r[j] > cut) {
      g[j] <- g[j - 1] + 1
    } else {
      g[j] <- g[j - 1]
    }
  }
  
  group_sizes <- table(g)
  largest_group <- as.numeric(names(group_sizes)[group_sizes == max(group_sizes)])[1]
  idx <- which(g == largest_group)
  
  if (length(idx) == 0) return(data.frame(Gene = NA, BMC = NA, N = NA))
  
  data.frame(Gene = p[idx[1]], BMC = x[idx[1]], N = idx[1] + 1)
}

#bmc = bmd 

bmc_24BPF <- as.list(df_24BPF$`Best.BMD`)
bmc_v_24BPF <- unlist(bmc_24BPF)
probe_24BPF <-as.list(df_24BPF$`Probe.ID`)
probe_v_24BPF <- unlist(probe_24BPF)

LCRD2_24BPF_Tpod_table <- LCRD.2(bmc_v_24BPF,probe_v_24BPF)
LCRD2_24BPF_Tpod <- as.numeric(LCRD.2(bmc_v_24BPF,probe_v_24BPF)[2])

##################################################################
#FORMULA FOR FIRST MODE-Test run

First.Mode <- function(bmc, log_base = 10, bw = "nrd0", min.size = 0.1) {
  # Clean and filter input data
  valid_idx <- which(bmc > 0 & is.finite(bmc))
  if (length(valid_idx) < 2) {
    return(data.frame(First_Mode = NA, Peak_Size = NA, Total_Modes = 0))
  }
  x <- bmc[valid_idx]
  
  # Step 1: Log Transformation
  log_x <- log(x, base = log_base)
  
  # Step 2: Density Estimation (Smoothing)
  dens <- density(log_x, bw = bw)
  sumdens <- sum(dens$y)
  
  # Step 3: Map the slopes (1 for uphill, 0 for downhill)
  y.diff <- diff(dens$y)
  incr <- rep(0, length(y.diff))
  incr[which(y.diff > 0)] <- 1
  
  # Find boundaries where direction changes
  begin <- 1
  count <- 1
  for (i in 2:length(incr)) {
    if (incr[i] != incr[i - 1]) {
      count <- count + 1
      begin <- c(begin, i)
    }
  }
  begin <- c(begin, length(incr))
  
  # Step 4: Extract Modes (Peaks)
  modes <- numeric()
  sizes <- numeric()
  
  # Skip the first boundary if the curve starts by going downhill
  init <- ifelse(incr[1] == 0, 2, 1)
  
  j <- init
  while (j <= (length(begin) - 2)) {
    temp.x <- dens$x[begin[j]:begin[j + 2]]
    temp.y <- dens$y[begin[j]:begin[j + 2]]
    
    # Locate highest point in this segment
    high.point <- median(which(temp.y == max(temp.y)))
    modes <- c(modes, temp.x[high.point])
    
    # Calculate area/size of this peak
    sizes <- c(sizes, sum(dens$y[begin[j]:begin[j + 2]]) / sumdens)
    
    j <- j + 2
  }
  
  # Step 5: Filter by minimum size threshold
  valid_peaks <- which(sizes >= min.size)
  if (length(valid_peaks) == 0) {
    return(data.frame(First_Mode = NA, Peak_Size = NA, Total_Modes = 0))
  }
  
  modes <- modes[valid_peaks]
  sizes <- sizes[valid_peaks]
  
  # Step 6: Select the first mode and apply linear reversion
  first_mode_log <- modes[1]
  first_mode_linear <- log_base^first_mode_log
  
  data.frame(First_Mode = first_mode_linear, Peak_Size = sizes[1], Total_Modes = length(modes))
}

First_mode_24BPF_table<-First.Mode(bmc_v_24BPF)
First_mode_24BPF_Tpod <- First_mode_24BPF_table$First_Mode
##################################################################
#tpod from 20th gene:

nth_gene_bootstrap <- function(x, seed = 1, nth_gene = 20, repeats = 2000) {
  set.seed(seed)
  boot_nth_gene <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    sort(sampleData)[nth_gene]
  })
  return(quantile(boot_nth_gene, probs = c(0.025, 0.5, 0.975)))
}


nth_gene_bootstrap(bmc_v_24BPF, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_24BPF_table <- nth_gene_bootstrap(bmc_v_24BPF, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_24BPF_Tpod <- as.numeric(twenty_gene_24BPF_table[2])

##################################################################
#tpod from 10 percentile


nth_percent_bootstrap <- function(x, seed = 1, nth_percent = 10, repeats = 2000) {
  set.seed(seed)
  boot_nth_percent <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    quantile(sort(sampleData), probs = (nth_percent / 100))
  })
  return(quantile(boot_nth_percent, probs = c(0.025, 0.5, 0.975)))
}

nth_percent_bootstrap (bmc_v_24BPF, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_24BPF_table <- nth_percent_bootstrap (bmc_v_24BPF, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_24BPF_Tpod <- as.numeric(tenth_percentile_24BPF_table[2])
##################################################################
#Defined category Analysis Tpod (lowest median BMD)


setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Defined Category analysis")


df_24BPF_defCAT <- read.delim(
  "24BPF_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_BMD_null_DEFINED-CatMapFile- Gallus_true_true_conf0.5_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 70
)

D_24BPF_Cat_tpod <- min(df_24BPF_defCAT$BMD.Median, na.rm = TRUE)


###########################################################################
#tpod summary

D_24BPF_tpod_summary <- data.frame(LCRD2 = LCRD2_24BPF_Tpod, First_Mode = First_mode_24BPF_Tpod,twentieth_gene_lower_bound = twenty_gene_24BPF_table[1], twentieth_gene = twenty_gene_24BPF_Tpod,twentieth_gene_upper_bound = twenty_gene_24BPF_table[3], tenth_percentile_lower_bound = tenth_percentile_24BPF_table[1],tenth_percentile = tenth_percentile_24BPF_Tpod, tenth_percentile_upper_bound = tenth_percentile_24BPF_table[3], Lowest_BMC_median_Categorical = D_24BPF_Cat_tpod )

##########################################################################
##########################################################################
#plot

# Calculate the 10th percentile
p10 <- quantile(df_24BPF$Best.BMD, 0.10, na.rm = TRUE)

# Density plot with percentile marker
ggplot(df_24BPF, aes(x = Best.BMD)) +
  geom_density(fill = "lightblue", alpha = 0.4) +
  geom_vline(xintercept = p10,
             linetype = "dashed",
             color = "red") +
  annotate("text",
           x = p10,
           y = 0.02,
           label = "10th percentile",
           hjust = -0.1) +
  labs(x = "Best.BMD",
       y = "Density")



df <- df_24BPF %>%
  mutate(logBestBMD = log10(Best.BMD))

ggplot(df, aes(x = logBestBMD)) +
  geom_density() +
  labs(x = "log10(Best.BMD)")






#tpod variables: LCRD2_24BPF_Tpod, LCRD2_24BPF_Tpod, twenty_gene_24BPF_Tpod, tenth_perecentile_24BPF_Tpod, 24BPF_Cat_tpod

df <- df_24BPF %>%
  mutate(logBestBMD = log10(Best.BMD))

ggplot(df, aes(x = logBestBMD)) +
  geom_density(fill = "lightblue", alpha = 0.4) +
  
  geom_vline(xintercept = log10(p10),
             linetype = "dashed",
             color = "red") +
  
  geom_vline(xintercept = log10(LCRD2_24BPF_Tpod),
             linetype = "dashed",
             color = "blue") +
  
  geom_vline(xintercept = log10(First_mode_24BPF_Tpod),
             linetype = "dashed",
             color = "green") +
  
  geom_vline(xintercept = log10(twenty_gene_24BPF_Tpod),
             linetype = "dashed",
             color = "purple") +
  
  geom_vline(xintercept = log10(D_24BPF_Cat_tpod),
             linetype = "dashed",
             color = "pink") +
  
  labs(x = "log10(Best.BMD)",
       y = "Density")

######################################################
#Separate lines
#dev.off()

library(dplyr)
library(ggplot2)

df <- df_24BPF %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    p10,
    LCRD2_24BPF_Tpod,
    First_mode_24BPF_Tpod,
    twenty_gene_24BPF_Tpod,
    D_24BPF_Cat_tpod
  )),
  
  label = c(
    paste0("10th Percentile = ", round(p10, 3)),
    paste0("LCRD = ", round((LCRD2_24BPF_Tpod), 3)),
    paste0("First Mode = ", round((First_mode_24BPF_Tpod), 3)),
    paste0("20th Gene = ", round((twenty_gene_24BPF_Tpod), 3)),
    paste0("Category = ", round((D_24BPF_Cat_tpod), 3))
  )
)

ggplot(df, aes(x = logBestBMD)) +
  
  geom_density(fill = "lightblue", alpha = 0.4) +
  
  geom_vline(data = tpod_df,
             aes(xintercept = value,
                 color = label),
             linetype = "dashed",
             linewidth = 1) +
  
  geom_text(data = tpod_df,
            aes(x = value + 0.05,
                y = c(0.02, 0.04, 0.06, 0.08, 0.10),
                label = round(value, 3),
                color = label),
            angle = 90,
            hjust = 0,
            show.legend = FALSE) +
  
  labs(x = "log10(Best.BMD)",
       y = "Density",
       color = "TPoD Methods",
       title = "                                            Distribution of Best BMD Values for 24BPF")


