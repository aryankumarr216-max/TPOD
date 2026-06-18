####################################################################
#BPZ_Analysis
#Load in Libraries
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


#Set Working Directory
setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic")



df_BPZ <- read.delim(
  "BPZ_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_BMD_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 56
)

#TDR = top dose removed

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/top_dose_removed_BDA")

df_BPZ_TDR <- read.delim(
  "BPZ_bmdexpress_input_log2_transformed_top_dose_removed_williams_0.05_NOMTC_foldfilter1.5_BMD_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 59
)
############### DRG Calculation #############################
DRG_BPZ <-as.numeric(nrow(df_BPZ))

DRG_BPZ_TDR <-as.numeric(nrow(df_BPZ_TDR))


############### Top Dose and Lose Dose Calculation ########################

setwd("C:/Users/KumarA/Downloads/Gene_count_matrices/Gene_count_matrices")

dose_BPZ <- read.delim(
  "BPZ_bmdexpress_input_log2_transformed.txt",
  sep = "\t",
  header = TRUE)

low_dose_BPZ <- min(dose_BPZ[1, ][dose_BPZ[1, ] != 0], na.rm = TRUE)
high_dose_BPZ <-max(as.numeric(dose_BPZ[1, ]), na.rm = TRUE)

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic")


#TDR
setwd("C:/Users/KumarA/Downloads/Gene_count_matrices/Gene_count_matrices/Top_Dose_removed")

dose_D_BPZ_TDR <- read.delim(
  "BPZ_bmdexpress_input_log2_transformed_top_dose_removed.txt",
  sep = "\t",
  header = TRUE)

low_dose_BPZ_TDR <- min(dose_D_BPZ_TDR[1, ][dose_D_BPZ_TDR[1, ] != 0], na.rm = TRUE)
high_dose_BPZ_TDR <-max(as.numeric(dose_D_BPZ_TDR[1, ]), na.rm = TRUE)

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic")

####################################################################
#Boot strap function (summarizes results )
bootstrap_results <-
  function(logBMDvalues, func, endpoint_name, ...) {
    results <- as.data.frame(t(func(logBMDvalues, ...))) %>%
      dplyr::mutate(endpoint = endpoint_name) %>%
      dplyr::mutate(chemical = chemname) %>%
      relocate(chemical, endpoint)
    
    return(results)
  }



# ============================================================================
# Function: LCRD (non-bootstrap)
# ============================================================================

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

bmc_BPZ <- as.list(df_BPZ$`Best.BMD`)
bmc_v_BPZ <- unlist(bmc_BPZ)
probe_BPZ <-as.list(df_BPZ$`Probe.ID`)
probe_v_BPZ <- unlist(probe_BPZ)


#TDR
bmc_BPZ_TDR <- as.list(df_BPZ_TDR$`Best.BMD`)
bmc_v_BPZ_TDR <- unlist(bmc_BPZ_TDR)
probe_BPZ_TDR <-as.list(df_BPZ_TDR$`Probe.ID`)
probe_v_BPZ_TDR <- unlist(probe_BPZ_TDR)


#LCRD2_BPZ_Tpod_table_non_boot_strap <-LCRD.2(bmc_v_BPZ,probe_v_BPZ)
#LCRD2_BPZ_Tpod_non_boot_strap <- LCRD2_BPZ_Tpod_table$BMC

# ============================================================================
# Function: LCRD (bootstrap)
# ============================================================================

LCRD <- function(bmc, probe, cut = 1.778, logbase = 10) {
  # Order BMC values and their corresponding probes
  ordered_indices <- order(bmc)
  sorted_bmc <- bmc[ordered_indices]
  sorted_probe <- probe[ordered_indices]
  
  # Calculate ratios
  ratios <- sorted_bmc[-1] / sorted_bmc[-length(sorted_bmc)]
  
  # Initialize variables to store CRGB information
  crgb_list <- list()
  current_start <- 1
  
  # Iterate through the ratios to find CRGBs
  for (i in seq_along(ratios)) {
    if (ratios[i] >= cut) {
      # Append each CRGB together into a list
      crgb_list <- append(crgb_list, list(list(start = current_start, end = i)))
      current_start <- i + 1
    }
  }
  # Append the last CRGB
  crgb_list <- append(crgb_list, list(list(start = current_start, end = length(sorted_bmc))))
  
  # Identify the largest CRGB
  crgb_lengths <- sapply(crgb_list, function(group) group$end - group$start + 1)
  largest_crgb_index <- which.max(crgb_lengths)
  largest_crgb <- crgb_list[[largest_crgb_index]]
  
  # Extract BMCs and probes of the largest CRGB
  largest_crgb_bmcs <- sorted_bmc[largest_crgb$start:largest_crgb$end]
  largest_crgb_probes <- sorted_probe[largest_crgb$start:largest_crgb$end]
  
  # Identify the lowest BMC in the largest CRGB
  lcrd_bmc <- min(largest_crgb_bmcs)
  lcrd_probe <- largest_crgb_probes[which.min(largest_crgb_bmcs)]
  
  # Extract the lowest BMC for each CRGB
  crgb_bmcs_df <- do.call(rbind, lapply(seq_along(crgb_list), function(i) {
    group <- crgb_list[[i]]
    crgb_bmcs <- sorted_bmc[group$start:group$end]
    crgb_probes <- sorted_probe[group$start:group$end]
    data.frame(
      Gene = crgb_probes[which.min(crgb_bmcs)],
      BMC = min(crgb_bmcs),
      logBMC = log(min(crgb_bmcs), base = logbase),
      CRGB_Number = i,
      CRGB_Size = length(crgb_bmcs)
    )
  }))
  
  # Return the result
  list(
    LCRD_Result = data.frame(
      Gene = lcrd_probe,
      BMC = lcrd_bmc,
      logBMC = log(lcrd_bmc, base = logbase),
      CRGB_Number = largest_crgb_index,
      CRGB_Size = length(largest_crgb_bmcs)
    ),
    CRGB_BMCs = crgb_bmcs_df
  )
}





lcrd_bootstrap <- function(x, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10) {
  set.seed(seed)
  bmds <- x[["Best.BMD"]]
  names(bmds) <- x[["Probe.ID"]]
  if (all(names(bmds) == x[["Probe.ID"]], na.rm = TRUE) &&
      all(bmds == x[["Best.BMD"]])) {
    boot_lcrd <- replicate(repeats, {
      sampleData <- sample(unlist(bmds), length(unlist(bmds)), replace = TRUE)
      result <-
        LCRD(bmc = sampleData,
             probe = names(sampleData),
             cut = lcrdratiocut,
             logbase = lcrdlogbase)
      result$LCRD_Result$logBMC
    })
    return(quantile(boot_lcrd, probs = c(0.025, 0.5, 0.975)))
  } else {
    print("Probe IDs do not match BMDs, please review Bootstrapping_Funtions.R")
  }
}

lcrd_bootstrap(df_BPZ, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10)
library(dplyr)
chemname <- "BPZ"

LCRD_BPZ_Tpod_table<- lcrd_bootstrap(df_BPZ, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10)
LCRD2_BPZ_Tpod_median <-10^(as.numeric(LCRD_BPZ_Tpod_table[2]))
LCRD2_BPZ_Tpod_lower<-10^(as.numeric(LCRD_BPZ_Tpod_table[1]))
LCRD2_BPZ_Tpod_upper<-10^(as.numeric(LCRD_BPZ_Tpod_table[3]))
  

#TDR LCRD
LCRD_BPZ_TDR_Tpod_table<- lcrd_bootstrap(df_BPZ_TDR, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10)
LCRD2_BPZ_TDR_Tpod_median <-10^(as.numeric(LCRD_BPZ_TDR_Tpod_table[2]))
LCRD2_BPZ_TDR_Tpod_lower<-10^(as.numeric(LCRD_BPZ_TDR_Tpod_table[1]))
LCRD2_BPZ_TDR_Tpod_upper<-10^(as.numeric(LCRD_BPZ_TDR_Tpod_table[3]))


# ============================================================================
# Function: First mode (non-bootstrap)
# ============================================================================

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

First.Mode(bmc_v_BPZ)
First_mode_BPZ_table_non_boot_strap <- First.Mode(bmc_v_BPZ)
First_mode_BPZ_Tpod_non_boot_strap <- First_mode_BPZ_table$First_Mode
# ============================================================================
# Function: First mode (bootstrap)
# ============================================================================

mode.antimode<-function (x, min.size = 0.1, bw="nrd0", min.bw=NULL) 
{
  #checks
  if (missing(x)) 
    stop("The x argument is required.")
  x <- as.vector(as.numeric(as.character(x)))
  x <- x[is.finite(x)]
  if (length(unique(x))<=1) #checks if x is a constant
    return(list(modes = NA, mode.dens = NA, size = 1))
  
  #density function
  dens<-density(x,bw=bw)
  if(!is.null(min.bw)){
    if(dens$bw < min.bw){
      dens<-density(x,bw=min.bw)
    }
  }
  
  #save bw
  mode.bw<-dens$bw
  
  #the difference between each sequential y value
  y.diff<-diff(dens$y)      
  
  #points along the desntiy distribution where y is increasing
  incr<-rep(0,length(y.diff))
  incr[which(y.diff>0)]<-1
  
  #identify points where y changes direction
  begin <- 1
  count <- 1
  for (i in 2:length(incr)) {
    if (incr[i] != incr[i - 1]) {
      count <- count + 1
      begin <- c(begin, i)
    }
  }
  begin <- c(begin, length(incr))
  
  #placeholders for results
  size <- modes <- mode.dens <- rep(0, count/2)
  anti.modes<-rep(0,length(size-1))
  
  #sum of all y values
  sumdens <- sum(dens$y)
  
  ###identify modes
  init <- 1
  
  #if density plot begins on a local max (i.e. a mode) 
  if (incr[1] == 0) {
    size[1] <- sum(dens$y[1:begin[2]])/sumdens
    init <- 2
  }
  
  #identfy all modes, size and max density
  j <- init
  for (i in init:length(size)) {
    size[i] <- sum(dens$y[begin[j]:begin[j + 2]])/sumdens
    temp.x<- dens$x[begin[j]:begin[j + 2]]
    temp.y<- dens$y[begin[j]:begin[j + 2]]
    highs<-which(temp.y==max(temp.y))
    if(median(highs)%%1>0){
      high.point<-median(highs[-length(highs)])
    }else{
      high.point<-median(highs)
    }
    modes[i] <- temp.x[high.point]
    mode.dens[i] <- temp.y[high.point]
    j <- j + 2
  }
  
  #identify and remove modes smaller than min.size
  if (any(size < min.size)) {
    modes <- modes[-which(size < min.size)]
    mode.dens <- mode.dens[-which(size < min.size)]
    size <- size[-which(size < min.size)]
  }
  
  #identify anti-mode between modes
  if(length(modes)>1){
    anti.modes<-rep(0,length(size)-1)
    for(i in 1:length(anti.modes)){
      m1<-which(dens$y==mode.dens[i])
      m1<-m1[m1%in%begin]
      m2<-which(dens$y==mode.dens[i+1])
      m2<-m2[m2%in%begin]
      temp.x <- dens$x[m1:m2]
      temp.y <- dens$y[m1:m2]
      lows<-which(temp.y==min(temp.y))
      if(median(lows)%%1>0){
        low.point<-median(lows[-length(lows)])
      }else{
        low.point<-median(lows)
      }
      anti.modes[[i]]<-temp.x[low.point]
    }
  }else{
    anti.modes<-NULL
  }
  
  return(list(modes = modes, mode.dens = mode.dens, size = size, anti.modes=anti.modes, bw=mode.bw))
}

mode_bootstrap <- function(x, seed = 1, repeats = 2000, min_size = 0.06, min_bw = 0.15) {
  #source("mode_antimode.R")
  set.seed(seed)
  boot_mode <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    dataMode <- mode.antimode(sampleData, min.size = min_size, bw = "SJ", min.bw = min_bw)
    dataMode$modes[[1]]
  })
  return(quantile(boot_mode, probs = c(0.025, 0.5, 0.975)))
}

mode_bootstrap <- function(x, seed = 1, repeats = 2000, min_size = 0.06, min_bw = 0.01) {
  #source("mode_antimode.R")
  set.seed(seed)
  boot_mode <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    dataMode <- mode.antimode(sampleData, min.size = min_size, bw = "SJ", min.bw = min_bw)
    dataMode$modes[[1]]
  })
  return(quantile(boot_mode, probs = c(0.025, 0.5, 0.975)))
}


mode_bootstrap <- function(x, seed = 1, repeats = 2000) {
  set.seed(seed)
  
  #A) Pre-calc the log10 values exactly like the non-bootstrap does
  valid_idx <- which(x > 0 & is.finite(x))
  log_x <- log10(unlist(x)[valid_idx])
  
  boot_mode <- replicate(repeats, {
    # Sample the log10'd data
    sampleData <- sample(log_x, length(log_x), replace = TRUE)
    
    # Run mode.antimode on the log10'd data, matching First.Mode parameters so that the value actually computes properly
    # bw = "nrd0" or "0.01" (nrd0 will more likely have an effect if any), min.size = 0.1, and turn off min.bw
    dataMode <- mode.antimode(sampleData, min.size = 0.1, bw = "nrd0", min.bw = NULL)
    
    # Return the first mode after that:
    dataMode$modes[[1]]
  })
  
  # B) Convert the results back out of log space (linearize/exponentiate/ 10^...)
  # just like the end of First.Mode does (first_mode_linear <- log_base^first_mode_log)
  boot_mode_linear <- 10^(boot_mode)
  
  return(quantile(boot_mode_linear, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}

#mode_bootstrap(bmc_v_BPZ, seed = 1, repeats = 2000, min_size = 0.06, min_bw = 0.15)
#First_mode_BPZ_table <-mode_bootstrap(bmc_v_BPZ, seed = 1, repeats = 2000, min_size = 0.06, min_bw = 0.15)
First_mode_BPZ_table <-mode_bootstrap(bmc_v_BPZ, seed = 1, repeats = 2000)

First_mode_BPZ_Tpod_median <- as.numeric(First_mode_BPZ_table[2])
First_mode_BPZ_Tpod_lower <-  as.numeric(First_mode_BPZ_table[1])
First_mode_BPZ_Tpod_upper <-  as.numeric(First_mode_BPZ_table[3])

#TDR First mode
First_mode_BPZ_TDR_table <-mode_bootstrap(bmc_v_BPZ_TDR, seed = 1, repeats = 2000)
First_mode_BPZ_TDR_Tpod_median <- as.numeric(First_mode_BPZ_TDR_table[2])
First_mode_BPZ_TDR_Tpod_lower <-  as.numeric(First_mode_BPZ_TDR_table[1])
First_mode_BPZ_TDR_Tpod_upper <-  as.numeric(First_mode_BPZ_TDR_table[3])


# ============================================================================
# Function: 20th Gene (bootstrap)
# ============================================================================

nth_gene_bootstrap <- function(x, seed = 1, nth_gene = 20, repeats = 2000) {
  set.seed(seed)
  boot_nth_gene <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    sort(sampleData)[nth_gene]
  })
  return(quantile(boot_nth_gene, probs = c(0.025, 0.5, 0.975)))
}


nth_gene_bootstrap(bmc_v_BPZ, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_BPZ_table <- nth_gene_bootstrap(bmc_v_BPZ, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_BPZ_Tpod_median <- as.numeric(twenty_gene_BPZ_table[2])
twentieth_gene_lower_bound_BPZ = as.numeric(twenty_gene_BPZ_table[1])
twentieth_gene_upper_bound_BPZ = as.numeric(twenty_gene_BPZ_table[3])

#TDR 20th Gene
nth_gene_bootstrap(bmc_v_BPZ_TDR, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_BPZ_TDR_table <- nth_gene_bootstrap(bmc_v_BPZ_TDR, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_BPZ_TDR_Tpod_median <- as.numeric(twenty_gene_BPZ_TDR_table[2])
twentieth_gene_lower_bound_BPZ_TDR = as.numeric(twenty_gene_BPZ_TDR_table[1])
twentieth_gene_upper_bound_BPZ_TDR = as.numeric(twenty_gene_BPZ_TDR_table[3])

# ============================================================================
# Function: 10th Percentile (bootstrap)
# ============================================================================


nth_percent_bootstrap <- function(x, seed = 1, nth_percent = 10, repeats = 2000) {
  set.seed(seed)
  boot_nth_percent <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    quantile(sort(sampleData), probs = (nth_percent / 100))
  })
  return(quantile(boot_nth_percent, probs = c(0.025, 0.5, 0.975)))
}

nth_percent_bootstrap (bmc_v_BPZ, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_BPZ_table <- nth_percent_bootstrap (bmc_v_BPZ, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_BPZ_Tpod_median <- as.numeric(tenth_percentile_BPZ_table[2])

tenth_percentile_lower_bound_BPZ = as.numeric(tenth_percentile_BPZ_table[1])
tenth_percentile_upper_bound_BPZ = as.numeric(tenth_percentile_BPZ_table[3])

#TDR 10th percentile
nth_percent_bootstrap (bmc_v_BPZ_TDR, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_BPZ_TDR_table <- nth_percent_bootstrap (bmc_v_BPZ_TDR, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_BPZ_TDR_Tpod_median <- as.numeric(tenth_percentile_BPZ_TDR_table[2])

tenth_percentile_lower_bound_BPZ_TDR = as.numeric(tenth_percentile_BPZ_TDR_table[1])
tenth_percentile_upper_bound_BPZ_TDR = as.numeric(tenth_percentile_BPZ_TDR_table[3])

# ============================================================================
#Defined category Analysis Tpod (lowest median BMD)
# ============================================================================


setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Defined Category analysis")


df_BPZ_defCAT <- read.delim(
  "BPZ_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_BMD_null_DEFINED-CatMapFile- Gallus_true_true_conf0.5_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 70
)

BPZ_Cat_tpod_min <- min(df_BPZ_defCAT$BMD.Median, na.rm = TRUE)

BPZ_Cat_tpod_path_name <- df_BPZ_defCAT$GO.Pathway.Gene.Set.Gene.Name[df_BPZ_defCAT$BMD.Median == BPZ_Cat_tpod_min]

BPZ_Cat_tpod_list_bootstrap <- df_BPZ_defCAT$BMD.List[df_BPZ_defCAT$BMD.Median == BPZ_Cat_tpod_min]

BPZ_Cat_tpod_list_bootstrap_num <- as.numeric(strsplit(BPZ_Cat_tpod_list_bootstrap, ";")[[1]])
pathway_bootstrap <- function(x, seed = 1, repeats = 2000) {
  set.seed(seed)
  
  # Unlist the data once to avoid doing it repeatedly
  data_vector <- unlist(x)
  data_length <- length(data_vector)
  
  # Handle the case where there's only one value
  if (data_length == 1) {
    # Option 1: Simply return the original value with appropriate quantiles
    # Since bootstrapping a single value will always return the same value
    return(setNames(rep(data_vector, 3), c("2.5%", "50%", "97.5%")))
    
    # Option 2: If you prefer sampling behavior even for single values:
    # data_vector <- c(data_vector)  # Force it to be treated as a vector
  }
  
  # Perform the bootstrap
  boot_pathway <- replicate(repeats, {
    sampleData <- sample(data_vector, data_length, replace = TRUE)
    median(sampleData)
  })
  
  return(quantile(boot_pathway, probs = c(0.025, 0.5, 0.975)))
}

BPZ_Def_Cat_tpod_table <- pathway_bootstrap(BPZ_Cat_tpod_list_bootstrap_num, seed = 1, repeats = 2000)
BPZ_Cat_tpod_median <-as.numeric(BPZ_Def_Cat_tpod_table[2])
BPZ_Cat_tpod_lower <-as.numeric(BPZ_Def_Cat_tpod_table[1])
BPZ_Cat_tpod_upper <-as.numeric(BPZ_Def_Cat_tpod_table[3])

#TDR Defind cat
setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Defined Category analysis/top_dose_removed")
list.files()
df_D_BPZ_TDR_defCAT <- read.delim(
  "BPZ_bmdexpress_input_log2_transformed_top_dose_removed_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 83
)

D_BPZ_TDR_Cat_tpod_min <- min(df_D_BPZ_TDR_defCAT$BMD.Median, na.rm = TRUE)

D_BPZ_TDR_Cat_tpod_path_name <- df_D_BPZ_TDR_defCAT$GO.Pathway.Gene.Set.Gene.Name[df_D_BPZ_TDR_defCAT$BMD.Median == D_BPZ_TDR_Cat_tpod_min]

D_BPZ_TDR_Cat_tpod_list_bootstrap <- df_D_BPZ_TDR_defCAT$BMD.List[df_D_BPZ_TDR_defCAT$BMD.Median == D_BPZ_TDR_Cat_tpod_min]

D_BPZ_TDR_Cat_tpod_list_bootstrap_num <- as.numeric(strsplit(D_BPZ_TDR_Cat_tpod_list_bootstrap, ";")[[1]])


D_BPZ_TDR_Def_Cat_tpod_table <- pathway_bootstrap(D_BPZ_TDR_Cat_tpod_list_bootstrap_num, seed = 1, repeats = 2000)
D_BPZ_TDR_Cat_tpod_median <-as.numeric(D_BPZ_TDR_Def_Cat_tpod_table[2])
D_BPZ_TDR_Cat_tpod_lower <-as.numeric(D_BPZ_TDR_Def_Cat_tpod_table[1])
D_BPZ_TDR_Cat_tpod_upper <-as.numeric(D_BPZ_TDR_Def_Cat_tpod_table[3])

###########################################################################
#tpod summary

#BPZ_tpod_summary <- data.frame(LCRD2 = LCRD2_BPZ_Tpod, First_Mode = First_mode_BPZ_Tpod,twentieth_gene_lower_bound = twenty_gene_BPZ_table[1], twentieth_gene = twenty_gene_BPZ_Tpod,twentieth_gene_upper_bound = twenty_gene_BPZ_table[3], tenth_percentile_lower_bound = tenth_percentile_BPZ_table[1],tenth_percentile = tenth_percentile_BPZ_Tpod, tenth_percentile_upper_bound = tenth_percentile_BPZ_table[3], Lowest_BMC_median_Categorical = BPZ_Cat_tpod )
twentieth_gene_lower_bound_BPZ = as.numeric(twenty_gene_BPZ_table[1])
twentieth_gene_upper_bound_BPZ = as.numeric(twenty_gene_BPZ_table[3])
tenth_percentile_lower_bound_BPZ = as.numeric(tenth_percentile_BPZ_table[1])
tenth_percentile_upper_bound_BPZ = as.numeric(tenth_percentile_BPZ_table[3])
##########################################################################
#LC50 values:
LC50_BPZ_Lower <- (16.2-7.0)
LC50_BPZ_median <- (16.2)
LC50_BPZ_upper  <- (16.2+7.0)

#TDR LC50
LC50_BPZ_TDR_Lower <- (16.2-7.0)
LC50_BPZ_TDR_median <- (16.2)
LC50_BPZ_TDR_upper  <- (16.2+7.0)
##########################################################################
#Table
BPZ_table_full <- data.frame(
  Chemical = c("BPZ", "", "", "", "", ""),
  
  Endpoints = c(
    "LCRD",
    "First mode",
    paste0("Category (", BPZ_Cat_tpod_path_name, ")"),
    "20th gene",
    "10th percentile",
    "LC50"
  ),
  
  Lower = c(
    LCRD2_BPZ_Tpod_lower,
    First_mode_BPZ_Tpod_lower,
    BPZ_Cat_tpod_lower,
    twentieth_gene_lower_bound_BPZ,
    tenth_percentile_lower_bound_BPZ,
    LC50_BPZ_Lower
  ),
  
  Median = c(
    LCRD2_BPZ_Tpod_median,
    First_mode_BPZ_Tpod_median,
    BPZ_Cat_tpod_median,
    twenty_gene_BPZ_Tpod_median,
    tenth_percentile_BPZ_Tpod_median,
    LC50_BPZ_median
  ),
  
  Upper = c(
    LCRD2_BPZ_Tpod_upper,
    First_mode_BPZ_Tpod_upper,
    BPZ_Cat_tpod_upper,
    twentieth_gene_upper_bound_BPZ,
    tenth_percentile_upper_bound_BPZ,
    LC50_BPZ_upper
  ),
  
  Range = c(
    LCRD2_BPZ_Tpod_upper - LCRD2_BPZ_Tpod_lower,
    First_mode_BPZ_Tpod_upper - First_mode_BPZ_Tpod_lower,
    BPZ_Cat_tpod_upper - BPZ_Cat_tpod_lower,
    twentieth_gene_upper_bound_BPZ - twentieth_gene_lower_bound_BPZ,
    tenth_percentile_upper_bound_BPZ - tenth_percentile_lower_bound_BPZ,
    LC50_BPZ_upper - LC50_BPZ_Lower
  ),
  
  `DRG (Dose Responsive Genes)` = c(DRG_BPZ, "", "", "", "", ""),
  Top_Dose = c(high_dose_BPZ, "", "", "", "", ""),
  Low_Dose = c(low_dose_BPZ, "", "", "", "", "")
)

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Plots")

write.csv(BPZ_table_full, "BPZ_table_full.csv", row.names = FALSE)


#TDR table 
D_BPZ_TDR_table_full <- data.frame(
  Chemical = c("BPZ_TDR", "", "", "", "", ""),
  
  Endpoints = c(
    "LCRD",
    "First mode",
    paste0("Category (", D_BPZ_TDR_Cat_tpod_path_name, ")"),
    "20th gene",
    "10th percentile",
    "LC50"
  ),
  
  Lower = c(
    LCRD2_BPZ_TDR_Tpod_lower,
    First_mode_BPZ_TDR_Tpod_lower,
    D_BPZ_TDR_Cat_tpod_lower,
    twentieth_gene_lower_bound_BPZ_TDR,
    tenth_percentile_lower_bound_BPZ_TDR,
    LC50_BPZ_TDR_Lower
  ),
  
  Median = c(
    LCRD2_BPZ_TDR_Tpod_median,
    First_mode_BPZ_TDR_Tpod_median,
    D_BPZ_TDR_Cat_tpod_median,
    twenty_gene_BPZ_TDR_Tpod_median,
    tenth_percentile_BPZ_TDR_Tpod_median,
    LC50_BPZ_TDR_median
  ),
  
  Upper = c(
    LCRD2_BPZ_TDR_Tpod_upper,
    First_mode_BPZ_TDR_Tpod_upper,
    D_BPZ_TDR_Cat_tpod_upper,
    twentieth_gene_upper_bound_BPZ_TDR,
    tenth_percentile_upper_bound_BPZ_TDR,
    LC50_BPZ_TDR_upper
  ),
  
  Range = c(
    LCRD2_BPZ_TDR_Tpod_upper - LCRD2_BPZ_TDR_Tpod_lower,
    First_mode_BPZ_TDR_Tpod_upper - First_mode_BPZ_TDR_Tpod_lower,
    D_BPZ_TDR_Cat_tpod_upper - D_BPZ_TDR_Cat_tpod_lower,
    twentieth_gene_upper_bound_BPZ_TDR - twentieth_gene_lower_bound_BPZ_TDR,
    tenth_percentile_upper_bound_BPZ_TDR - tenth_percentile_lower_bound_BPZ_TDR,
    LC50_BPZ_TDR_upper - LC50_BPZ_TDR_Lower
  ),
  
  `DRG (Dose Responsive Genes)` = c(DRG_BPZ_TDR, "", "", "", "", ""),
  Top_Dose = c(high_dose_BPZ_TDR, "", "", "", "", ""),
  Low_Dose = c(low_dose_BPZ_TDR, "", "", "", "", ""))

#Full table merged
library(dplyr)

empty_row <- as.data.frame(
  matrix(" ", nrow = 1, ncol = ncol(BPZ_table_full))
)

merged_BPZ_FULL_Table_1 <- bind_rows(
  BPZ_table_full,
  empty_row,
  D_BPZ_TDR_table_full
)

merged_BPZ_FULL_Table <- merged_BPZ_FULL_Table_1[,1:9]

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Plots")

write.csv(merged_BPZ_FULL_Table, "BPZ_table_full_TDR_merged.csv", row.names = FALSE)


######################################################
#Density plot generation

library(dplyr)
library(ggplot2)

df <- df_BPZ %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    p10,
    LCRD2_BPZ_Tpod_median,
    First_mode_BPZ_Tpod_median,
    twenty_gene_BPZ_Tpod_median,
    BPZ_Cat_tpod_median,
    as.numeric(high_dose_BPZ),
    as.numeric(low_dose_BPZ)
  )),
  
  label = c(
    paste0("10th Percentile = ", round(p10, 3)),
    paste0("LCRD = ", round((LCRD2_BPZ_Tpod_median), 3)),
    paste0("First Mode = ", round((First_mode_BPZ_Tpod_median), 3)),
    paste0("20th Gene = ", round((twenty_gene_BPZ_Tpod_median), 3)),
    paste0("Category = ", round((BPZ_Cat_tpod_median), 3)),
    paste0("Top dose = ", as.numeric(high_dose_BPZ)),
    paste0("Low dose = ", as.numeric(low_dose_BPZ))
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
                y = c(0.02, 0.04, 0.06, 0.08, 0.10,0.12,0.14),
                label = round(value, 3),
                color = label),
            angle = 90,
            hjust = 0,
            show.legend = FALSE) +
  
  labs(x = "log10(Best.BMD)",
       y = "Density",
       color = "TPoD Methods",
       title = "                                            Distribution of Best BMD Values for BPZ",  
       subtitle = paste0(
         "                                                  DRG = ", round(DRG_BPZ, 3),
         "   |    Top Dose = ", high_dose_BPZ,
         "   |    Low Dose = ", low_dose_BPZ
       )
  )


######################################################
#TDR Plot

######################################################
#Density plot generation

library(dplyr)
library(ggplot2)

df <- df_BPZ_TDR %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    p10,
    LCRD2_BPZ_TDR_Tpod_median,
    First_mode_BPZ_TDR_Tpod_median,
    twenty_gene_BPZ_TDR_Tpod_median,
    D_BPZ_TDR_Cat_tpod_median,
    as.numeric(high_dose_BPZ_TDR),
    as.numeric(low_dose_BPZ_TDR)
  )),
  
  label = c(
    paste0("10th Percentile = ", round(p10, 3)),
    paste0("LCRD = ", round((LCRD2_BPZ_TDR_Tpod_median), 3)),
    paste0("First Mode = ", round((First_mode_BPZ_TDR_Tpod_median), 3)),
    paste0("20th Gene = ", round((twenty_gene_BPZ_TDR_Tpod_median), 3)),
    paste0("Category = ", round((D_BPZ_TDR_Cat_tpod_median), 3)),
    paste0("Top dose = ", as.numeric(high_dose_BPZ_TDR)),
    paste0("Low dose = ", as.numeric(low_dose_BPZ_TDR))
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
                y = c(0.02, 0.04, 0.06, 0.08, 0.10,0.12,0.14),
                label = round(value, 3),
                color = label),
            angle = 90,
            hjust = 0,
            show.legend = FALSE) +
  
  labs(x = "log10(Best.BMD)",
       y = "Density",
       color = "TPoD Methods",
       title = "                                            Distribution of Best BMD Values for BPZ_TDR",  
       subtitle = paste0(
         "                                                  DRG = ", round(DRG_BPZ_TDR, 3),
         "   |    Top Dose = ", high_dose_BPZ_TDR,
         "   |    Low Dose = ", low_dose_BPZ_TDR
       )
  )


















####################################################################################################
#Summary Plot

library(dplyr)
library(ggplot2)

summary_df <- data.frame(
  endpoint = c(
    "Category",
    "First Mode",
    "20th Gene",
    "10th Percentile",
    "LCRD"
  ),
  
  lowerCI = c(
    BPZ_Cat_tpod_lower,
    First_mode_BPZ_Tpod_lower,
    twentieth_gene_lower_bound_BPZ,
    tenth_percentile_lower_bound_BPZ,
    LCRD2_BPZ_Tpod_lower
  ),
  
  median = c(
    BPZ_Cat_tpod_median,
    First_mode_BPZ_Tpod_median,
    twenty_gene_BPZ_Tpod_median,
    tenth_percentile_BPZ_Tpod_median,
    LCRD2_BPZ_Tpod_median
  ),
  
  upperCI = c(
    BPZ_Cat_tpod_upper,
    First_mode_BPZ_Tpod_upper,
    twentieth_gene_upper_bound_BPZ,
    tenth_percentile_upper_bound_BPZ,
    LCRD2_BPZ_Tpod_upper
  )
)

summary_df <- summary_df %>%
  mutate(
    endpoint = factor(
      endpoint,
      levels = c(
        "Category",
        "First Mode",
        "20th Gene",
        "10th Percentile",
        "LCRD"
      )
    )
  )

ggplot(summary_df,
       aes(x = median,
           y = endpoint)) +
  
  geom_point(
    aes(colour = endpoint),
    size = 3,
    show.legend = FALSE
  ) +
  
  geom_errorbarh(
    aes(
      xmin = lowerCI,
      xmax = upperCI,
      colour = endpoint
    ),
    height = 0.3,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  theme_classic() +
  
  labs(
    title = "BPZ tPOD Summary",
    x = "tPOD",
    y = "Endpoint"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    )
  )

#################################################################################
#accumalation plot

accumulation_plot <- ggplot(
  accumulation_data,
  aes(x = logBMD, y = CumulativeCount)
)

# Add remaining plot elements
accumulation_plot <- accumulation_plot +
  
  scale_y_continuous(limits = c(0, max_y), expand = c(0, 0)) +
  
  # Accumulation curve
  geom_line(linewidth = 1.5, color = "black") +
  
  # TPOD lines
  geom_vline(
    aes(
      xintercept = log(First_mode_BPZ_Tpod_median,
                       base = params$logtransformationscale),
      color = "First Mode",
      linetype = "First Mode"
    ),
    linewidth = 0.8
  ) +
  
  geom_vline(
    aes(
      xintercept = log(BPZ_Cat_tpod_median,
                       base = params$logtransformationscale),
      color = "Category TPOD",
      linetype = "Category TPOD"
    ),
    linewidth = 0.8
  ) +
  
  geom_vline(
    aes(
      xintercept = log(twenty_gene_BPZ_Tpod_median,
                       base = params$logtransformationscale),
      color = "20th Gene",
      linetype = "20th Gene"
    ),
    linewidth = 0.8
  ) +
  
  geom_vline(
    aes(
      xintercept = log(tenth_percentile_BPZ_Tpod_median,
                       base = params$logtransformationscale),
      color = "10th Percentile",
      linetype = "10th Percentile"
    ),
    linewidth = 0.8
  ) +
  
  geom_vline(
    aes(
      xintercept = log(LCRD2_BPZ_Tpod_median,
                       base = params$logtransformationscale),
      color = "LCRD",
      linetype = "LCRD"
    ),
    linewidth = 0.8
  ) +
  
  geom_vline(
    aes(
      xintercept = log(low_dose_BPZ,
                       base = params$logtransformationscale),
      color = "Low Dose",
      linetype = "Low Dose"
    ),
    linewidth = 0.8
  ) +
  
  # Labels
  labs(
    title = paste0("BPZ (n=", nrow(raw_data_filtered), " genes)"),
    x = paste0("logBMD (", unique(metadata$DoseUnits), ")"),
    y = "Cumulative Count"
  ) +
  
  scale_color_manual(
    name = "Lines",
    values = c(
      "First Mode" = "green",
      "Category TPOD" = "orange",
      "20th Gene" = "blue",
      "10th Percentile" = "turquoise",
      "LCRD" = "hotpink",
      "Low Dose" = "grey"
    )
  ) +
  
  scale_linetype_manual(
    name = "Lines",
    values = c(
      "First Mode" = "dashed",
      "Category TPOD" = "twodash",
      "20th Gene" = "dotdash",
      "10th Percentile" = "dotted",
      "LCRD" = "longdash",
      "Low Dose" = "solid"
    )
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "top"
  )
  