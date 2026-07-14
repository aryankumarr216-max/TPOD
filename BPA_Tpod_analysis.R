####################################################################
#BPA_Analysis
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


setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic")



df_BPA <- read.delim(
  "BPA_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_BMD_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 56
)

#TDR = top dose removed

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/top_dose_removed_BDA")

df_BPA_TDR <- read.delim(
  "BPA_bmdexpress_input_log2_transformed_top_dose_removed_williams_0.05_NOMTC_foldfilter1.5_BMD_1_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 55
)

############### DRG Calculation #############################
DRG_BPA <-as.numeric(nrow(df_BPA))

DRG_BPA_TDR <-as.numeric(nrow(df_BPA_TDR))

############### Top Dose and Lose Dose Calculation ########################

setwd("C:/Users/KumarA/Downloads/Gene_count_matrices/Gene_count_matrices")

dose_BPA <- read.delim(
  "BPA_bmdexpress_input_log2_transformed.txt",
  sep = "\t",
  header = TRUE)

low_dose_BPA <- min(dose_BPA[1, ][dose_BPA[1, ] != 0], na.rm = TRUE)
high_dose_BPA <-max(as.numeric(dose_BPA[1, ]), na.rm = TRUE)

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic")

#TDR
setwd("C:/Users/KumarA/Downloads/Gene_count_matrices/Gene_count_matrices/Top_Dose_removed")

dose_D_BPA_TDR <- read.delim(
  "BPA_bmdexpress_input_log2_transformed_top_dose_removed.txt",
  sep = "\t",
  header = TRUE)

low_dose_BPA_TDR <- min(dose_D_BPA_TDR[1, ][dose_D_BPA_TDR[1, ] != 0], na.rm = TRUE)
high_dose_BPA_TDR <-max(as.numeric(dose_D_BPA_TDR[1, ]), na.rm = TRUE)

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic")


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

bmc_BPA <- as.list(df_BPA$`Best.BMD`)
bmc_v_BPA <- unlist(bmc_BPA)
probe_BPA <-as.list(df_BPA$`Probe.ID`)
probe_v_BPA <- unlist(probe_BPA)

#TDR
bmc_BPA_TDR <- as.list(df_BPA_TDR$`Best.BMD`)
bmc_v_BPA_TDR <- unlist(bmc_BPA_TDR)
probe_BPA_TDR <-as.list(df_BPA_TDR$`Probe.ID`)
probe_v_BPA_TDR <- unlist(probe_BPA_TDR)



LCRD2_BPA_Tpod_table_non_boot_strap <-LCRD.2(bmc_v_BPA,probe_v_BPA)
LCRD2_BPA_Tpod_non_boot_strap <- LCRD2_BPA_Tpod_table$BMC


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

bmds <- df_BPA[["Best.BMD"]]
names(bmds) <- df_BPA[["Probe.ID"]]

repeats <- 2000
lcrdratiocut <- 1.778
lcrdlogbase <- 10

boot_lcrd <- replicate(repeats, {
  sampleData <- sample(unlist(bmds), length(unlist(bmds)), replace = TRUE)
  result <- LCRD(
    bmc = sampleData,
    probe = names(sampleData),
    cut = lcrdratiocut,
    logbase = lcrdlogbase
  )
  result$LCRD_Result$logBMC
})

summary(boot_lcrd)
table(boot_lcrd)
length(unique(boot_lcrd))


LCRD_BPA_Tpod_table<- lcrd_bootstrap(df_BPA, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10)
LCRD2_BPA_Tpod_median <-10^(as.numeric(LCRD_BPA_Tpod_table[2]))
LCRD2_BPA_Tpod_lower<-10^(as.numeric(LCRD_BPA_Tpod_table[1]))
LCRD2_BPA_Tpod_upper<-10^(as.numeric(LCRD_BPA_Tpod_table[3]))


#TDR LCRD
LCRD_BPA_TDR_Tpod_table<- lcrd_bootstrap(df_BPA_TDR, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10)
LCRD2_BPA_TDR_Tpod_median <-10^(as.numeric(LCRD_BPA_TDR_Tpod_table[2]))
LCRD2_BPA_TDR_Tpod_lower<-10^(as.numeric(LCRD_BPA_TDR_Tpod_table[1]))
LCRD2_BPA_TDR_Tpod_upper<-10^(as.numeric(LCRD_BPA_TDR_Tpod_table[3]))

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

First.Mode(bmc_v_BPA)
First_mode_BPA_table_non_boot_strap <- First.Mode(bmc_v_BPA)
First_mode_BPA_Tpod_non_boot_strap <- First_mode_BPA_table$First_Mode


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



First_mode_BPA_table <-mode_bootstrap(bmc_v_BPA, seed = 1, repeats = 2000)

First_mode_BPA_Tpod_median <- as.numeric(First_mode_BPA_table[2])
First_mode_BPA_Tpod_lower <-  as.numeric(First_mode_BPA_table[1])
First_mode_BPA_Tpod_upper <-  as.numeric(First_mode_BPA_table[3])

#TDR First mode
First_mode_BPA_TDR_table <-mode_bootstrap(bmc_v_BPA_TDR, seed = 1, repeats = 2000)
First_mode_BPA_TDR_Tpod_median <- as.numeric(First_mode_BPA_TDR_table[2])
First_mode_BPA_TDR_Tpod_lower <-  as.numeric(First_mode_BPA_TDR_table[1])
First_mode_BPA_TDR_Tpod_upper <-  as.numeric(First_mode_BPA_TDR_table[3])

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


nth_gene_bootstrap(bmc_v_BPA, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_BPA_table <- nth_gene_bootstrap(bmc_v_BPA, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_BPA_Tpod_median <- as.numeric(twenty_gene_BPA_table[2])
twentieth_gene_lower_bound_BPA = as.numeric(twenty_gene_BPA_table[1])
twentieth_gene_upper_bound_BPA = as.numeric(twenty_gene_BPA_table[3])

#TDR 20th Gene
nth_gene_bootstrap(bmc_v_BPA_TDR, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_BPA_TDR_table <- nth_gene_bootstrap(bmc_v_BPA_TDR, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_BPA_TDR_Tpod_median <- as.numeric(twenty_gene_BPA_TDR_table[2])
twentieth_gene_lower_bound_BPA_TDR = as.numeric(twenty_gene_BPA_TDR_table[1])
twentieth_gene_upper_bound_BPA_TDR = as.numeric(twenty_gene_BPA_TDR_table[3])


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

nth_percent_bootstrap (bmc_v_BPA, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_BPA_table <- nth_percent_bootstrap (bmc_v_BPA, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_BPA_Tpod_median <- as.numeric(tenth_percentile_BPA_table[2])

tenth_percentile_lower_bound_BPA = as.numeric(tenth_percentile_BPA_table[1])
tenth_percentile_upper_bound_BPA = as.numeric(tenth_percentile_BPA_table[3])

#TDR 10th percentile
nth_percent_bootstrap (bmc_v_BPA_TDR, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_BPA_TDR_table <- nth_percent_bootstrap (bmc_v_BPA_TDR, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_BPA_TDR_Tpod_median <- as.numeric(tenth_percentile_BPA_TDR_table[2])

tenth_percentile_lower_bound_BPA_TDR = as.numeric(tenth_percentile_BPA_TDR_table[1])
tenth_percentile_upper_bound_BPA_TDR = as.numeric(tenth_percentile_BPA_TDR_table[3])

# ============================================================================
#Defined category Analysis Tpod (lowest median BMD)
# ============================================================================

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Defined Category analysis")


df_BPA_defCAT <- read.delim(
  "BPA_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_BMD_null_DEFINED-CatMapFile- Gallus_true_true_conf0.5_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 70
)

BPA_Cat_tpod_min <- min(df_BPA_defCAT$BMD.Median, na.rm = TRUE)

BPA_Cat_tpod_path_name <- df_BPA_defCAT$GO.Pathway.Gene.Set.Gene.Name[df_BPA_defCAT$BMD.Median == BPA_Cat_tpod_min]

BPA_Cat_tpod_list_bootstrap <- df_BPA_defCAT$BMD.List[df_BPA_defCAT$BMD.Median == BPA_Cat_tpod_min]

BPA_Cat_tpod_list_bootstrap_num <- as.numeric(strsplit(BPA_Cat_tpod_list_bootstrap, ";")[[1]])
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

BPA_Def_Cat_tpod_table <- pathway_bootstrap(BPA_Cat_tpod_list_bootstrap_num, seed = 1, repeats = 2000)
BPA_Cat_tpod_median <-as.numeric(BPA_Def_Cat_tpod_table[2])
BPA_Cat_tpod_lower <-as.numeric(BPA_Def_Cat_tpod_table[1])
BPA_Cat_tpod_upper <-as.numeric(BPA_Def_Cat_tpod_table[3])


#TDR Defind cat
setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Defined Category analysis/top_dose_removed")
list.files()
df_D_BPA_TDR_defCAT <- read.delim(
  "BPA_bmdexpress_input_log2_transformed_top_dose_removed_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 79
)

D_BPA_TDR_Cat_tpod_min <- min(df_D_BPA_TDR_defCAT$BMD.Median, na.rm = TRUE)

D_BPA_TDR_Cat_tpod_path_name <- df_D_BPA_TDR_defCAT$GO.Pathway.Gene.Set.Gene.Name[df_D_BPA_TDR_defCAT$BMD.Median == D_BPA_TDR_Cat_tpod_min]

D_BPA_TDR_Cat_tpod_list_bootstrap <- df_D_BPA_TDR_defCAT$BMD.List[df_D_BPA_TDR_defCAT$BMD.Median == D_BPA_TDR_Cat_tpod_min]

D_BPA_TDR_Cat_tpod_list_bootstrap_num <- as.numeric(strsplit(D_BPA_TDR_Cat_tpod_list_bootstrap, ";")[[1]])


D_BPA_TDR_Def_Cat_tpod_table <- pathway_bootstrap(D_BPA_TDR_Cat_tpod_list_bootstrap_num, seed = 1, repeats = 2000)
D_BPA_TDR_Cat_tpod_median <-as.numeric(D_BPA_TDR_Def_Cat_tpod_table[2])
D_BPA_TDR_Cat_tpod_lower <-as.numeric(D_BPA_TDR_Def_Cat_tpod_table[1])
D_BPA_TDR_Cat_tpod_upper <-as.numeric(D_BPA_TDR_Def_Cat_tpod_table[3])
#####################################################################################
#LC50 values:
LC50_BPA_Lower <- (61.7-43.1)
LC50_BPA_median <- (61.7)
LC50_BPA_upper  <- (61.7+43.1)

#TDR LC50
LC50_BPA_TDR_Lower <- (61.7-43.1)
LC50_BPA_TDR_median <- (61.7)
LC50_BPA_TDR_upper  <- (61.7+43.1)
##########################################################################
#Table
BPA_table_full <- data.frame(
  Chemical = c("BPA", "", "", "", "", ""),
  
  Endpoints = c(
    "LCRD",
    "First mode",
    paste0("Category (", BPA_Cat_tpod_path_name, ")"),
    "20th gene",
    "10th percentile",
    "LC50"
  ),
  
  Lower = c(
    LCRD2_BPA_Tpod_lower,
    First_mode_BPA_Tpod_lower,
    BPA_Cat_tpod_lower,
    twentieth_gene_lower_bound_BPA,
    tenth_percentile_lower_bound_BPA,
    LC50_BPA_Lower
  ),
  
  Median = c(
    LCRD2_BPA_Tpod_median,
    First_mode_BPA_Tpod_median,
    BPA_Cat_tpod_median,
    twenty_gene_BPA_Tpod_median,
    tenth_percentile_BPA_Tpod_median,
    LC50_BPA_median
  ),
  
  Upper = c(
    LCRD2_BPA_Tpod_upper,
    First_mode_BPA_Tpod_upper,
    BPA_Cat_tpod_upper,
    twentieth_gene_upper_bound_BPA,
    tenth_percentile_upper_bound_BPA,
    LC50_BPA_upper
  ),
  
  Range = c(
    LCRD2_BPA_Tpod_upper - LCRD2_BPA_Tpod_lower,
    First_mode_BPA_Tpod_upper - First_mode_BPA_Tpod_lower,
    BPA_Cat_tpod_upper - BPA_Cat_tpod_lower,
    twentieth_gene_upper_bound_BPA - twentieth_gene_lower_bound_BPA,
    tenth_percentile_upper_bound_BPA - tenth_percentile_lower_bound_BPA,
    LC50_BPA_upper - LC50_BPA_Lower
  ),
  
  `DRG (Dose Responsive Genes)` = c(DRG_BPA, "", "", "", "", ""),
  Top_Dose = c(high_dose_BPA, "", "", "", "", ""),
  Low_Dose = c(low_dose_BPA, "", "", "", "", "")
)

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Plots")

write.csv(BPA_table_full, "BPA_table_full.csv", row.names = FALSE)



#TDR table 
D_BPA_TDR_table_full <- data.frame(
  Chemical = c("BPA_TDR", "", "", "", "", ""),
  
  Endpoints = c(
    "LCRD",
    "First mode",
    paste0("Category (", D_BPA_TDR_Cat_tpod_path_name, ")"),
    "20th gene",
    "10th percentile",
    "LC50"
  ),
  
  Lower = c(
    LCRD2_BPA_TDR_Tpod_lower,
    First_mode_BPA_TDR_Tpod_lower,
    D_BPA_TDR_Cat_tpod_lower,
    twentieth_gene_lower_bound_BPA_TDR,
    tenth_percentile_lower_bound_BPA_TDR,
    LC50_BPA_TDR_Lower
  ),
  
  Median = c(
    LCRD2_BPA_TDR_Tpod_median,
    First_mode_BPA_TDR_Tpod_median,
    D_BPA_TDR_Cat_tpod_median,
    twenty_gene_BPA_TDR_Tpod_median,
    tenth_percentile_BPA_TDR_Tpod_median,
    LC50_BPA_TDR_median
  ),
  
  Upper = c(
    LCRD2_BPA_TDR_Tpod_upper,
    First_mode_BPA_TDR_Tpod_upper,
    D_BPA_TDR_Cat_tpod_upper,
    twentieth_gene_upper_bound_BPA_TDR,
    tenth_percentile_upper_bound_BPA_TDR,
    LC50_BPA_TDR_upper
  ),
  
  Range = c(
    LCRD2_BPA_TDR_Tpod_upper - LCRD2_BPA_TDR_Tpod_lower,
    First_mode_BPA_TDR_Tpod_upper - First_mode_BPA_TDR_Tpod_lower,
    D_BPA_TDR_Cat_tpod_upper - D_BPA_TDR_Cat_tpod_lower,
    twentieth_gene_upper_bound_BPA_TDR - twentieth_gene_lower_bound_BPA_TDR,
    tenth_percentile_upper_bound_BPA_TDR - tenth_percentile_lower_bound_BPA_TDR,
    LC50_BPA_TDR_upper - LC50_BPA_TDR_Lower
  ),
  
  `DRG (Dose Responsive Genes)` = c(DRG_BPA_TDR, "", "", "", "", ""),
  Top_Dose = c(high_dose_BPA_TDR, "", "", "", "", ""),
  Low_Dose = c(low_dose_BPA_TDR, "", "", "", "", ""))

#Full table merged
library(dplyr)

empty_row <- as.data.frame(
  matrix(" ", nrow = 1, ncol = ncol(BPA_table_full))
)

merged_BPA_FULL_Table_1 <- bind_rows(
  BPA_table_full,
  empty_row,
  D_BPA_TDR_table_full
)

merged_BPA_FULL_Table <- merged_BPA_FULL_Table_1[,1:9]

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Plots")

write.csv(merged_BPA_FULL_Table, "BPA_table_full_TDR_merged.csv", row.names = FALSE)

##########################################################################
#plot

######################################################
#Density plot generation

library(dplyr)
library(ggplot2)

p10_BPA <- tenth_percentile_BPA_Tpod_median

df <- df_BPA %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    p10_BPA,
    LCRD2_BPA_Tpod_median,
    First_mode_BPA_Tpod_median,
    twenty_gene_BPA_Tpod_median,
    BPA_Cat_tpod_median,
    as.numeric(high_dose_BPA),
    as.numeric(low_dose_BPA)
  )),
  
  label = c(
    paste0("10th Percentile = ", round(p10_BPA, 3)),
    paste0("LCRD = ", round((LCRD2_BPA_Tpod_median), 3)),
    paste0("First Mode = ", round((First_mode_BPA_Tpod_median), 3)),
    paste0("20th Gene = ", round((twenty_gene_BPA_Tpod_median), 3)),
    paste0("Category = ", round((D_BPA_TDR_Cat_tpod_median), 3)),
    paste0("Top dose = ", as.numeric(high_dose_BPA)),
    paste0("Low dose = ", as.numeric(low_dose_BPA))
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
       title = "                                            Distribution of Best BMD Values for BPA",  
       subtitle = paste0(
         "                                                  DRG = ", round(DRG_BPA, 3),
         "   |    Top Dose = ", high_dose_BPA,
         "   |    Low Dose = ", low_dose_BPA
       )
  )




###############################################################################

#TDR Plot

######################################################
#Density plot generation

library(dplyr)
library(ggplot2)

p10_BPA_TDR <- tenth_percentile_BPA_TDR_Tpod_median


df <- df_BPA_TDR %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    p10_BPA_TDR,
    LCRD2_BPA_TDR_Tpod_median,
    First_mode_BPA_TDR_Tpod_median,
    twenty_gene_BPA_TDR_Tpod_median,
    D_BPA_TDR_Cat_tpod_median,
    as.numeric(high_dose_BPA_TDR),
    as.numeric(low_dose_BPA_TDR)
  )),
  
  label = c(
    paste0("10th Percentile = ", round(p10_BPA_TDR, 3)),
    paste0("LCRD = ", round((LCRD2_BPA_TDR_Tpod_median), 3)),
    paste0("First Mode = ", round((First_mode_BPA_TDR_Tpod_median), 3)),
    paste0("20th Gene = ", round((twenty_gene_BPA_TDR_Tpod_median), 3)),
    paste0("Category = ", round((D_BPA_TDR_Cat_tpod_median), 3)),
    paste0("Top dose = ", as.numeric(high_dose_BPA_TDR)),
    paste0("Low dose = ", as.numeric(low_dose_BPA_TDR))
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
       title = "                                            Distribution of Best BMD Values for BPA_TDR",  
       subtitle = paste0(
         "                                                  DRG = ", round(DRG_BPA_TDR, 3),
         "   |    Top Dose = ", high_dose_BPA_TDR,
         "   |    Low Dose = ", low_dose_BPA_TDR
       )
  )


###############################################################################
#Density plot generation - updated tpods - normilized 


setwd("C:/Users/KumarA/OneDrive - EC-EC/Documents/Updated_version2_BMD_RUN_FILES/BPA_Updated_norm_version2!Count_files_TDR!Updated_version_2_BPA_Alt_median_ratio_BPA!full_data_4.4m")
list.files()

df_updated_all_topds_BPA <- read.delim(
  "outliers_removed_lowcount_filtered_Iter_1_tpod.txt",
)



#first mdoe
#updated_BPA_first_mode_lower_tpod <- min(df_updated_all_topds_BPA$first_mode, na.rm = TRUE)
updated_BPA_first_mode_median_tpod <- df_updated_all_topds_BPA$first_mode
#updated_BPA_first_mode_upper_tpod <- max(df_updated_all_topds_BPA$first_mode, na.rm = TRUE)

# tenth_percentile
#updated_BPA_tenth_percentile_lower_tpod <- min(df_updated_all_topds_BPA$tenth_percentile, na.rm = TRUE)
updated_BPA_tenth_percentile_median_tpod <-df_updated_all_topds_BPA$tenth_percentile
#updated_BPA_tenth_percentile_upper_tpod <- max(df_updated_all_topds_BPA$tenth_percentile, na.rm = TRUE)

# twentieth_gene
#updated_BPA_twentieth_gene_lower_tpod <- min(df_updated_all_topds_BPA$twentieth_gene, na.rm = TRUE)
updated_BPA_twentieth_gene_median_tpod <- df_updated_all_topds_BPA$twentieth_gene
#updated_BPA_twentieth_gene_upper_tpod <- max(df_updated_all_topds_BPA$twentieth_gene, na.rm = TRUE)

# LCRD2
#updated_BPA_LCRD2_lower_tpod <- min(df_updated_all_topds_BPA$lcrd2, na.rm = TRUE)
updated_BPA_LCRD2_median_tpod <- df_updated_all_topds_BPA$lcrd2
#updated_BPA_LCRD2_upper_tpod <- max(df_updated_all_topds_BPA$lcrd2, na.rm = TRUE)

updated_drg_mean_BPA <- as.numeric(df_updated_all_topds_BPA$drgs)
############### TDR #####################

setwd("C:/Users/KumarA/OneDrive - EC-EC/Documents/Updated_version2_BMD_RUN_FILES/TDR/TDR_BPA_Updated_norm_version2!Count_files_TDR!Updated_version_2_BPA_Alt_median_ratio_BPA!full_data_4.4m")
list.files()

df_updated_all_topds_TDR_BPA <- read.delim(
  "outliers_removed_lowcount_filtered_Iter_1_tpod.txt",
)



#first mdoe
updated_TDR_BPA_first_mode_median_tpod <- df_updated_all_topds_TDR_BPA$first_mode

# tenth_percentile
updated_TDR_BPA_tenth_percentile_median_tpod <- df_updated_all_topds_TDR_BPA$tenth_percentile

# twentieth_gene
updated_TDR_BPA_twentieth_gene_median_tpod <- df_updated_all_topds_TDR_BPA$twentieth_gene

# LCRD2
updated_TDR_BPA_LCRD2_median_tpod <- df_updated_all_topds_TDR_BPA$lcrd2

updated_TDR_drg_mean_BPA <- as.numeric(df_updated_all_topds_TDR_BPA$drgs)


############### updated bmd file #####################


df_BPA <- read.delim(
  "outliers_removed_lowcount_filtered_Iter_1_Bmds_Analysis.txt",
  sep = "\t",
  header = TRUE,
)

############### TDR updated bmd file #####################

setwd("C:/Users/KumarA/OneDrive - EC-EC/Documents/Updated_version2_BMD_RUN_FILES/TDR/TDR_BPA_Updated_norm_version2!Count_files_TDR!Updated_version_2_BPA_Alt_median_ratio_BPA!full_data_4.4m")

df_TDR_BPA <- read.delim(
  "outliers_removed_lowcount_filtered_Iter_1_Bmds_Analysis.txt",
  sep = "\t",
  header = TRUE,
)

###################################
#normal data plot
library(dplyr)
library(ggplot2)

D_BPA_Cat_tpod_median <- 0
p10_BPA <- tenth_percentile_BPA_Tpod_median

df <- df_BPA %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    as.numeric(updated_BPA_tenth_percentile_median_tpod),
    as.numeric(updated_BPA_LCRD2_median_tpod),
    as.numeric(updated_BPA_first_mode_median_tpod),
    as.numeric(updated_BPA_twentieth_gene_median_tpod),
    D_BPA_Cat_tpod_median,
    as.numeric(high_dose_BPA),
    as.numeric(low_dose_BPA)
  )),
  
  label = c(
    paste0("10th Percentile = ", round(updated_BPA_tenth_percentile_median_tpod, 3)),
    paste0("LCRD = ", round(updated_BPA_LCRD2_median_tpod, 3)),
    paste0("First Mode = ", round(updated_BPA_first_mode_median_tpod, 3)),
    paste0("20th Gene = ", round(updated_BPA_twentieth_gene_median_tpod, 3)),
    paste0("Category = ", round(D_BPA_Cat_tpod_median, 3)),
    paste0("Top dose = ", as.numeric(high_dose_BPA)),
    paste0("Low dose = ", as.numeric(low_dose_BPA))
  )
)

# Histogram bin width
bw <- 0.2

ggplot(df, aes(x = logBestBMD)) +
  
  # Histogram (counts)
  geom_histogram(
    aes(y = after_stat(count)),
    binwidth = bw,
    fill = "lightblue",
    color = "black",
    alpha = 0.5
  ) +
  
  # Density curve scaled to counts
  geom_density(
    aes(y = after_stat(density * nrow(df) * bw)),
    color = "blue",
    linewidth = 1.2
  ) +
  
  # TPoD lines
  geom_vline(
    data = tpod_df,
    aes(
      xintercept = value,
      color = label
    ),
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # TPoD labels
  geom_text(
    data = tpod_df,
    aes(
      x = value + 0.03,
      y = Inf,
      label = round(value, 3),
      color = label
    ),
    angle = 90,
    vjust = 1.2,
    show.legend = FALSE
  ) +
  
  labs(
    x = "log10(Best.BMD)",
    y = "Count",
    color = "TPoD Methods",
    title = "Distribution of Best BMD Values for BPA",
    subtitle = paste0(
      "DRG = ", round(updated_drg_mean_BPA, 3),
      " | Top Dose = ", high_dose_BPA,
      " | Low Dose = ", low_dose_BPA
    )
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

####################################################
#TDR Plot
#normal data plot
library(dplyr)
library(ggplot2)

D_TDR_BPA_Cat_tpod_median <- 0
p10_TDR_BPA <- tenth_percentile_TDR_BPA_Tpod_median

df <- df_TDR_BPA %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    as.numeric(updated_TDR_BPA_tenth_percentile_median_tpod),
    as.numeric(updated_TDR_BPA_LCRD2_median_tpod),
    as.numeric(updated_TDR_BPA_first_mode_median_tpod),
    as.numeric(updated_TDR_BPA_twentieth_gene_median_tpod),
    D_TDR_BPA_Cat_tpod_median,
    as.numeric(high_dose_BPA_TDR),
    as.numeric(low_dose_BPA_TDR)
  )),
  
  label = c(
    paste0("10th Percentile = ", round(updated_TDR_BPA_tenth_percentile_median_tpod, 3)),
    paste0("LCRD = ", round(updated_TDR_BPA_LCRD2_median_tpod, 3)),
    paste0("First Mode = ", round(updated_TDR_BPA_first_mode_median_tpod, 3)),
    paste0("20th Gene = ", round(updated_TDR_BPA_twentieth_gene_median_tpod, 3)),
    paste0("Category = ", round(D_TDR_BPA_Cat_tpod_median, 3)),
    paste0("Top dose = ", as.numeric(high_dose_BPA_TDR)),
    paste0("Low dose = ", as.numeric(low_dose_BPA_TDR))
  )
)

# Histogram bin width
bw <- 0.2

ggplot(df, aes(x = logBestBMD)) +
  
  # Histogram (counts)
  geom_histogram(
    aes(y = after_stat(count)),
    binwidth = bw,
    fill = "lightblue",
    color = "black",
    alpha = 0.5
  ) +
  
  # Density curve scaled to counts
  geom_density(
    aes(y = after_stat(density * nrow(df) * bw)),
    color = "blue",
    linewidth = 1.2
  ) +
  
  # TPoD lines
  geom_vline(
    data = tpod_df,
    aes(
      xintercept = value,
      color = label
    ),
    linetype = "dashed",
    linewidth = 1
  ) +
  
  # TPoD labels
  geom_text(
    data = tpod_df,
    aes(
      x = value + 0.03,
      y = Inf,
      label = round(value, 3),
      color = label
    ),
    angle = 90,
    vjust = 1.2,
    show.legend = FALSE
  ) +
  
  labs(
    x = "log10(Best.BMD)",
    y = "Count",
    color = "TPoD Methods",
    title = "Distribution of Best BMD Values for TDR BPA",
    subtitle = paste0(
      "DRG = ", round(updated_TDR_drg_mean_BPA, 3),
      " | Top Dose = ", high_dose_BPA_TDR,
      " | Low Dose = ", low_dose_BPA_TDR
    )
  ) +
  
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

updated_drg_mean_BPA
updated_TDR_drg_mean_BPA






