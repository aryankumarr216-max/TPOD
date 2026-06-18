####################################################################
#24BPF_Analysis
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



df_24BPF <- read.delim(
  "24BPF_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_BMD_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 56
)
#TDR = top dose removed

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/top_dose_removed_BDA")

df_24BPF_TDR <- read.delim(
  "24BPF_bmdexpress_input_log2_transformed_top_dose_removed_williams_0.05_NOMTC_foldfilter1.5_BMD.txt",
  sep = "\t",
  header = TRUE,
  skip = 50
)

############### DRG Calculation #############################
DRG_24BPF <-as.numeric(nrow(df_24BPF))

DRG_24BPF_TDR <-as.numeric(nrow(df_24BPF_TDR))


############### Top Dose and Lose Dose Calculation ########################
setwd("C:/Users/KumarA/Downloads/Gene_count_matrices/Gene_count_matrices")

dose_D_24BPF <- read.delim(
  "24BPF_bmdexpress_input_log2_transformed.txt",
  sep = "\t",
  header = TRUE)

low_dose_24BPF <- min(dose_D_24BPF[1, ][dose_D_24BPF[1, ] != 0], na.rm = TRUE)
high_dose_24BPF <-max(as.numeric(dose_D_24BPF[1, ]), na.rm = TRUE)

setwd("C:/Users/KumarA/Downloads/Gene_count_matrices/Gene_count_matrices/Top_Dose_removed")

dose_D_24BPF_TDR <- read.delim(
  "24BPF_bmdexpress_input_log2_transformed_top_dose_removed.txt",
  sep = "\t",
  header = TRUE)

low_dose_24BPF_TDR <- min(dose_D_24BPF_TDR[1, ][dose_D_24BPF_TDR[1, ] != 0], na.rm = TRUE)
high_dose_24BPF_TDR <-max(as.numeric(dose_D_24BPF_TDR[1, ]), na.rm = TRUE)

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



bmc_24BPF <- as.list(df_24BPF$`Best.BMD`)
bmc_v_24BPF <- unlist(bmc_24BPF)
probe_24BPF <-as.list(df_24BPF$`Probe.ID`)
probe_v_24BPF <- unlist(probe_24BPF)

#TDR
bmc_24BPF_TDR <- as.list(df_24BPF_TDR$`Best.BMD`)
bmc_v_24BPF_TDR <- unlist(bmc_24BPF_TDR)
probe_24BPF_TDR <-as.list(df_24BPF_TDR$`Probe.ID`)
probe_v_24BPF_TDR <- unlist(probe_24BPF_TDR)


LCRD2_24BPF_Tpod_table_non_boot_strap <- LCRD.2(bmc_v_24BPF,probe_v_24BPF)
LCRD2_24BPF_Tpod_non_boot_strap <- as.numeric(LCRD.2(bmc_v_24BPF,probe_v_24BPF)[2])



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

bmds <- df_24BPF[["Best.BMD"]]
names(bmds) <- df_24BPF[["Probe.ID"]]

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

lcrd_bootstrap(df_24BPF, seed = 1, repeats = 1, lcrdratiocut = 1.778, lcrdlogbase = 10)
LCRD_24BPF_Tpod_table<- lcrd_bootstrap(df_24BPF, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10)
LCRD2_24BPF_Tpod_median <-10^(as.numeric(LCRD_24BPF_Tpod_table[2]))
LCRD2_24BPF_Tpod_lower<-10^(as.numeric(LCRD_24BPF_Tpod_table[1]))
LCRD2_24BPF_Tpod_upper<-10^(as.numeric(LCRD_24BPF_Tpod_table[3]))

#TDR LCRD
LCRD_24BPF_TDR_Tpod_table<- lcrd_bootstrap(df_24BPF_TDR, seed = 1, repeats = 2000, lcrdratiocut = 1.778, lcrdlogbase = 10)
LCRD2_24BPF_TDR_Tpod_median <-10^(as.numeric(LCRD_24BPF_TDR_Tpod_table[2]))
LCRD2_24BPF_TDR_Tpod_lower<-10^(as.numeric(LCRD_24BPF_TDR_Tpod_table[1]))
LCRD2_24BPF_TDR_Tpod_upper<-10^(as.numeric(LCRD_24BPF_TDR_Tpod_table[3]))

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

First_mode_24BPF_table_non_boot_strap<-First.Mode(bmc_v_24BPF)
First_mode_24BPF_Tpod_non_boot_strap <- First_mode_24BPF_table$First_Mode

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

#Old function
#mode_bootstrap <- function(x, seed = 1, repeats = 2000, min_size = 0.06, min_bw = 0.01) {
  #source("mode_antimode.R")
  set.seed(seed)
  boot_mode <- replicate(repeats, {
    sampleData <- sample(unlist(x), length(unlist(x)), replace = TRUE)
    dataMode <- mode.antimode(sampleData, min.size = min_size, bw = "SJ", min.bw = min_bw)
    dataMode$modes[[1]]
  })
  return(quantile(boot_mode, probs = c(0.025, 0.5, 0.975)))



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



First_mode_24BPF_table <-mode_bootstrap(bmc_v_24BPF, seed = 1, repeats = 2000)
First_mode_24BPF_Tpod_median <- as.numeric(First_mode_24BPF_table[2])
First_mode_24BPF_Tpod_lower <-  as.numeric(First_mode_24BPF_table[1])
First_mode_24BPF_Tpod_upper <-  as.numeric(First_mode_24BPF_table[3])

#TDR First mode
First_mode_24BPF_TDR_table <-mode_bootstrap(bmc_v_24BPF_TDR, seed = 1, repeats = 2000)
First_mode_24BPF_TDR_Tpod_median <- as.numeric(First_mode_24BPF_TDR_table[2])
First_mode_24BPF_TDR_Tpod_lower <-  as.numeric(First_mode_24BPF_TDR_table[1])
First_mode_24BPF_TDR_Tpod_upper <-  as.numeric(First_mode_24BPF_TDR_table[3])

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


nth_gene_bootstrap(bmc_v_24BPF, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_24BPF_table <- nth_gene_bootstrap(bmc_v_24BPF, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_24BPF_Tpod_median <- as.numeric(twenty_gene_24BPF_table[2])
twentieth_gene_lower_bound_24BPF = as.numeric(twenty_gene_24BPF_table[1])
twentieth_gene_upper_bound_24BPF = as.numeric(twenty_gene_24BPF_table[3])

#TDR 20th Gene
nth_gene_bootstrap(bmc_v_24BPF_TDR, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_24BPF_TDR_table <- nth_gene_bootstrap(bmc_v_24BPF_TDR, seed = 1, nth_gene = 20, repeats = 2000)
twenty_gene_24BPF_TDR_Tpod_median <- as.numeric(twenty_gene_24BPF_TDR_table[2])
twentieth_gene_lower_bound_24BPF_TDR = as.numeric(twenty_gene_24BPF_TDR_table[1])
twentieth_gene_upper_bound_24BPF_TDR = as.numeric(twenty_gene_24BPF_TDR_table[3])

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

nth_percent_bootstrap (bmc_v_24BPF, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_24BPF_table <- nth_percent_bootstrap (bmc_v_24BPF, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_24BPF_Tpod_median <- as.numeric(tenth_percentile_24BPF_table[2])

tenth_percentile_lower_bound_24BPF = as.numeric(tenth_percentile_24BPF_table[1])
tenth_percentile_upper_bound_24BPF = as.numeric(tenth_percentile_24BPF_table[3])

#TDR 10th percentile
nth_percent_bootstrap (bmc_v_24BPF_TDR, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_24BPF_TDR_table <- nth_percent_bootstrap (bmc_v_24BPF_TDR, seed = 1, nth_percent = 10, repeats = 2000)
tenth_percentile_24BPF_TDR_Tpod_median <- as.numeric(tenth_percentile_24BPF_TDR_table[2])

tenth_percentile_lower_bound_24BPF_TDR = as.numeric(tenth_percentile_24BPF_TDR_table[1])
tenth_percentile_upper_bound_24BPF_TDR = as.numeric(tenth_percentile_24BPF_TDR_table[3])

# ============================================================================
#Defined category Analysis Tpod (lowest median BMD)
# ============================================================================

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Defined Category analysis")


df_D_24BPF_defCAT <- read.delim(
  "24BPF_bmdexpress_input_log2_transformed_williams_0.05_NOMTC_foldfilter1.5_BMD_null_DEFINED-CatMapFile- Gallus_true_true_conf0.5_filtered.txt",
  sep = "\t",
  header = TRUE,
  skip = 70
)


D_24BPF_Cat_tpod_min <- min(df_D_24BPF_defCAT$BMD.Median, na.rm = TRUE)

D_24BPF_Cat_tpod_path_name <- df_D_24BPF_defCAT$GO.Pathway.Gene.Set.Gene.Name[df_D_24BPF_defCAT$BMD.Median == D_24BPF_Cat_tpod_min]

D_24BPF_Cat_tpod_list_bootstrap <- df_D_24BPF_defCAT$BMD.List[df_D_24BPF_defCAT$BMD.Median == D_24BPF_Cat_tpod_min]

D_24BPF_Cat_tpod_list_bootstrap_num <- as.numeric(strsplit(D_24BPF_Cat_tpod_list_bootstrap, ";")[[1]])
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

D_24BPF_Def_Cat_tpod_table <- pathway_bootstrap(D_24BPF_Cat_tpod_list_bootstrap_num, seed = 1, repeats = 2000)
D_24BPF_Cat_tpod_median <-as.numeric(D_24BPF_Def_Cat_tpod_table[2])
D_24BPF_Cat_tpod_lower <-as.numeric(D_24BPF_Def_Cat_tpod_table[1])
D_24BPF_Cat_tpod_upper <-as.numeric(D_24BPF_Def_Cat_tpod_table[3])

#TDR Defind cat
setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Defined Category analysis/top_dose_removed")

df_D_24BPF_TDR_defCAT <- read.delim(
  "24BPF_bmdexpress_input_log2_transformed_top_dose_removed_williams_0.05_NOMTC_foldfilter1.5_BMD_null_DEFINED-CatMapFile.txt",
  sep = "\t",
  header = TRUE,
  skip = 79
)

D_24BPF_TDR_Cat_tpod_min <- min(df_D_24BPF_TDR_defCAT$BMD.Median, na.rm = TRUE)

D_24BPF_TDR_Cat_tpod_path_name <- df_D_24BPF_TDR_defCAT$GO.Pathway.Gene.Set.Gene.Name[df_D_24BPF_TDR_defCAT$BMD.Median == D_24BPF_TDR_Cat_tpod_min]

D_24BPF_TDR_Cat_tpod_list_bootstrap <- df_D_24BPF_TDR_defCAT$BMD.List[df_D_24BPF_TDR_defCAT$BMD.Median == D_24BPF_TDR_Cat_tpod_min]

D_24BPF_TDR_Cat_tpod_list_bootstrap_num <- as.numeric(strsplit(D_24BPF_TDR_Cat_tpod_list_bootstrap, ";")[[1]])


D_24BPF_TDR_Def_Cat_tpod_table <- pathway_bootstrap(D_24BPF_TDR_Cat_tpod_list_bootstrap_num, seed = 1, repeats = 2000)
D_24BPF_TDR_Cat_tpod_median <-as.numeric(D_24BPF_TDR_Def_Cat_tpod_table[2])
D_24BPF_TDR_Cat_tpod_lower <-as.numeric(D_24BPF_TDR_Def_Cat_tpod_table[1])
D_24BPF_TDR_Cat_tpod_upper <-as.numeric(D_24BPF_TDR_Def_Cat_tpod_table[3])
###############################################################################
#LC50 values:
LC50_24BPF_median <- 122
LC50_24BPF_Lower <- 0
LC50_24BPF_upper <- 0

#TDR LC50
LC50_24BPF_TDR_median <- 122
LC50_24BPF_TDR_Lower <- 0
LC50_24BPF_TDR_upper <- 0
##########################################################################
#Table
D_24BPF_table_full <- data.frame(
  Chemical = c("24BPF", "", "", "", "", ""),
  
  Endpoints = c(
    "LCRD",
    "First mode",
    paste0("Category (", D_24BPF_Cat_tpod_path_name, ")"),
    "20th gene",
    "10th percentile",
    "LC50"
  ),
  
  Lower = c(
    LCRD2_24BPF_Tpod_lower,
    First_mode_24BPF_Tpod_lower,
    D_24BPF_Cat_tpod_lower,
    twentieth_gene_lower_bound_24BPF,
    tenth_percentile_lower_bound_24BPF,
    LC50_24BPF_Lower
  ),
  
  Median = c(
    LCRD2_24BPF_Tpod_median,
    First_mode_24BPF_Tpod_median,
    D_24BPF_Cat_tpod_median,
    twenty_gene_24BPF_Tpod_median,
    tenth_percentile_24BPF_Tpod_median,
    LC50_24BPF_median
  ),
  
  Upper = c(
    LCRD2_24BPF_Tpod_upper,
    First_mode_24BPF_Tpod_upper,
    D_24BPF_Cat_tpod_upper,
    twentieth_gene_upper_bound_24BPF,
    tenth_percentile_upper_bound_24BPF,
    LC50_24BPF_upper
  ),
  
  Range = c(
    LCRD2_24BPF_Tpod_upper - LCRD2_24BPF_Tpod_lower,
    First_mode_24BPF_Tpod_upper - First_mode_24BPF_Tpod_lower,
    D_24BPF_Cat_tpod_upper - D_24BPF_Cat_tpod_lower,
    twentieth_gene_upper_bound_24BPF - twentieth_gene_lower_bound_24BPF,
    tenth_percentile_upper_bound_24BPF - tenth_percentile_lower_bound_24BPF,
    LC50_24BPF_upper - LC50_24BPF_Lower
  ),
  
  `DRG (Dose Responsive Genes)` = c(DRG_24BPF, "", "", "", "", ""),
  Top_Dose = c(high_dose_24BPF, "", "", "", "", ""),
  Low_Dose = c(low_dose_24BPF, "", "", "", "", "")
)

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Plots")

write.csv(D_24BPF_table_full, "24BPF_table_full.csv", row.names = FALSE)


#TDR table 
D_24BPF_TDR_table_full <- data.frame(
  Chemical = c("24BPF_TDR", "", "", "", "", ""),
  
  Endpoints = c(
    "LCRD",
    "First mode",
    paste0("Category (", D_24BPF_TDR_Cat_tpod_path_name, ")"),
    "20th gene",
    "10th percentile",
    "LC50"
  ),
  
  Lower = c(
    LCRD2_24BPF_TDR_Tpod_lower,
    First_mode_24BPF_TDR_Tpod_lower,
    D_24BPF_TDR_Cat_tpod_lower,
    twentieth_gene_lower_bound_24BPF_TDR,
    tenth_percentile_lower_bound_24BPF_TDR,
    LC50_24BPF_TDR_Lower
  ),
  
  Median = c(
    LCRD2_24BPF_TDR_Tpod_median,
    First_mode_24BPF_TDR_Tpod_median,
    D_24BPF_TDR_Cat_tpod_median,
    twenty_gene_24BPF_TDR_Tpod_median,
    tenth_percentile_24BPF_TDR_Tpod_median,
    LC50_24BPF_TDR_median
  ),
  
  Upper = c(
    LCRD2_24BPF_TDR_Tpod_upper,
    First_mode_24BPF_TDR_Tpod_upper,
    D_24BPF_TDR_Cat_tpod_upper,
    twentieth_gene_upper_bound_24BPF_TDR,
    tenth_percentile_upper_bound_24BPF_TDR,
    LC50_24BPF_TDR_upper
  ),
  
  Range = c(
    LCRD2_24BPF_TDR_Tpod_upper - LCRD2_24BPF_TDR_Tpod_lower,
    First_mode_24BPF_TDR_Tpod_upper - First_mode_24BPF_TDR_Tpod_lower,
    D_24BPF_TDR_Cat_tpod_upper - D_24BPF_TDR_Cat_tpod_lower,
    twentieth_gene_upper_bound_24BPF_TDR - twentieth_gene_lower_bound_24BPF_TDR,
    tenth_percentile_upper_bound_24BPF_TDR - tenth_percentile_lower_bound_24BPF_TDR,
    LC50_24BPF_TDR_upper - LC50_24BPF_TDR_Lower
  ),
  
  `DRG (Dose Responsive Genes)` = c(DRG_24BPF_TDR, "", "", "", "", ""),
  Top_Dose = c(high_dose_24BPF_TDR, "", "", "", "", ""),
  Low_Dose = c(low_dose_24BPF_TDR, "", "", "", "", ""))

#Full table merged
library(dplyr)

empty_row <- as.data.frame(
  matrix(" ", nrow = 1, ncol = ncol(D_24BPF_table_full))
)

merged_24BPF_FULL_Table_1 <- bind_rows(
  D_24BPF_table_full,
  empty_row,
  D_24BPF_TDR_table_full
)

merged_24BPF_FULL_Table <- merged_24BPF_FULL_Table_1[,1:9]

setwd("C:/Users/KumarA/Downloads/R_Stuff/Actual Analysis/FINAL_tpod/Generic/Plots")

write.csv(merged_24BPF_FULL_Table, "24BPF_table_full_TDR_merged.csv", row.names = FALSE)


##########################################################################
#Density plot generation

library(dplyr)
library(ggplot2)

df <- df_24BPF %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    as.numeric(p10),
    as.numeric(LCRD2_24BPF_Tpod_median),
    as.numeric(First_mode_24BPF_Tpod_median),
    twenty_gene_24BPF_Tpod_median,
    D_24BPF_Cat_tpod_median,
    as.numeric(high_dose_24BPF),
    as.numeric(low_dose_24BPF)
  )),
  
  label = c(
    paste0("10th Percentile = ", round(as.numeric(p10), 3)),
    paste0("LCRD = ", round(LCRD2_24BPF_Tpod_median, 3)),
    paste0("First Mode = ", round(First_mode_24BPF_Tpod_median, 3)),
    paste0("20th Gene = ", round(twenty_gene_24BPF_Tpod_median, 3)),
    paste0("Category = ", round(D_24BPF_Cat_tpod_median, 3)),
    paste0("Top dose = ", as.numeric(high_dose_24BPF)),
    paste0("Low dose = ", as.numeric(low_dose_24BPF))
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
                y = c(0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14),
                label = round(value, 3),
                color = label),
            angle = 90,
            hjust = 0,
            show.legend = FALSE) +
  
  labs(x = "log10(Best.BMD)",
       y = "Density",
       color = "TPoD Methods",
       title = "                                            Distribution of Best BMD Values for 24BPF",  
       subtitle = paste0(
         "                                                  DRG = ", round(DRG_24BPF, 3),
         "   |    Top Dose = ", high_dose_24BPF,
         "   |    Low Dose = ", low_dose_24BPF
       )
  )
###############################################################################

#TDR Plot

###############################################################################
#Density plot generation

library(dplyr)
library(ggplot2)

df <- df_24BPF_TDR %>%
  mutate(logBestBMD = log10(Best.BMD))

# Dataframe for TPoD lines
tpod_df <- data.frame(
  value = log10(c(
    p10,
    LCRD2_24BPF_TDR_Tpod_median,
    First_mode_24BPF_TDR_Tpod_median,
    twenty_gene_24BPF_TDR_Tpod_median,
    D_24BPF_TDR_Cat_tpod_median,
    as.numeric(high_dose_24BPF_TDR),
    as.numeric(low_dose_24BPF_TDR)
  )),
  
  label = c(
    paste0("10th Percentile = ", round(p10, 3)),
    paste0("LCRD = ", round((LCRD2_24BPF_TDR_Tpod_median), 3)),
    paste0("First Mode = ", round((First_mode_24BPF_TDR_Tpod_median), 3)),
    paste0("20th Gene = ", round((twenty_gene_24BPF_TDR_Tpod_median), 3)),
    paste0("Category = ", round((D_24BPF_TDR_Cat_tpod_median), 3)),
    paste0("Top dose = ", as.numeric(high_dose_24BPF_TDR)),
    paste0("Low dose = ", as.numeric(low_dose_24BPF_TDR))
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
       title = "                                            Distribution of Best BMD Values for 24BPF_TDR",  
       subtitle = paste0(
         "                                                  DRG = ", round(DRG_24BPF_TDR, 3),
         "   |    Top Dose = ", high_dose_24BPF_TDR,
         "   |    Low Dose = ", low_dose_24BPF_TDR
       )
  )



