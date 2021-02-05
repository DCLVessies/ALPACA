ALPACA <- function(FPR=FalsePositiveRates) {
  require(dplyr)
  
  # Helper function to set the LoB based on nearest FPR
  ALPACA.adaptiveLoB <- function(FPRtable=FalsePositiveRates, assay, concentration, sample, nr.molecules) {
    a <- FPRtable[FPRtable$Assay == assay,]
    if (nrow(a) == 0) {LoB <- NA}
    
    else {
      FPR <- a$FPR[abs(a$cpul-concentration) == min(abs(a$cpul-concentration))]
      if (length(FPR) == 2) {
        print(paste0("Notification: ", sample, " has averaged FPR"))
        FPR <- mean(FPR)
      }
      LoB <- which.max(sapply(X=1:50, function(FPR, nr.molecules, cutoff) sum(dbinom(x=0:cutoff-1, size=nr.molecules, prob=FPR)), FPR=FPR, nr.molecules=nr.molecules) >= 0.999)
    }
    LoB
  }
  
  # Set the correct working directory
  wd <- readline(prompt="Please enter the path to the directory that contains your ddPCR .csv file(s): ")
  setwd(wd)
  
  # Select the correct CSV files
  csv_list <- select.list(list.files(pattern="*.csv"), multiple=TRUE, graphics=FALSE)

  dff <- lapply(X=csv_list, FUN=function(csvs) read.csv(paste0(getwd(), "/", csvs), row.names=NULL, stringsAsFactors=FALSE)[,c(4:6,17:20)])
  df <- bind_rows(dff)
  df <- df[df[,2]=="Ch1Unknown",-2]
  colnames(df) <- c("Sample", "Assay", "Ch1Ch2", "Ch1", "Ch2", "Empty")

  # Merge data by sample and assay combination
  df$Sample <- paste(df$Sample, df$Assay)
  df$nwell <- rep(1, nrow(df))
  dfmerge <- aggregate(cbind(df$nwell, df$Ch1Ch2, df$Ch1, df$Ch2, df$Empty), by=list(Category=df$Sample), FUN=sum)
  colnames(dfmerge) <- c("Sample", "nr.Wells", "Ch1Ch2", "Ch1", "Ch2", "Empty")
  dfmerge$N <- rowSums(dfmerge[,3:6])
  dfmerge$Assay <- df$Assay[match(dfmerge$Sample, df$Sample)]
  
  # Check No Template Control (NTC) results
  NTC <- dfmerge[grep("!NTC", dfmerge$Sample, ignore.case=TRUE), ]
  
  if (nrow(NTC) != length(unique(dfmerge$Assay))) {
    assays <- unique(dfmerge$Assay)
    missing <- assays[!assays %in% NTC$Assay]
    dfmerge <- dfmerge[!dfmerge$Assay %in% missing,]
    print(paste0("WARNING! No Template Control (NTC) results for assay(s) ", paste(missing, collapse=" and "), " not found. SAMPLES NOT ANALYSED! -- Please give NTC wells the name '!NTC' in QuantaSoft so that ALPACA can recognise them."))
  }
  
  # Perform the PIF correction
  dfmerge$nPIFs <- mapply(FUN=function(Ch1, Ch2, Ch1Ch2, Empty, N) which.max(sapply(X=Ch2:(Ch2+Ch1Ch2), FUN=function(Ch1, Ch2, Empty, N) sum(dbinom(x=(N-sum(Ch1, Ch2, Empty)):(N-sum(Ch2, Empty)), size=(N-sum(Ch2, Empty)), prob=Ch2/(Ch2+Empty))), Ch1=Ch1, Empty=Empty, N=N) >= 0.1)-1, Ch1=dfmerge$Ch1, Ch2=dfmerge$Ch2, Ch1Ch2=dfmerge$Ch1Ch2, Empty=dfmerge$Empty, N=dfmerge$N)
  
  # Calculate number of mutant and wildtype molecules based on PIF corrected droplet counts (Poisson)
  dfmerge$nr.Mutant <- round(-log((dfmerge$Ch2+dfmerge$Empty+dfmerge$nPIFs)/dfmerge$N)*dfmerge$N)
  dfmerge$nr.Wildtype <- round(-log((dfmerge$Ch1+dfmerge$Empty)/dfmerge$N)*dfmerge$N)
  dfmerge$Concentration.cpul <- round((dfmerge$nr.Mutant+dfmerge$nr.Wildtype)/dfmerge$N/0.85*1000)
  dfmerge$VAF <- round(dfmerge$nr.Mutant/rowSums(cbind(dfmerge$nr.Mutant, dfmerge$nr.Wildtype)), digits=5)
  dfmerge$VAF[is.na(dfmerge$VAF)] <- 0

  # Apply dynamic LoB based on FPR
  dfmerge$LoB <- unlist(mapply(FUN=ALPACA.adaptiveLoB, assay=dfmerge$Assay, concentration=dfmerge$Concentration.cpul, sample=dfmerge$Sample, nr.molecules=dfmerge$nr.Mutant+dfmerge$nr.Wildtype))
  
  # If some assay names used in QuantaSoft are not in the FalsePositiveRates list, give a warning
  noLoB <- unique(dfmerge$Assay[is.na(dfmerge$LoB)])
  if (length(noLoB) >0) {
    print(paste0("ERROR: Assay name(s) ", paste(noLoB, collapse=", "), " not recognised. LoB not applied!"))
    print(paste0("Accepted assay names are: ", paste(unique(FPR$Assay), collapse=", ")))
  }
  dfmerge$MutationCalled <- (dfmerge$Ch1 + dfmerge$Ch1Ch2 - dfmerge$nPIFs) >= dfmerge$LoB
  dfmerge$VAF[!dfmerge$MutationCalled] <- 0
  
  RESULT <- dfmerge[order(dfmerge$Assay), c(1:6, 9, 10, 14, 15, 12, 13)]
  RESULT
}

FalsePositiveRates <- data.frame(matrix(c(
  "1","BRAF V600E",19,6e-05,
  "2","BRAF V600E",297,4e-05,
  "3","BRAF V600E",472,5e-06,
  "4","EGFR E746-A750del",34,8e-05,
  "5","EGFR E746-A750del",210,6e-05,
  "6","EGFR E746-A750del",270,5e-06,
  "7","EGFR exon19 del",33,1e-04,
  "8","EGFR exon19 del",77,3e-05,
  "9","EGFR exon19 del",235,4e-05,
  "10","EGFR L858R",34,4e-05,
  "11","EGFR L858R",173,2e-05,
  "12","EGFR L858R",309,4e-06,
  "13","EGFR L858R",441,3e-6,
  "14","EGFR T790M",5,3e-04,
  "15","EGFR T790M",34,3e-05,
  "16","EGFR T790M",144,9e-05,
  "17","EGFR T790M",471,4e-05,
  "18","KRAS G12/G13",25,3e-04,
  "19","KRAS G12/G13",250,4e-05
), ncol=4, byrow=TRUE), stringsAsFactors = FALSE)

colnames(FalsePositiveRates) <- c(" ", "Assay","cpul","FPR")
rownames(FalsePositiveRates) <- FalsePositiveRates[,1]
FalsePositiveRates <- FalsePositiveRates[,-1]
FalsePositiveRates$Assay <- as.character(FalsePositiveRates$Assay)
FalsePositiveRates$cpul <- as.numeric(FalsePositiveRates$cpul)
FalsePositiveRates$FPR <- as.numeric(FalsePositiveRates$FPR)
