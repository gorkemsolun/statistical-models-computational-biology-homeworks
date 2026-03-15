##########################################################################
###########             Project 4:  PROFILE HMM                 ##########
##########################################################################

#
# Parse an alignment from a file
#
parseAlignment <- function(alignmentFile) {
  alignment <- scan(file = alignmentFile, what = character(0))
  # Split rows and convert to matrix
  alignment.mat <- matrix(nrow = length(alignment), ncol = nchar(alignment[[1]]))
  for (i in 1:length(alignment)) {
    alignment.mat[i, ] <- unlist(strsplit(alignment[i], ""))
  }
  return(alignment.mat)
}
getwd()

GTP_binding_proteins <- parseAlignment("../data/GTP_binding_proteins.txt")
ATPases <- parseAlignment("../data/ATPases.txt")

View(ATPases)
#
# Learn HMM from data
#
learnHMM <- function(alignment, alphabet = NULL, counts.only = FALSE) {
  # Lists with all amino acids as symbols and the state transitions
  transitions <- c("MM", "MD", "MI", "IM", "ID", "II", "DM", "DD", "DI")
  if (is.null(alphabet)) {
    alphabet <- c(
      "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
      "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
    )
  }
  alRows <- nrow(alignment)
  alCols <- ncol(alignment)
  L <- 0 # Length of the model
  stateNum <- array(0, dim = alCols) # Keeping track of state number
  assignedPos <- array(FALSE, dim = alCols)
  # Getting the model states and length
  for (i in 1:alCols) {
    hyphenfreq <- sum(alignment[, i] == "-") # Hyphen frequencies
    if (hyphenfreq < alRows / 2) {
      L <- L + 1
      assignedPos[i] <- TRUE
    }
    stateNum[i] <- L
  }
  # Initialise matrices for: (index shift because of begin state)
  # match emissions (mE)
  # insertion emssions (iE)
  # transitions (T)
  mE <- matrix(0, nrow = length(alphabet), ncol = L + 1, dimnames = list(alphabet, paste(0:L, sep = "")))
  iE <- matrix(0, nrow = length(alphabet), ncol = L + 1, dimnames = list(alphabet, NULL))
  T <- matrix(0, nrow = length(transitions), ncol = L + 1, dimnames = list(transitions, NULL))

  mE[, 1] <- NA

  for (i in 1:alRows) {
    prev <- "M" # Assume starting with a match
    prevStateNum <- 0 # Assume a begin state
    for (j in 1:alCols) {
      aa <- alignment[i, j]
      if (assignedPos[j] == TRUE) {
        if (aa == "-") {
          newState <- "D"
        } else {
          newState <- "M"
          mE[aa, stateNum[j] + 1] <- mE[aa, stateNum[j] + 1] + 1
        }
      } else {
        if (aa == "-") {
          next
        } else {
          iE[aa, stateNum[j] + 1] <- iE[aa, stateNum[j] + 1] + 1
          newState <- "I"
        }
      }
      # Update transition matrix
      # j -> j+1 save at entry j
      transition <- paste(prev, newState, sep = "")
      T[transition, prevStateNum + 1] <- T[transition, prevStateNum + 1] + 1
      prevStateNum <- stateNum[j]
      prev <- newState
    }

    # Last state is always match
    transition <- paste(prev, "M", sep = "")
    T[transition, prevStateNum + 1] <- T[transition, prevStateNum + 1] + 1
  }

  # If !counts.only, add a pseudocount of 1
  if (!counts.only) {
    # Transition matrix
    T <- T + 1
    for (i in 1:ncol(T)) {
      T[1:3, i] <- T[1:3, i] / sum(T[1:3, i])
      T[4:6, i] <- T[4:6, i] / sum(T[4:6, i])
      T[7:9, i] <- T[7:9, i] / sum(T[7:9, i])
    }
    # Match emission matrix
    mE <- mE + 1
    mE <- apply(mE, 2, function(x) x / sum(x))
    # Insertion emission matrix
    iE <- iE + 1
    iE <- apply(iE, 2, function(x) x / sum(x))
  }
  return(list(
    T = T, mE = mE, iE = iE, L = L, alphabet = alphabet, transitions = transitions,
    stateNum = stateNum
  ))
}

ATPase_profileHMM <- learnHMM(ATPases)
GTPbinding_profileHMM <- learnHMM(GTP_binding_proteins)

# Helper function
analyzeHMM <- function(hmm, family_name) {
  mE <- hmm$mE
  iE <- hmm$iE

  max_match_per_pos <- apply(mE, 2, max, na.rm = TRUE)
  best_match_pos <- which(max_match_per_pos == max(max_match_per_pos))

  max_insert_per_pos <- apply(iE, 2, max)
  best_insert_pos <- which(max_insert_per_pos == max(max_insert_per_pos))

  cat("\n============================\n")
  cat("Family:", family_name, "\n")
  cat("============================\n")

  cat("Best match position(s):", paste(best_match_pos - 1, collapse = ", "), "\n")
  for (pos in best_match_pos) {
    sym <- rownames(mE)[which(mE[, pos] == max(mE[, pos], na.rm = TRUE))]
    freq <- max(mE[, pos], na.rm = TRUE)

    cat("M", pos - 1, ": ", paste(sym, collapse = ", "),
      " with frequency ", round(freq, 6), "\n",
      sep = ""
    )

    cat("Full match emission distribution:\n")
    print(round(mE[, pos], 6))
  }

  cat("Best insert position(s):", paste(best_insert_pos - 1, collapse = ", "), "\n")
  for (pos in best_insert_pos) {
    sym <- rownames(iE)[which(iE[, pos] == max(iE[, pos]))]
    freq <- max(iE[, pos])

    cat("I", pos - 1, ": ", paste(sym, collapse = ", "),
      " with frequency ", round(freq, 6), "\n",
      sep = ""
    )

    cat("Full insert emission distribution:\n")
    print(round(iE[, pos], 6))
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  n_plots <- max(length(best_match_pos), length(best_insert_pos))
  par(mfrow = c(n_plots, 2), mar = c(5, 4, 4, 1))

  for (k in seq_len(n_plots)) {
    if (k <= length(best_match_pos)) {
      pos <- best_match_pos[k]

      bp <- barplot(mE[, pos],
        ylim = c(0, 1),
        col = "steelblue",
        xaxt = "n",
        main = paste(family_name, "- Match emissions at M", pos - 1),
        ylab = "Emission frequency",
        xlab = "Amino acid"
      )

      axis(1, at = bp, labels = rownames(mE), las = 2, cex.axis = 0.8)
    } else {
      plot.new()
    }

    if (k <= length(best_insert_pos)) {
      pos <- best_insert_pos[k]

      bp <- barplot(iE[, pos],
        ylim = c(0, 1),
        col = "darkorange",
        xaxt = "n",
        main = paste(family_name, "- Insert emissions at I", pos - 1),
        ylab = "Emission frequency",
        xlab = "Amino acid"
      )

      axis(1, at = bp, labels = rownames(iE), las = 2, cex.axis = 0.8)
    } else {
      plot.new()
    }
  }
}

# Run for both families
analyzeHMM(ATPase_profileHMM, "ATPase")
analyzeHMM(GTPbinding_profileHMM, "GTP binding")


#
# parse file with one protein per line
#
parseProteins <- function(proteinsFile) {
  proteins <- scan(file = proteinsFile, what = character(0))
  proteinList <- as.list(strsplit(proteins, ""))
  return(proteinList)
}

unclassifiedProteins <- parseProteins("../data/Unclassified_proteins.txt")

#
# The Forward algorithm
#
forward <- function(HMM, seq) {
  L <- HMM$L
  mE <- HMM$mE
  iE <- HMM$iE
  T <- HMM$T
  len <- length(seq)

  # Random model: same emission probability for all letters in the alphabet
  qr <- 1 / nrow(mE)
  # Initialize score matrices for match, insert and delete
  Fm <- Fi <- Fd <- matrix(-Inf, nrow = len + 1, ncol = L + 1)
  # First column for the begin state (M0) (indices get shifted by 1)
  Fm[1, 1] <- 0

  # We then allow transitions to I0 and D1 (start with insertion or deletion)
  for (i in 2:len + 1) {
    Fi[i, 1] <- log(iE[seq[i - 1], 1] / qr) +
      log(T["MI", 1] * exp(Fm[i - 1, 1]) +
        T["II", 1] * exp(Fi[i - 1, 1]))
  }

  # paths starting with deletions should also be available
  for (j in 2:(L + 1)) {
    Fd[1, j] <- log(T["MD", j - 1] * exp(Fm[1, j - 1]) +
      T["ID", j - 1] * exp(Fi[1, j - 1]) +
      T["DD", j - 1] * exp(Fd[1, j - 1]))
  }

  # i is the pos in the sequence, j is the profile state
  for (i in 2:(len + 1)) {
    for (j in 2:(L + 1)) {
      # Match state M_{j-1}
      Fm[i, j] <- log(mE[seq[i - 1], j] / qr) +
        log(T["MM", j - 1] * exp(Fm[i - 1, j - 1]) +
          T["IM", j - 1] * exp(Fi[i - 1, j - 1]) +
          T["DM", j - 1] * exp(Fd[i - 1, j - 1]))

      # Insert state I_{j-1}
      Fi[i, j] <- log(iE[seq[i - 1], j] / qr) +
        log(T["MI", j] * exp(Fm[i - 1, j]) +
          T["II", j] * exp(Fi[i - 1, j]) +
          T["DI", j] * exp(Fd[i - 1, j]))

      # Delete state D_{j-1}
      Fd[i, j] <- log(T["MD", j - 1] * exp(Fm[i, j - 1]) +
        T["ID", j - 1] * exp(Fi[i, j - 1]) +
        T["DD", j - 1] * exp(Fd[i, j - 1]))
    }
  }

  # Termination step: Fe = Fm[len,L+2] = log(P(x|M)/P(x|R))
  Fe <- log(T["MM", L + 1] * exp(Fm[len + 1, L + 1]) +
    T["IM", L + 1] * exp(Fi[len + 1, L + 1]) +
    T["DM", L + 1] * exp(Fd[len + 1, L + 1]))
  return(Fe)
}


# M1 = GTP binding, M2 = ATPase
M1 <- GTPbinding_profileHMM
M2 <- ATPase_profileHMM

classifyProteins <- function(protein_list, M1, M2,
                             name_M1 = "GTP binding",
                             name_M2 = "ATPase") {
  n <- length(protein_list)
  scores_M1 <- numeric(n)
  scores_M2 <- numeric(n)
  q <- numeric(n)
  predicted_family <- character(n)

  for (i in seq_len(n)) {
    x <- protein_list[[i]]

    scores_M1[i] <- forward(M1, x)
    scores_M2[i] <- forward(M2, x)
    q[i] <- scores_M1[i] - scores_M2[i]

    if (q[i] > 0) {
      predicted_family[i] <- name_M1
    } else if (q[i] < 0) {
      predicted_family[i] <- name_M2
    } else {
      predicted_family[i] <- "Undecided"
    }
  }

  results <- data.frame(
    protein_index = seq_len(n),
    score_M1 = scores_M1,
    score_M2 = scores_M2,
    q = q,
    predicted_family = predicted_family
  )

  print(results)

  barplot(q,
    names.arg = seq_len(n),
    col = ifelse(q > 0, "steelblue", "darkorange"),
    main = expression(q(x[i]) == log(frac(P(x[i] ~ "|" ~ M[1]), P(x[i] ~ "|" ~ M[2])))),
    xlab = "Protein index",
    ylab = "Log odds ratio"
  )
  abline(h = 0, lwd = 2)

  return(results)
}

results <- classifyProteins(unclassifiedProteins, M1, M2)
