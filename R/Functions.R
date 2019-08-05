
# find density gap
findGapSmooth <- function(vec) {
    # bandwidth <- kedd::h.ucv(vec, 4)$h
    # bandwidth <- density(vec)$bw
    # dens <- KernSmooth::bkde(vec, canonical=TRUE) #, bandwidth=bandwidth
    suppressWarnings(dens <- density(vec, kernel="gaussian", adjust=1, bw="ucv"))  #"nrd0"
    tryCatch({
        idy <- which(diff(sign(diff(dens$y))) == 2) + 1
        idx <- idy[which.min(dens$y[idy])]
        d <- ifelse(length(idx) != 0, dens$x[idx], NA)
    },
    error = function(e) {
        warning('No minimum was found between two modes.')
        return(NULL) })
    return(d)
}
# *****************************************************************************

# *****************************************************************************
# epsilon calculator by "infer"
infer <- function(vec) {
    vec <- sort(vec)
    # vec[1] <- vec[2]/2
    n <- length(vec)
    d <- NA
    # upper level search
    if (n > 2) {
        d <- findGapSmooth(vec=vec)    
        if (!is.na(d)) { 
            # d <- max(vec[vec <= d]) 
            d <- ifelse(max(vec[vec <= d]) == 0, d, max(vec[vec <= d]))
        }
    }
    # lower level search
    if (is.na(d)) {
        diffVec <- diff(vec)
        if (length(unique(diffVec[diffVec > 0])) == 1) {
            d <- ceiling(mean(vec))
        } else {
            x <- which.max(diffVec)
            d <- ifelse(vec[x] == 0, mean(c(vec[x], vec[vec>0][1])), vec[x])
        }   
    }
    return(d)
}
# *****************************************************************************

# *****************************************************************************
# kernel matrix calculator
krnlMtxGenerator <- function(mtx) {
    # Radial basis function kernel: In the Gaussian Kernel if two points are
    # close then K_ij≈1 and when two points are far apart then Kij≈0
    n <- nrow(mtx)
    # calculate epsilons
    epsilon <- rep(0, length=n)
    for (i in 1:n) {
        epsilon[i] <- infer(vec=mtx[i,])
    }    
    # calculate kernel matrix
    krnl_mtx <- matrix(data=1, nrow=n, ncol=n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            krnl_mtx[i,j] <- exp(-mtx[i,j]^2/(epsilon[i]*epsilon[j]))  #2*
            krnl_mtx[j,i] <- krnl_mtx[i,j]
        }
    }
    krnl_mtx[is.nan(krnl_mtx)] <- 1  # if mtx[i,j] and epsilon == 0
    return(krnl_mtx)
}
# *****************************************************************************

# *****************************************************************************
# affinity matrix calculator
# Disconnect those edges with distance larger than threshold
makeAffinity <- function(mtx_o, mtx_k, thd) {
    mtx_k[mtx_o > thd] <- 0
    return(mtx_k)
}
# *****************************************************************************

# *****************************************************************************
# laplacian matrix calculator
laplacianMtx <- function(entry) {
    # Calculate unnormalised Laplacian matrix and its eigenfunctions
    D <- diag(apply(entry, 1, sum))
    L <- D - entry
    return(L)
}
# *****************************************************************************

# *****************************************************************************
# range a vector from a to b
rangeAtoB <- function(x, a, b){
    return((b-a)*(x-min(x))/(max(x)-min(x)) + a)
}
# *****************************************************************************

# *****************************************************************************
# calculate likelihoods
likelihoods <- function(tot_mtx, sh_mtx, mutab_mtx) {
    sigma_sh <- sd(sh_mtx[upper.tri(sh_mtx)])
    sigma_tot <- sd(tot_mtx[upper.tri(tot_mtx)])
    if (sigma_sh %in% c(NA, 0)) { # there is no any prefences among pairs of sequences, therefore all likelihoods are zero
        z_mtx <- matrix(0, nrow=nrow(sh_mtx), ncol=nrow(sh_mtx))
    } else if (sigma_tot %in% c(NA, 0)) { # pairs with more shared mutations are more likly to belong to the same clone
        x_mtx <- 1.0 - exp(-sh_mtx^2/(2.0*sigma_sh^2))
        z_mtx <- mutab_mtx*x_mtx
    } else {
        x_mtx <- 1.0 - exp(-sh_mtx^2/(2.0*sigma_sh^2))
        y_mtx <- exp(-(tot_mtx-sh_mtx)^2/(2.0*sigma_tot^2))
        z_mtx <- mutab_mtx*x_mtx*y_mtx
    }
    diag(z_mtx) <- 1
    return(z_mtx)
}
# *****************************************************************************

# *****************************************************************************
pairwiseMutions <- function(germ_imgt, 
                            seq_imgt, 
                            junc_length, 
                            len_limit = NULL,
                            cdr3 = FALSE,
                            mutabs = NULL) {
    
    ##### get number of seqs
    n <- unique(c(length(seq_imgt), length(germ_imgt)))
    ##### check number of sequences
    if (length(n) > 1) stop("germ_imgt and seq_imgt number should be the same")
    if (n == 1) stop("there should be at least two seqs")
    # check consensus length
    if (!is.null(len_limit)) {
        lenConsensus <- len_limit@seqLength
        seq_imgt  <- substr(seq_imgt, start = 1, stop = lenConsensus)
        germ_imgt <- substr(germ_imgt, start = 1, stop = lenConsensus)
        eff_germ <- ifelse(length(unique(germ_imgt)) == 1, 
                           unique(germ_imgt), 
                           consensusSequence(sequences = unique(germ_imgt),
                                             muFreqColumn = NULL, 
                                             lenLimit = lenConsensus,
                                             method = "catchAll",
                                             minFreq = NULL,
                                             includeAmbiguous = FALSE,
                                             breakTiesStochastic = FALSE,
                                             breakTiesByColumns = NULL, 
                                             db = NULL)$cons)
    } else {
        ##### constants
        lv <- ifelse(cdr3, shazam::IMGT_V@seqLength, shazam::IMGT_V@seqLength - 3)
        trim_l <- junc_length
        ##### trim out junction/cdr3 segments from seq_imgt
        seq_imgt <- sapply(1:length(seq_imgt), function(i){
            x <- s2c(seq_imgt[i])         # x <- strsplit(seq_imgt[i], split="")[[1]]
            x[(lv+1):(lv+trim_l)] <- ""   # x[(lv+1):(lv+trim_l[i])] <- ""
            return(paste(x, collapse=""))
        })
        ##### Pads ragged ends
        l <- unique(nchar(seq_imgt))
        if (length(l) > 1) {
            seq_imgt <- padSeqEnds(seq = seq_imgt, len = NULL, start = FALSE, pad_char = "N")
        }
        ##### trim out junction/cdr3 segments from germ_imgt
        germ_imgt <- sapply(1:length(germ_imgt), function(i){
            x <- s2c(germ_imgt[i])       # x <- strsplit(germ_imgt[i], split="")[[1]]
            x[(lv+1):(lv+trim_l)] <- ""  # x[(lv+1):(lv+trim_l[i])] <- ""
            return(paste(x, collapse=""))
        })
        ##### Pads ragged ends
        l <- unique(nchar(germ_imgt))
        if (length(l) > 1) {
            germ_imgt <- padSeqEnds(seq = germ_imgt, len = NULL, start = FALSE, pad_char = "N")
        }
        ##### find consensus germline (allel level grouping)
        # see arg "method" from shazam::collapseClones function
        eff_germ <- ifelse(length(unique(germ_imgt)) == 1, 
                           unique(germ_imgt), 
                           consensusSequence(sequences = unique(germ_imgt),
                                             muFreqColumn = NULL, 
                                             lenLimit = NULL,
                                             method = "catchAll",
                                             minFreq = NULL,
                                             includeAmbiguous = FALSE,
                                             breakTiesStochastic = FALSE,
                                             breakTiesByColumns = NULL, 
                                             db = NULL)$cons)
        ##### check germ and seqs lengths
        seq_imgt_lent <- unique(nchar(seq_imgt))
        germ_imgt_lent <- unique(nchar(germ_imgt))
        eff_germ_lent <- nchar(eff_germ)
        lenConsensus <- min(seq_imgt_lent, germ_imgt_lent, eff_germ_lent) 
        ##### trim extra characters
        if ( seq_imgt_lent > lenConsensus)  { seq_imgt <- substr( seq_imgt, start = 1, stop = lenConsensus) }
        if (germ_imgt_lent > lenConsensus) { germ_imgt <- substr(germ_imgt, start = 1, stop = lenConsensus) }
        if ( eff_germ_lent > lenConsensus)  { eff_germ <- substr( eff_germ, start = 1, stop = lenConsensus) }  
    }
    ##### count informative positions
    informative_pos <- sapply(1:n, function(x){ sum(stri_count(seq_imgt[x], fixed = c("A","C","G","T"))) })
    ##### convert eff_germ and seq_imgt to matrices
    effMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
    seqsMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
    for (i in 1:n) {
        seqsMtx[i, ] <- s2c(seq_imgt[i])[1:lenConsensus]
        effMtx[i, ] <- s2c(eff_germ)[1:lenConsensus]
    }
    ##### make a distance matrix
    dnaMtx <- getDNAMatrix(gap = 0)
    mutMtx <- matrix(NA, nrow=n, ncol=lenConsensus)
    for (i in 1:n) {
        mutMtx[i, ] <- sapply(1:lenConsensus, function(j) {
            return(dnaMtx[effMtx[i,j], seqsMtx[i,j]]) 
        })
    }
    ##### make a mutation matrix
    mutMtx <- matrix(paste0(effMtx, mutMtx, seqsMtx), nrow=n, ncol=lenConsensus)
    ##### clean non-mutated elements
    mutMtx[grepl(pattern="0", mutMtx)] <- NA
    ##### check mutabilities
    ##### make a motif matrix
    motifMtx <- matrix(0, nrow=n, ncol=lenConsensus)
    if (!is.null(mutabs)) {
        for (i in 1:n) {
            for (j in 3:(lenConsensus-2)) {
                motifMtx[i, j] <- mutabs[substr(germ_imgt[i], start = j-2, stop = j+2)]
            }
        }
        motifMtx[is.na(motifMtx)] <- 0
    }
    ##### calculate mutation matrix
    results <- pairwiseMutMatrix(informative_pos = informative_pos, 
                                 mutMtx = mutMtx, 
                                 motifMtx = motifMtx)
    sh_mtx <- results$sh_mtx
    tot_mtx <- results$tot_mtx
    mutab_mtx <- results$mutab_mtx
    ##### make symmetric matrix
    sh_mtx[lower.tri(sh_mtx)] <- t(sh_mtx)[lower.tri(sh_mtx)]
    tot_mtx[lower.tri(tot_mtx)] <- t(tot_mtx)[lower.tri(tot_mtx)]
    mutab_mtx[lower.tri(mutab_mtx)] <- t(mutab_mtx)[lower.tri(mutab_mtx)]
    # return results
    return_list <- list("pairWiseSharedMut" = sh_mtx,
                        "pairWiseTotalMut" = tot_mtx,
                        "pairWiseMutability" = mutab_mtx)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
# inter-clone-distance vs intra-clone-distance
calculateInterVsIntra <- function(db,
                                  clone_col,
                                  vjl_groups,
                                  junction_col = "JUNCTION",
                                  cdr3 = FALSE,
                                  cdr3_col = NA,
                                  nproc = 1,
                                  verbose = FALSE) {
    ### Create cluster of nproc size and export namespaces
    if(nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK", outfile = "")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    n_groups <- nrow(vjl_groups)  
    ### check the progressbar
    if (verbose) {
        pb <- progressBar(n_groups)
    }
    
    k <- NULL
    # open dataframes
    vec_ff <- foreach(k=1:n_groups,
                      .combine="c",
                      .errorhandling='stop') %dopar% {
                          
                          # *********************************************************************************
                          clones <- strsplit(vjl_groups$CLONE_ID[k], split=",")[[1]]
                          l <- vjl_groups$JUNCTION_LENGTH[k]
                          n_clones <- length(clones)
                          seqs <- db[[ifelse(cdr3, cdr3_col, junction_col)]][db[[clone_col]] %in% clones]
                          names(seqs) <- db[[clone_col]][db[[clone_col]] %in% clones]
                          ### make a dataframe of unique seqs in each clone
                          seqs_db <- data_frame(value = seqs, name = names(seqs)) %>%
                              dplyr::group_by(!!!rlang::syms(c("name", "value"))) %>% # alternatively: group_by(name) if name value pair is always unique
                              slice(1) %>%
                              ungroup()
                          seqs <- seqs_db$value
                          names(seqs) <- seqs_db$name
                          ### calculate distance matrix among all seqs
                          dist_mtx <- pairwiseDist(seqs, dist_mat=getDNAMatrix(gap = 0))
                          ### prealoocate a vector = no. of max-dist in each clone (intra) + no. of min-dist between clones (inter)
                          nrow_f <- n_clones + n_clones*(n_clones-1)/2
                          vec_f <- rep(NA, nrow_f)
                          ### calculate minimum and maximum distance in each clone
                          n <- 0
                          if (n_clones == 1) {
                              n <- n + 1
                              vec_f[n] <- max(dist_mtx)/l
                              names(vec_f)[n] <- paste(clones[1], "NA", "intra", sep="_")
                          } else {
                              for (i in 1:(n_clones-1)) {
                                  xx <- dist_mtx[rownames(dist_mtx) == clones[i], colnames(dist_mtx) == clones[i]]
                                  n <- n + 1
                                  vec_f[n] <- max(xx)/l
                                  names(vec_f)[n] <- paste(clones[i], "NA", "intra", sep="_")
                                  for (j in (i+1):n_clones) {
                                      xy <- dist_mtx[rownames(dist_mtx) == clones[i], colnames(dist_mtx) == clones[j]]
                                      n <- n + 1
                                      vec_f[n] <- min(xy)/l
                                      names(vec_f)[n] <- paste(clones[i], clones[j], "inter", sep="_")
                                  }
                              }
                              yy <- dist_mtx[rownames(dist_mtx) == clones[j], colnames(dist_mtx) == clones[j]]
                              n <- n + 1
                              vec_f[n] <- max(yy)/l
                              names(vec_f)[n] <- paste(clones[j], "NA", "intra", sep="_")
                          }
                          
                          # Update progress
                          if (verbose) { pb$tick() }
                          
                          return(vec_f)
                          # *********************************************************************************
                      }
    
    ### Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    # convert to a data.frame
    db_dff <- data.frame(keyName = names(vec_ff), 
                         DISTANCE = vec_ff, 
                         row.names=NULL)
    db_dff$LABEL <- "intra"
    db_dff$LABEL[grepl("inter", db_dff$keyName)] <- "inter"
    clones_xy <- data.frame(matrix(unlist(stri_split_fixed(db_dff$keyName, "_", n=3)),
                                   nrow=nrow(db_dff),
                                   byrow=T),
                            stringsAsFactors=FALSE)
    db_dff <- cbind(clones_xy, db_dff)
    db_dff$keyName <- NULL
    colnames(db_dff)[colnames(db_dff) == "X1"] <- "CLONE_X"
    colnames(db_dff)[colnames(db_dff) == "X2"] <- "CLONE_Y"
    db_dff$X3 <- NULL
    
    # return results list
    return_list <- list("inter_intra" = db_dff)
    
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
### plot inter-clone-distance vs intra-clone-distance
plotInterVsIntra <- function(data) {
    
    data <- select(data, c("DISTANCE", "LABEL"))
    data_intra <- data %>% 
        dplyr::filter(!!rlang::sym("LABEL") == "intra", !!rlang::sym("DISTANCE") > 0)
    data_inter <- data %>%
        dplyr::filter(!!rlang::sym("LABEL") == "inter", !!rlang::sym("DISTANCE") > 0)
    
    # fill color
    fill_manual <- c("intra"="grey30",
                     "inter"="grey60")
    
    hline_size <- 0.75
    
    ### find effective threshold
    eff_threshold <- NA
    if (nrow(data_intra) > 0 & nrow(data_inter) > 0) {
        a <- data_intra$DISTANCE
        b <- data_inter$DISTANCE
        xlim = c(min(c(a,b)), max(c(a,b)))
        df <- merge(
            as.data.frame(density(a, from = xlim[1], to = xlim[2])[c("x", "y")]),
            as.data.frame(density(b, from = xlim[1], to = xlim[2])[c("x", "y")]),
            by = "x", suffixes = c(".a", ".b")
        )
        df$comp <- as.numeric(df$y.a > df$y.b)
        df$cross <- c(NA, diff(df$comp))
        df <- df[which(df$cross != 0), c("x", "y.a")]
        if (nrow(df) > 0) {
            eff_th <- df$x
            eff_th <- eff_th[mean(a) - sd(a) < eff_th & eff_th < mean(b) + sd(b)]
            if (length(eff_th) > 0) {
                eff_threshold <- round(mean(eff_th), 2)
            } 
        }    
    }
    
    ### plot p1
    p1 <- ggplot() +
        baseTheme() +
        theme(plot.title = element_text(size = 13, hjust = 0.5),
              axis.title.x = element_blank(),
              axis.title = element_text(size=13),
              axis.text.x=element_text(size=12),
              axis.text.y=element_text(size=12),
              legend.position = "bottom",
              legend.text=element_text(size=12)) +
        xlab("Normalized hamming distance") +
        ylab("Density") +
        scale_y_continuous(labels = abs) +
        scale_fill_manual(name="",
                          values=fill_manual,
                          labels=c("intra"="maximum-distance within clones  ",
                                   "inter"="minimum-distance between clones  "))
    
    if (nrow(data_intra) > 0) {
        p1 <- p1 + 
            geom_histogram(data = data_intra,
                           aes_string(x = "DISTANCE", y = "..density..", fill = "LABEL"),
                           binwidth = 0.02, color = "white", alpha = 0.85) +
            geom_density(data = data_intra, 
                         aes_string(x = "DISTANCE"),
                         size = hline_size, color = "grey30")
    }
    
    if (nrow(data_inter) > 0) {
        p1 <- p1 + 
            geom_histogram(data = data_inter,
                           aes_string(x = "DISTANCE", y = "..density..", fill = "LABEL"),
                           binwidth=0.02, color = "white", alpha = 0.75) +
            geom_density(data = data_inter, 
                         aes_string(x = "DISTANCE"),
                         size = hline_size, color = "grey60")
    } else {
        warning("No inter clonal distance is detected. Each group of sequences with same V-gene, J-gene, and junction length may contain only one clone.")
    }
    
    if (!is.na(eff_threshold)) {
        p1 <- p1 + 
            ggtitle(paste("Effective threshold=", eff_threshold)) +
            geom_vline(xintercept=eff_threshold, color="grey30", linetype=2, size=hline_size)
    } else {
        p1 <- p1 + 
            ggtitle(paste("Effective threshold not found"))
    }
    
    # return results list
    return_list <- list("eff_threshold" = eff_threshold,
                        "plot_inter_intra" = p1)
    
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
### Define verbose reporting function
printVerbose <- function(n_groups, vjl_group, model, method, cdr3,
                         gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) {
    cat("     TOTAL GROUPS> ", n_groups,  "\n", sep=" ")
    cat("            GROUP> ", vjl_group, "\n", sep=" ")
    cat("             SIZE> ", gp_size,   "\n", sep=" ")
    cat("        V CALL(s)> ", gp_vcall,  "\n", sep=" ")
    cat("        J CALL(s)> ", gp_jcall,  "\n", sep=" ")
    cat("  JUNCTION LENGTH> ", gp_lent,   "\n", sep=" ") 
    cat("            MODEL> ", model,     "\n", sep=" ")
    cat("           METHOD> ", method,    "\n", sep=" ")
    cat("             CDR3> ", cdr3,      "\n", sep=" ")
    cat("         CLONE(s)> ", n_cluster, "\n", sep=" ")
    cat("", "\n", sep=" ")
}
# *****************************************************************************

# *****************************************************************************
logVerbose <- function(out_dir, log_verbose_name,
                       n_groups, vjl_group, model, method, cdr3,
                       gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) {
    cat("     TOTAL GROUPS> ", n_groups,  "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("            GROUP> ", vjl_group, "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("             SIZE> ", gp_size,   "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("        V CALL(s)> ", gp_vcall,  "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("        J CALL(s)> ", gp_jcall,  "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("  JUNCTION LENGTH> ", gp_lent,   "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)   
    cat("            MODEL> ", model,     "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)   
    cat("           METHOD> ", method,    "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("             CDR3> ", cdr3,      "\n", sep=" ", file = file.path(out_dir, log_verbose_name), append=TRUE)
    cat("         CLONE(s)> ", n_cluster, "\n", sep=" ", file=file.path(out_dir, log_verbose_name), append=TRUE)
    cat("", "\n", sep=" ", file=file.path(out_dir, log_verbose_name), append=TRUE)
}
# *****************************************************************************

# *****************************************************************************
# @export
pairwiseMutMatrix <- function(informative_pos, mutMtx, motifMtx) {
    pairwiseMutMatrixRcpp(informative_pos, mutMtx, motifMtx)
}
# *****************************************************************************


#### defineClonesScoper ####

#' Assigning Ig sequences into clonal groups
#'
#' The \code{defineClonesScoper} function provides a computational pipline for assigning Ig 
#' sequences into clonal groups sharing same V gene, J gene, and junction length.
#'
#' @param    db                 data.frame containing sequence data.
#' @param    model              one of the \code{"identical"}, \code{"hierarchical"}, or \code{"spectral"}. 
#'                              See Details for description.
#' @param    method             one of the \code{"nt"}, \code{"aa"}, \code{"single"}, \code{"average"}, 
#'                              \code{"complete"}, \code{"novj"}, or \code{"vj"}. See Details for description.
#' @param    germline_col       character name of the column containing the germline or reference sequence.
#' @param    sequence_col       character name of the column containing input sequences. 
#' @param    junction_col       character name of the column containing junction sequences.
#'                              Also used to determine sequence length for grouping.
#' @param    v_call_col         character name of the column containing the V-segment allele calls.
#' @param    j_call_col         character name of the column containing the J-segment allele calls.
#' @param    clone_col          one of the \code{"CLONE"} or \code{"clone_id"} for the output column name 
#'                              containing the clone ids.
#' @param    targeting_model    \link{TargetingModel} object. Only applicable if \code{model} = \code{"spectral"} 
#'                              and \code{method} = \code{"vj"}. See Details for description.
#' @param    len_limit          \link{IMGT_V} object defining the regions and boundaries of the Ig 
#'                              sequences. If NULL, mutations are counted for entire sequence. Only 
#'                              applicable if \code{model} = \code{"spectral"} and \code{method} = \code{"vj"}.
#' @param    first              specifies how to handle multiple V(D)J assignments for initial grouping. 
#'                              If \code{TRUE} only the first call of the gene assignments is used. 
#'                              If \code{FALSE} the union of ambiguous gene assignments is used to 
#'                              group all sequences with any overlapping gene calls.
#' @param    cdr3               if \code{TRUE} removes 3 nts from both ends of \code{"junction_col"}
#'                              (converts IMGT junction to CDR3 region). if \code{TRUE} remove 
#'                              \code{junction_col}(s) with length less than 7 nts.
#' @param    mod3               if \code{TRUE} removes \code{junction_col}(s) with number of nucleotides not 
#'                              modulus of 3.
#' @param    max_n              The maximum number of N's to permit in the junction sequence before excluding the 
#'                              record from clonal assignment. Note, under model \code{"hierarchical"} and method 
#'                              \code{"single"} non-informative positions can create artifactual links between 
#'                              unrelated sequences. Use with caution. Default is set to be \code{"NULL"} for no action.
#' @param    threshold          the distance threshold for clonal grouping if \code{model} = \code{"hierarchical"}; or 
#'                              the upper-limit cut-off if \code{model} = \code{"spectral"}.
#' @param    base_sim           required similarity cut-off for sequences in equal distances from each other.
#'                              Only applicable if \code{model} = \code{"spectral"}.
#' @param    iter_max	        the maximum number of iterations allowed for kmean clustering step.
#' @param    nstart	            the number of random sets chosen for kmean clustering initialization.
#' @param    nproc              number of cores to distribute the function over.
#' @param    verbose            if \code{TRUE} report a summary of each step cloning process;
#'                              if \code{FALSE} process cloning silently.
#' @param    log_verbose        if \code{TRUE} write verbose logging to a file in \code{out_dir}.
#' @param    out_dir            specify the output directory to save \code{log_verbose}. The input 
#'                              file directory is used if this is not specified.
#' @param    summerize_clones     if \code{TRUE} performs a series of analysis to assess the clonal landscape.
#'                              See Value for description.
#'
#' @return
#' For \code{summerize_clones} = \code{FALSE}, a modified data.frame with clone identifiers in the \code{clone_col} column. 
#' For \code{summerize_clones} = \code{TRUE} returns a list containing:
#' \itemize{
#'      \item   \code{db}:                   modified \code{db} data.frame with clone identifiers in the \code{clone_col} column. 
#'      \item   \code{vjl_group_summ}:       data.frame of clones summary, e.g. size, V-gene, J-gene, junction lentgh,
#'                                           and so on.
#'      \item   \code{inter_intra}:          data.frame containing minimum inter (between) and maximum intra (within) 
#'                                           clonal distances.
#'      \item   \code{eff_threshold}:        effective cut-off separating the inter (between) and intra (within) clonal 
#'                                           distances.
#'      \item   \code{plot_inter_intra}:     ggplot histogram of inter (between) versus intra (within) clonal distances. The 
#'                                           effective threshold is shown with a horizental dashed-line.
#' }
#' If \code{log_verbose} = \code{TRUE}, it will write verbose logging to a file in the current directory or 
#' the specified \code{out_dir}.
#'
#' @details
#' \code{defineClonesScoper} provides a computational platform to explore the B cell clonal 
#' relationships in high-throughput Adaptive Immune Receptor Repertoire sequencing (AIRR-seq) 
#' data sets. Three models are included which perform clustering among sequences of B cell receptors 
#' (BCRs, also referred to as Immunoglobulins, (Igs)) that share the same V gene, J gene and junction length: 
#' \itemize{
#'       \item \code{model} = \code{"identical"}: defines clones among identical junctions. Available \code{method}(s) are:
#'       (1) \code{"nt"} (nucleotide based clustering) and (2) \code{"aa"} (amino acid based clustering).
#'       \item \code{model} = \code{"hierarchical"}: hierarchical clustering-based method for partitioning sequences 
#'       into clones. Availabe agglomeration \code{method}(s) are: (1) \code{"single"}, (2) \code{"average"}, and (3) 
#'       \code{"complete"}. The fixed \code{threshold} (a numeric scalar where the tree should be cut) must be provided.
#'       \item \code{model} = \code{"spectral"}: provides an unsupervised pipline for assigning Ig sequences into clonal 
#'       groups. If \code{method} = \code{"novj"}, clonal relationships are inferred using an adaptive threshold that 
#'       indicates the level of similarity among junction sequences in a local neighborhood. If \code{method} = \code{"vj"}: 
#'       clonal relationships are inferred not only based on the junction region homology, but also takes into account 
#'       the mutation profiles in the V and J segments. \code{germline_col} and \code{sequence_col} must be provided. 
#'       Mutation counts are determined by comparing the input sequences (in the column specified by \code{sequence_col}) 
#'       to the effective germline sequence (calculated from sequences in the column specified by \code{germline_col}). 
#'       Not mandatory, but the influence of SHM hot- and cold-spot biases in the clonal inference process will be noted 
#'       if a SHM targeting model is provided through argument \code{targeting_model} (see \link{createTargetingModel} 
#'       for more technical details).
#' }
#'
#' @examples
#' results <- defineClonesScoper(ExampleDb, 
#'                               model="hierarchical", method="single", 
#'                               threshold=0.15, summerize_clones=TRUE)
#' @export
defineClonesScoper <- function(db,
                               model = c("identical", "hierarchical", "spectral"),
                               method = c("nt", "aa", "single", "average", "complete", "novj", "vj"),
                               germline_col = "GERMLINE_IMGT",
                               sequence_col = "SEQUENCE_IMGT",
                               junction_col = "JUNCTION",
                               v_call_col = "V_CALL",
                               j_call_col = "J_CALL",
                               clone_col = c("clone_id", "CLONE"),
                               targeting_model = NULL,
                               len_limit = NULL,
                               first = FALSE, 
                               cdr3 = FALSE, 
                               mod3 = FALSE,
                               max_n = NULL,
                               threshold = NULL,
                               base_sim = 0.95,
                               iter_max = 1000, 
                               nstart = 1000, 
                               nproc = 1,
                               verbose = FALSE,
                               log_verbose = FALSE,
                               out_dir = ".",
                               summerize_clones = FALSE) {
    
    # Initial checks
    if (!is.data.frame(db)) {
        stop("'db' must be a data frame")
    }
    
    ### get model
    model <- match.arg(model)
    
    ### get method
    method <- match.arg(method)
    
    ### check model andmethod
    if (model == "identical") {
        if (!(method %in% c("nt", "aa"))) {
            stop(paste0("'method' should be one of 'nt' or 'aa' for model '", model, "'.")) 
        }
    } else if (model == "hierarchical") {
        method <- match.arg(method)
        if (!method %in% c("single", "average", "complete")) { 
            stop(paste0("'method' should be one of 'single', 'average', or 'complete' for model '", model, "'.")) 
        }
        if  (is.null(threshold)) {
            stop(paste0("'threshold' should be defined for model '", model, "'.")) 
        }
    } else if (model == "spectral") {
        if (!method %in% c("novj", "vj")) { 
            stop(paste0("'method' should be one of 'novj' or 'vj' for model '", model, "'.")) 
        }
    } else {
        stop(paste0("'model' shoubd be one of 'identical', 'hierarchical', or 'spectral'.")) 
    }
    
    ### get clone column name
    clone_col <- match.arg(clone_col)
    if (!(clone_col %in% c("CLONE", "clone_id"))) {
        stop(paste0("'clone_col' should be one of 'CLONE' or 'clone_id'", ".")) 
    }
    
    ### Check for invalid characters
    valid_chars <- colnames(getDNAMatrix(gap = 0))
    .validateSeq <- function(x) { all(unique(strsplit(x, "")[[1]]) %in% valid_chars) }
    valid_seq <- sapply(db[[junction_col]], .validateSeq)
    not_valid_seq <- which(!valid_seq)
    if (length(not_valid_seq) > 0) {
        stop("invalid sequence characters in the ", junction_col, " column. ",
             length(not_valid_seq)," sequence(s) found.", "\n Valid characters are: '",  valid_chars, "'")
    }
    
    ### check targeting model
    if (!is.null(targeting_model)) { 
        mutabs <- targeting_model@mutability 
    } else {
        mutabs <- NULL
    }
    
    ### check verbose and log
    verbose <- ifelse(verbose, 1, 0)
    log_verbose <- ifelse(log_verbose, 1, 0)
    
    ### temp columns
    temp_cols <- c("VJ_GROUP", "VJL_GROUP", "JUNCTION_L",  "CDR3", "CDR3_L", "CLONE_temp")
    
    ### check for invalid columns
    invalid_cols <- c(clone_col, temp_cols)
    if (any(invalid_cols %in% colnames(db))) {
        stop("Column(s) '", paste(invalid_cols[invalid_cols %in% colnames(db)], collapse = "', '"), "' already exist.",
             "\n Invalid column names are: '", paste(invalid_cols, collapse = "', '"), "'.")
    }
    
    ### Check general required columns
    columns <- c(junction_col, v_call_col, j_call_col) #, fields
    columns <- columns[!is.null(columns)]
    check <- checkColumns(db, columns)
    if (!check) { 
        stop("columns ", paste(columns, collapse = " and "), " are not found.") 
    }
    
    ### Check required columns for method "vj"
    if (model == "spectral" & method == "vj") {
        columns <- c(germline_col, sequence_col) #, fields
        columns <- columns[!is.null(columns)]
        check_ham_mut <- checkColumns(db, columns)
        if (!check_ham_mut) { 
            stop("columns ", paste(columns, collapse = " and "), " are not found.") 
        }
    } 
    
    ### check log verbose
    if (log_verbose) {
        if (is.null(out_dir)) stop("out_dir must be specified.")
        if (!dir.exists(out_dir)) stop("out_dir '", out_dir, "' does not exist.")
        log_verbose_name <- "log_verbose.dat"
        cat(file=file.path(out_dir, log_verbose_name), append=FALSE)
    }
    
    # add junction length column
    db$JUNCTION_L <- stri_length(db[[junction_col]])
    junction_l <- "JUNCTION_L"
    
    ### check for mod3
    # filter mod 3 junction lengths
    if (mod3) {
        n_rmv_mod3 <- sum(db[[junction_l]]%%3 != 0)
        db <- db %>% 
            dplyr::filter(!!rlang::sym(junction_l)%%3 == 0)
    }
    
    ### check for cdr3
    # filter junctions with length > 6
    if (cdr3) {
        n_rmv_cdr3 <- sum(db[[junction_l]] <= 6)
        db <- db %>% 
            dplyr::filter(!!rlang::sym(junction_l) > 6)
        # add cdr3 column
        db$CDR3 <- substr(db[[junction_col]], 4, db[[junction_l]]-3)
        cdr3_col <- "CDR3"
    }
    
    ### check for N's
    # Count the number of 'N's in junction
    if (!is.null(max_n)) {
        n_rmv_N <- sum(str_count(db[[junction_col]], "N") > max_n)
        db <- db %>% 
            dplyr::filter(str_count(!!rlang::sym(junction_col), "N") <= max_n)
    }

    ### Parse V and J columns to get gene
    if (verbose) { cat("ASSIGN VJL GROUPS>", "\n", sep=" ") }
    if (log_verbose) {
        cat("ASSIGN VJL GROUPS>", "\n", sep=" ",
            file = file.path(out_dir, log_verbose_name), append=TRUE) 
        cat("", "\n", sep=" ",
            file = file.path(out_dir, log_verbose_name), append=TRUE) 
    }
    db <- groupGenes(db,
                     v_call = v_call_col,
                     j_call = j_call_col,
                     first = first)
    
    ### groups to use
    groupBy <- c("VJ_GROUP", junction_l)
    
    ### assign group ids to db
    db$VJL_GROUP <- db %>%
        dplyr::group_by(!!!rlang::syms(groupBy)) %>%
        dplyr::group_indices()
    
    ### summary of the groups
    vjl_groups <- db %>% 
        dplyr::group_by(!!rlang::sym("VJL_GROUP")) %>% 
        dplyr::summarise(GROUP_V_CALL = paste(unique(!!rlang::sym(v_call_col)), collapse=","),
                         GROUP_J_CALL = paste(unique(!!rlang::sym(j_call_col)), collapse=","),
                         GROUP_JUNCTION_LENGTH = unique(!!rlang::sym(junction_l)),
                         GROUP_SIZE = n())
    vjl_groups$GROUP_V_CALL <- sapply(1:nrow(vjl_groups), 
                                      function(i){ paste(unique(strsplit(vjl_groups$GROUP_V_CALL[i], split=",")[[1]]), collapse=",") })
    vjl_groups$GROUP_J_CALL <- sapply(1:nrow(vjl_groups), 
                                      function(i){ paste(unique(strsplit(vjl_groups$GROUP_J_CALL[i], split=",")[[1]]), collapse=",") })
    n_groups <- nrow(vjl_groups)
    
    ### Create cluster of nproc size and export namespaces
    if(nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK", outfile = "")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }
    
    ### expoer function to clusters
    if (nproc > 1) { 
        export_functions <- list("passToClustering_lev1", "passToClustering_lev2", "passToClustering_lev3", "passToClustering_lev4",
                                 "findGapSmooth", "infer", "krnlMtxGenerator", "makeAffinity", "laplacianMtx", 
                                 "rangeAtoB", "likelihoods", "pairwiseMutions", "pairwiseMutMatrix",
                                 "printVerbose", "logVerbose")
        parallel::clusterExport(cluster, export_functions, envir=environment())
    }
    
    ### perform clustering for each group
    gp <- NULL
    db_cloned <- foreach(gp = 1:n_groups,
                         .combine = "rbind",
                         .inorder = TRUE,
                         .errorhandling='stop') %dopar% { 
                             # *********************************************************************************
                             # filter each group
                             vjl_group <- vjl_groups$VJL_GROUP[gp]
                             gp_vcall <- vjl_groups$GROUP_V_CALL[gp]
                             gp_jcall <- vjl_groups$GROUP_J_CALL[gp]
                             gp_lent <- vjl_groups$GROUP_JUNCTION_LENGTH[gp]
                             gp_size <- vjl_groups$GROUP_SIZE[gp]
                             db_gp <- db %>%
                                 dplyr::filter(!!rlang::sym("VJL_GROUP") == vjl_group)
                             
                             # pre-check each group
                             passToClust <- TRUE
                             if (method %in% c("novj", "single", "complete", "average")) {
                                 if (length(unique(db_gp[[ifelse(cdr3, cdr3_col, junction_col)]])) == 1) {
                                     idCluster <- rep(1, nrow(db_gp))
                                     n_cluster <- 1 
                                     passToClust <- FALSE
                                 }
                             }
                             if (method == "vj") {
                                 if (length(unique(db_gp[[sequence_col]])) == 1) {
                                     idCluster <- rep(1, nrow(db_gp))
                                     n_cluster <- 1  
                                     passToClust <- FALSE
                                 }
                             }
                             # pass the group for clustering
                             if (passToClust) {
                                 results <- passToClustering_lev1(db_gp,
                                                                  model = model,
                                                                  method = method,
                                                                  germline_col = germline_col,
                                                                  sequence_col = sequence_col,
                                                                  junction_col = junction_col,
                                                                  mutabs = mutabs,
                                                                  len_limit = len_limit,
                                                                  cdr3 = cdr3,
                                                                  cdr3_col = ifelse(cdr3, cdr3_col, NA),
                                                                  v_call_col = v_call_col,
                                                                  j_call_col = j_call_col,
                                                                  threshold = threshold,
                                                                  base_sim = base_sim,
                                                                  iter_max = iter_max, 
                                                                  nstart = nstart)
                                 idCluster <- results$idCluster
                                 n_cluster <- results$n_cluster
                             }
                             
                             if (length(idCluster) == 0 | any(is.na(idCluster))) stop("error in pass")
                             
                             # check verbose
                             if (verbose) { printVerbose(n_groups, vjl_group, model, method, cdr3,
                                                         gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) 
                             }
                             
                             # check log verbose
                             if (log_verbose) { logVerbose(out_dir, log_verbose_name,
                                                           n_groups, vjl_group, model, method, cdr3,
                                                           gp_vcall, gp_jcall, gp_lent, gp_size, n_cluster) }
                             
                             # attache clones
                             db_gp[[clone_col]] <- paste(vjl_group, idCluster, sep="_")   
                             
                             # return result from each proc
                             return(db_gp)
                             # *********************************************************************************
                         }
    
    ### Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }
    
    db_cloned$CLONE_temp <- db_cloned %>%
        dplyr::group_by(!!rlang::sym(clone_col)) %>%
        dplyr::group_indices()
    db_cloned[[clone_col]] <- db_cloned$CLONE_temp
    db_cloned <- db_cloned[order(db_cloned[[clone_col]]), ]
    db_cloned[[clone_col]] <- as.character(db_cloned[[clone_col]])
    
    ### report removed sequences
    if (mod3) {
        if (n_rmv_mod3 > 0) {
            cat("      MOD3 FILTER> ", n_rmv_mod3, "invalid junction length(s) (not mod3) in the", junction_col, "column removed.", "\n", sep=" ")
            if (log_verbose)  { 
                cat("      MOD3 FILTER> ", n_rmv_mod3, "invalid junction length(s) (not mod3) in the", junction_col, "column removed.", "\n", sep=" ",
                    file = file.path(out_dir, log_verbose_name), append=TRUE) 
            }
        }
    }
    if (cdr3) {
        if (n_rmv_cdr3 > 0) {
            cat("      CDR3 FILTER> ", n_rmv_cdr3, "invalid junction length(s) (< 7) in the", junction_col, "column removed.", "\n", sep=" ")
            if (log_verbose)  { 
                cat("      CDR3 FILTER> ", n_rmv_cdr3, "invalid junction length(s) (< 7) in the", junction_col, "column removed.", "\n", sep=" ",
                    file = file.path(out_dir, log_verbose_name), append=TRUE) 
            }
        }
    }
    if (!is.null(max_n)) {
        if (n_rmv_N > 0) {
            cat("     MAX N FILTER> ", n_rmv_N, "invalid junction(s) ( # of N >", max_n, ") in the", junction_col, "column removed.", "\n", sep=" ")
            if (log_verbose)  { 
                cat("     MAX N FILTER> ", n_rmv_N, "invalid junction(s) ( # of N >", max_n, ") in the", junction_col, "column removed.", "\n", sep=" ",
                    file = file.path(out_dir, log_verbose_name), append=TRUE) 
            }
        }
    }
    ### print ending message
    if (verbose) { cat("          CLONING> DONE.", "\n", sep=" ") }
    if (log_verbose)  { 
        cat("          CLONING> DONE.",  "\n", sep=" ", 
            file = file.path(out_dir, log_verbose_name), append=TRUE) 
    }
    
    ### make summary 
    if (summerize_clones) {
        if (verbose) { cat("   SUMMARY CLONES> ...", "\n") }
        if (log_verbose)  { 
            cat("   SUMMARY CLONES> ...",  "\n", sep=" ", 
                file = file.path(out_dir, log_verbose_name), append=TRUE) 
        }
        
        ### vjl group summary
        vjl_groups <- db_cloned %>%
            dplyr::group_by(!!rlang::sym("VJL_GROUP")) %>%
            dplyr::summarise(SIZE = n(),
                             V_CALL = paste(unique(!!rlang::sym(v_call_col)), collapse=","),
                             J_CALL = paste(unique(!!rlang::sym(j_call_col)), collapse=","),
                             JUNCTION_LENGTH = unique(!!rlang::sym(junction_l)),
                             NUMBER_OF_CLONE = length(unique(!!rlang::sym(clone_col))),
                             CLONE_ID = paste(unique(!!rlang::sym(clone_col)), collapse = ","))
        vjl_groups$V_CALL <- sapply(1:nrow(vjl_groups),
                                    function(i){ paste(unique(strsplit(vjl_groups$V_CALL[i], split=",")[[1]]), collapse=",") })
        vjl_groups$J_CALL <- sapply(1:nrow(vjl_groups),
                                    function(i){ paste(unique(strsplit(vjl_groups$J_CALL[i], split=",")[[1]]), collapse=",") })
        
        ### calculate inter and intra distances
        results <- calculateInterVsIntra(db = db_cloned,
                                         clone_col = clone_col,
                                         vjl_groups = vjl_groups,
                                         junction_col = junction_col,
                                         cdr3 = cdr3,
                                         cdr3_col = ifelse(cdr3, cdr3_col, NA),
                                         nproc = nproc,
                                         verbose = verbose)
        
        ### plot inter and intra distances
        p_results <- plotInterVsIntra(data = results$inter_intra)
        
        if (verbose)  { cat("   SUMMARY CLONES> DONE.", "\n") }
        if (log_verbose)  { 
            cat("   SUMMARY CLONES> DONE.",  "\n", sep=" ", 
                file = file.path(out_dir, log_verbose_name), append=TRUE) 
        }
    }
    
    ### remove extra columns
    db_cloned <- db_cloned[, !(names(db_cloned) %in% temp_cols)]
    
    # return results
    if (summerize_clones) {
        return_list <- list("db" = db_cloned,
                            "vjl_group_summ" = vjl_groups,
                            "inter_intra" = results$inter_intra,
                            "eff_threshold" = p_results$eff_threshold,
                            "plot_inter_intra" = p_results$plot_inter_intra)    
        return(return_list)
    } else {
        return(db_cloned)
    }
}
# *****************************************************************************

# *****************************************************************************
passToClustering_lev1 <- function (db_gp, 
                                   model = c("identical", "hierarchical", "spectral"),
                                   method = c("nt", "aa", "single", "average", "complete", "novj", "vj"),
                                   germline_col = "GERMLINE_IMGT",
                                   sequence_col = "SEQUENCE_IMGT",
                                   junction_col = "JUNCTION",
                                   mutabs = NULL,
                                   len_limit = NULL,
                                   cdr3 = FALSE,
                                   cdr3_col = NA,
                                   v_call_col = "V_CALL",
                                   j_call_col = "J_CALL",
                                   threshold = NULL,
                                   base_sim = 0.95,
                                   iter_max = 1000, 
                                   nstart = 1000) {
    ### get model
    model <- match.arg(model)
    
    ### get method
    method <- match.arg(method)
    
    ### number of sequences
    n <- nrow(db_gp)
    
    ### begin clustering
    if (model == "identical") {
        #*************************
        seq_col <- ifelse(cdr3, cdr3_col, junction_col)
        if (method == "aa") {
            db_gp[[seq_col]] <- translateDNA(db_gp[[seq_col]])
        }
        idCluster <- db_gp %>% 
            dplyr::group_by(!!rlang::sym(seq_col)) %>% 
            group_indices()
        n_cluster <- length(unique(idCluster))
        eigen_vals <- rep(0, n)
        #*************************
    } else if (model == "hierarchical") {
        #*************************
        # get sequences
        seqs <- db_gp[[ifelse(cdr3, cdr3_col, junction_col)]]
        # find unique seqs
        df <- as.data.table(seqs)[, list(list(.I)), by = seqs]
        n_unq <- nrow(df)
        ind_unq <- df$V1
        seqs_unq <- df$seqs
        # calculate distance matrix
        dist_mtx <- pairwiseDist(seq = seqs_unq, 
                                 dist_mat = getDNAMatrix(gap = 0))
        # calculate normalization factor
        junc_length <- unique(stri_length(seqs_unq))
        # perform hierarchical clustering
        hc <- hclust(as.dist(dist_mtx/junc_length), method = method)
        # cut the tree
        idCluster_unq <- cutree(hc, h = threshold)
        # back to reality
        idCluster <- rep(NA, n)
        for (i in 1:n_unq) {
            idCluster[ind_unq[[i]]] <- idCluster_unq[i]
        }
        n_cluster <- length(unique(idCluster))
        eigen_vals <- rep(0, n)
        #*************************
    } else if (model == "spectral") {
        #*************************
        if (method == "vj") {
            # get required info based on the method
            germs <- db_gp[[germline_col]]
            seqs <- db_gp[[sequence_col]]
            juncs <- db_gp[[ifelse(cdr3, cdr3_col, junction_col)]]
            junc_length <- unique(stri_length(juncs))
            # find unique seqs
            df <- as.data.table(seqs)[, list(list(.I)), by = seqs]
            n_unq <- nrow(df)
            ind_unq <- df$V1
            seqs_unq <- df$seqs
            # find corresponding unique germs and junctions
            germs_unq <- sapply(1:n_unq, function(x){ unique(germs[ind_unq[[x]]]) })
            juncs_unq <- sapply(1:n_unq, function(x){ unique(juncs[ind_unq[[x]]]) })
            # calculate unique junctions distance matrix
            dist_mtx <- pairwiseDist(seq = juncs_unq, 
                                     dist_mat = getDNAMatrix(gap = 0))
            # count mutations from unique sequence imgt
            results <- pairwiseMutions(germ_imgt = germs_unq, 
                                       seq_imgt = seqs_unq,
                                       junc_length = junc_length, 
                                       len_limit = len_limit,
                                       cdr3 = cdr3,
                                       mutabs = mutabs)
            tot_mtx <- results$pairWiseTotalMut
            sh_mtx <- results$pairWiseSharedMut
            mutab_mtx <- results$pairWiseMutability
            # calculate likelihhod matrix
            lkl_mtx <- likelihoods(tot_mtx = tot_mtx, 
                                   sh_mtx = sh_mtx, 
                                   mutab_mtx = mutab_mtx)
            # calculate weighted matrix  
            disim_mtx <- dist_mtx * (1.0 - lkl_mtx)
            # check if new method made any changes, otherwise go back to method "novj"
            if (all(disim_mtx == dist_mtx)) {
                # get required info based on the method
                seqs <- db_gp[[ifelse(cdr3, cdr3_col, junction_col)]]
                junc_length <- unique(stri_length(seqs))
                # find unique seqs
                df <- as.data.table(seqs)[, list(list(.I)), by = seqs]
                n_unq <- nrow(df)
                ind_unq <- df$V1
                seqs_unq <- df$seqs
                # calculate unique seuences distance matrix
                disim_mtx <- pairwiseDist(seq = seqs_unq, 
                                          dist_mat = getDNAMatrix(gap = 0))
            } 
        } else if (method == "novj") {
            # get required info based on the method
            seqs <- db_gp[[ifelse(cdr3, cdr3_col, junction_col)]]
            junc_length <- unique(stri_length(seqs))
            # find unique seqs
            df <- as.data.table(seqs)[, list(list(.I)), by = seqs]
            n_unq <- nrow(df)
            ind_unq <- df$V1
            seqs_unq <- df$seqs
            # calculate unique seuences distance matrix
            disim_mtx <- pairwiseDist(seq = seqs_unq, 
                                      dist_mat = getDNAMatrix(gap = 0))
        }
        ### pass to clustering pipeline
        result <- passToClustering_lev2(mtx = disim_mtx, 
                                        junc_length = junc_length,
                                        threshold = threshold, 
                                        base_sim = base_sim, 
                                        iter_max = iter_max, 
                                        nstart = nstart)
        idCluster_unq <- result$idCluster
        eigen_vals <- result$eigen_vals
        ### back to reality
        idCluster <- rep(NA, n)
        for (i in 1:n_unq) {
            idCluster[ind_unq[[i]]] <- idCluster_unq[i]
        }
        n_cluster <- length(unique(idCluster))
        #*************************
    }
    
    ### retrun results
    return_list <- list("model" = model,
                        "method" = method,
                        "n_cluster" = n_cluster, 
                        "idCluster" = idCluster, 
                        "eigen_vals" = eigen_vals)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
# check special case if all distances (off diagonal elements) are the same. 
# (this also includes cases with only two sequences)
passToClustering_lev2 <- function(mtx, 
                                  junc_length = NULL, 
                                  threshold = NULL, 
                                  base_sim = 0.95, 
                                  iter_max = 1000, 
                                  nstart = 1000) {
    ### constants
    n <- nrow(mtx)
    bs <- (1 - base_sim)*junc_length
    off_diags_nuq <- unique(mtx[row(mtx) != col(mtx)])
    ### check special cases
    if (n == 1) {
        idCluster <- 1    # all in same clone
        eigen_vals <- 0
    } else if (length(off_diags_nuq) == 1) {   # seqs have equal-distances from each other
        if (off_diags_nuq > bs) { 
            idCluster <- 1:n         # all singletons
            eigen_vals <- rep(0, n)
        } else {
            idCluster <- rep(1, n)    # all in same clone
            eigen_vals <- rep(0, n)
        }
    } else if (max(mtx) <= bs) {
        idCluster <- rep(1, n)    # all in same clone
        eigen_vals <- rep(0, n)
    } else {
        results <- passToClustering_lev3(mtx = mtx, 
                                         junc_length = junc_length,
                                         threshold = threshold, 
                                         iter_max = iter_max, 
                                         nstart = nstart)
        idCluster <- results$idCluster
        eigen_vals <- results$eigen_vals
    }
    ### return list
    return_list <- list("idCluster" = idCluster, 
                        "eigen_vals" = eigen_vals)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
passToClustering_lev3 <- function(mtx, 
                                  junc_length = NULL, 
                                  threshold = NULL, 
                                  iter_max = 1000, 
                                  nstart = 1000){
    n <- nrow(mtx)
    ### set seed for reproducibility
    set.seed(12345)
    ### check the minimum number of data points requirement
    if (n < 3) stop("SCOPer needs at least 3 unique data points")
    ### calculate the krnl matrix
    krnl_mtx <- krnlMtxGenerator(mtx = mtx)
    ### calculate the affinity matrix
    if (!is.null(threshold)) {
        aff_mtx <- makeAffinity(mtx_o = mtx, 
                                mtx_k = round(krnl_mtx, 3),
                                thd = threshold*junc_length) 
    } else {
        aff_mtx <- round(krnl_mtx, 3)
    }
    ### affinity matrix is diagonal. Each sequence belongs to a singlton clone.
    if (all(aff_mtx[!diag(nrow(aff_mtx))] == 0)) { 
        return_list <- list("idCluster" = c(1:n),
                            "eigen_vals" = rep(0,n))
        return(return_list)
    }
    ### calculate the laplacian matrix
    L <- laplacianMtx(entry = aff_mtx)
    ### scale
    # Each column forced to have zero-mean and unit-variance
    # Na is produced if L has a zero rows/columns, meaning aff_mtx has had unitary rows/columns.
    scaleL <- scale(L, center = TRUE, scale = TRUE)
    ### correlation
    # A value of near or equal to 0 implies little or no linear relationship.
    # In contrast, the closer comes to 1 or -1, the stronger the linear relationship.
    # Using "pairwise.complete.obs" cor produce Na for that particular observations (all elements of both rows and clolumns)
    # which(is.na(corL), arr.ind = T)
    corL <- cor(scaleL, method = "pearson", use = "pairwise.complete.obs")        
    corL[is.na(corL)] <- 0  
    ### correlations less than or equal 0.05% are ignored
    corL <- round(corL, 3)
    ### eigen-decomposition
    eigens <- eigen(corL, symmetric = TRUE)
    eigen_vals <- rev(eigens$values)
    eigen_vals <- rangeAtoB(eigen_vals, 0, 1) 
    if (all(eigen_vals == 0)) {  # Each sequence belongs to a singlton clone.
        return_list <- list("idCluster" = c(1:n),
                            "eigen_vals" = rep(0,n))
        return(return_list)
    }
    ### k upper bound
    eigenValDens <- density(eigen_vals)
    k_up <- sum(eigen_vals < eigenValDens$x[which.max(eigenValDens$y)]) + 1
    if (is.na(k_up) | k_up %in% c(1,2)) k_up <- n
    # determine the number of clusters using log distance
    logEigenValsDiff <- sapply(2:k_up, function(x){ log10(eigen_vals[x]) - log10(eigen_vals[x-1]) })
    nas <- sum(is.nan(logEigenValsDiff) | is.na(logEigenValsDiff) | is.infinite(logEigenValsDiff))
    if (nas == length(logEigenValsDiff)) {
        k <- nas
    } else {
        k <- nas + ifelse(nas > 0, which.max(logEigenValsDiff[-(1:nas)]), which.max(logEigenValsDiff))
    }
    if (is.na(k)) stop("level3:: problem at number of clusters determination")
    ### pick k smalest eigenvectors
    # The vectors are col-normalized to unit length
    eigenVecs <- eigens$vectors
    eigenVecs <- eigenVecs[, (n-k+1):(n)]
    ### kmeans clustering
    idCluster <- kmeans(x = round(eigenVecs, 6), 
                        centers = k, 
                        iter.max = iter_max, 
                        nstart = nstart)$cluster
    ### check if idclusters and affinity matrix agrees
    idCluster <- passToClustering_lev4(aff_mtx = aff_mtx, 
                                       idCluster = idCluster)
    ### return results
    return_list <- list("idCluster" = idCluster,
                        "eigen_vals" = eigen_vals)
    return(return_list)
}
# *****************************************************************************

# *****************************************************************************
# check affinity matrix and clusters id 
# aff_mtx_sub[which(aff_id_sub %in% gr_ls[[4]]),]
passToClustering_lev4 <- function(aff_mtx, idCluster) {
    n <- nrow(aff_mtx)
    new_idCluster <- rep(NA, n)
    new_k <- 0
    id_unq <- unique(idCluster)
    k <- length(unique(idCluster))
    for (z in 1:k) {
        aff_id_sub <- which(idCluster == id_unq[z])
        if (length(aff_id_sub) == 1) {
            new_k <- new_k + 1
            new_idCluster[aff_id_sub] <- new_k    
        } else {
            aff_mtx_sub <- aff_mtx[aff_id_sub, aff_id_sub]
            n <- nrow(aff_mtx_sub)
            rows_ls <- rep(list(NULL), n)
            for (i in 1:n) {
                rows_ls[[i]] <- aff_id_sub[aff_mtx_sub[, i] != 0]
            }
            gr_ls <- rep(list(NA), n)
            for (y in 1:n) {
                ids <- rows_ls[[y]]
                l <- length(gr_ls[!is.na(gr_ls)])
                for (x in 1:l) {
                    if (1 > l) break
                    if (length(intersect(ids, gr_ls[[x]])) > 0) {
                        ids <- union(ids, gr_ls[[x]])
                        gr_ls[[x]] <- NA
                    }
                }
                gr_ls[[y]] <- ids
                gr_ls <- c(gr_ls[!is.na(gr_ls)], gr_ls[is.na(gr_ls)])
            }
            gr_ls <- gr_ls[!is.na(gr_ls)]
            for (i in 1:length(gr_ls)) {
                new_k <- new_k + 1
                new_idCluster[gr_ls[[i]]] <- new_k    
            }
        }
    }
    return(new_idCluster)
}
# TO CHECK USE aff_mtx_sub[, which(aff_id_sub == ???)]
