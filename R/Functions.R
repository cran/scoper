#### Classes ####

#' Output of the analyzeClones function
#'
#' \code{ClonalAnalysis} contains output from the \link{analyzeClones} function.
#' It includes infromation to interpret clonal assignment performance.
#'
#' @slot   threshold              cut-off separating the inter (within) and intra (between)
#'                                clonal distances.
#' @slot   inter_intra           data.frame containing all inter and intra clonal distances.
#' @slot   plot_inter_intra       density plot of inter versus intra clonal distances. The threshold is
#'                                shown with a horizental dashed-line.
#' @slot   neighborhoods          a numeric vector containing scale parameters used in spectral
#'                                clustering process.
#' @slot   plot_neighborhoods      histogram of neighborhoods. The threshold is shown with a vertical
#'                                dashed-line.
#'
#' @seealso      \link{analyzeClones}
#'
#' @name         ClonalAnalysis-class
#' @rdname       ClonalAnalysis-class
#' @aliases      ClonalAnalysis
#' @exportClass  ClonalAnalysis
setClass("ClonalAnalysis",
         slots = c(threshold="numeric",
                   inter_intra="list",
                   plot_inter_intra="list",
                   neighborhoods="numeric",
                   plot_neighborhoods="list"))

####  sub-functions ####

# filter function
filterFunction <- function(groupBy, filterBy) {
    paste(paste(groupBy, paste("'",filterBy,"'", sep=""), sep="=="), collapse=" & ")
}

# epsilon calculator by "null"
epsilonNullD <- function(d, vec) {
    return(sd(vec))
}

# epsilon calculator by "assign"
assign <- function(d, vec) {
    vec <- vec[vec<=d]
    std <- ifelse (length(vec) == 1, 0, sd(vec))
    return(std)
}

# epsilon calculator by "infer"
infer <- function(d, vec) {
    vec <- sort(vec)
    diffVec <- diff(vec)
    if (length(unique(diffVec[diffVec > 0])) == 1) {
        d <- ceiling(mean(vec))
    } else {
        id <- which.max(diffVec)
        d <- vec[id]
    }
    vec <- vec[vec<=d]
    std <- ifelse (length(vec) == 1, 0, sd(vec))
    return(std)
}

# kernel matrix calculator
krnlMtxGenerator <- function(mtx, neighborhood=c("assign", "infer"), d=NULL) {
    # Radial basis function kernel: In the Gaussian Kernel if two points are
    # close then K_ij≈1 and when two points are far apart then Kij≈0
    fun <- match.arg(neighborhood)
    n <- dim(mtx)[1]
    # calculate epsilons
    epsilon <- rep(0, length=n)
    for (i in 1:n) {
        epsilon[i] <- apply(as.matrix(mtx[i,]), 2 , FUN=fun, d=d)
        #as.matrix produces a single-column matrix, therefore we use 2 in apply.
    }
    # calculate kernel matrix
    krnl_mtx <- matrix(data=1, nrow=n, ncol=n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            krnl_mtx[i,j] <- exp(-mtx[i,j]^2/(2*epsilon[i]*epsilon[j]))
            krnl_mtx[j,i] <- krnl_mtx[i,j]
        }
    }
    # if mtx[i,j] and epsilon == 0
    krnl_mtx[is.nan(krnl_mtx)] <- 1
    # Disconnect those edges with distance larger than d
    if (neighborhood == "assign") krnl_mtx[mtx > d] <- 0
    return(krnl_mtx)
}

 # affinity matrix calculator
makeAffinity <- function(mtx, n.neighboors=2) {
    n <- nrow(mtx)
    if (n.neighboors >= n) {  # fully connected
        aff_mtx <- mtx
    } else {
        aff_mtx <- matrix(data=0, ncol=n, nrow=n)
        for(i in 1:n) {
            # for each line only connect to those points with larger similarity
            best.similarities <- sort(mtx[i,], decreasing=TRUE)[1:n.neighboors]
            for (s in best.similarities) {
                j <- which(mtx[i,] == s)
                aff_mtx[i,j] <- mtx[i,j]
                aff_mtx[j,i] <- mtx[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
            }
        }
    }
    return(aff_mtx)
}

# laplacian matrix calculator
laplacian_mtx <- function(entry) {
    # Calculate unnormalised Laplacian matrix and its eigenfunctions
    D <- diag(apply(entry, 1, sum))
    L <- D - entry
    return(L)
}

# order of magnitude calculator
log10_ceiling <- function(x) {
    ceiling(log10(x))-1
}

# lowest integer
floor_dec <- function(x, level=1) {
    round(x - 5*10^(-level-1), level)
}

# normalizing function
normalize <- function(x) {
    return ((x-min(x)) / (max(x)-min(x)))
}

# Intersection Function
intersectPoint <- function(t,
                           func1.0 = 0, func1.1 = 0, func1.2 = 0,
                           func2.0 = 0, func2.1 = 0, func2.2 = 0) {
    fit1 <- func1.0*dnorm(t, mean = func1.1, sd = func1.2)
    fit2 <- func2.0*dnorm(t, mean = func2.1, sd = func2.2)
    return(fit1 - fit2)
}

# spectral clustering algorithm
spectralClustering <- function(entrySeq,
                               id,
                               similarity = NULL,
                               neighborhood = c("assign", "infer"),
                               iter_max = 25,
                               nstart = 25) {
    # constants
    nSeq <- length(entrySeq)
    uniqueSeq <- unique(entrySeq)
    nUniqueSeq <- length(uniqueSeq)
    l <- unique(nchar(uniqueSeq))
    uniqueId <-  unique(id)
    # scale parameteres
    neighborhood <- match.arg(neighborhood)
    if (neighborhood == "assign") {
        if (is.null(similarity)) stop("similarity needs to be assigned.")
        d <- round((1-similarity)*l)
    } else if (neighborhood == "infer") {
        d <- NULL
    }
    # calculate distance matrix
    dist_mtx <- pairwiseDist(seq=uniqueSeq, dist_mat=getDNAMatrix(gap=0))
    # calculate kernel matrix
    krnl_mtx <- krnlMtxGenerator(mtx=dist_mtx, neighborhood=neighborhood, d=d)
    # calculate affinity matrix. n.neighboors could be ceiling(sqrt(nUniqueSeq))
    aff_mtx <- makeAffinity(mtx=round(krnl_mtx, 2), n.neighboors=nUniqueSeq)
    # clustering
    unique_idCluster <- kMeanClustering(entry=aff_mtx, iter_max=iter_max, nstart=nstart)
    # back to reality
    idCluster <- rep(NA, nSeq)
    for (i in 1:nUniqueSeq) {
        j <- which(id == uniqueId[i])
        idCluster[j] <- unique_idCluster[i]
    }
    # return results
    return(idCluster)
}

# kmean clustering algorithm
kMeanClustering <- function (entry,
                             iter_max = 25,
                             nstart = 25) {
    n <- nrow(entry)
    if (all(entry[!diag(nrow(entry))] == 0)) {
        # affinity matrix is diagonal. Each sequence belongs to a singlton clone.
        return(c(1:n))
    } else {
        # every real symmetric matrix is Hermitian, and therefore all its eigenvalues are real.
        L <- laplacian_mtx(entry)
        evL <- eigen(L, symmetric=TRUE)
        eigenVecs <- round(evL$vectors, 5)
        eigenVals <- floor_dec(evL$values, 3)
        egvRev <- rev(eigenVals)
        nonZero.ID <- which(egvRev > 0)
        k <- min(nonZero.ID)-1
        # k <- which.max(diff(egvRev))   # slope: forward differencing
        kentry <- eigenVecs[ ,(n-k+1):n]
        clust <- kmeans(kentry, centers=k, iter.max=iter_max, nstart=nstart)
        return(clust$cluster)
    }
}

# inter-clone-distance vs intra-clone-distance
calculateInterVsIntra <- function(db,
                                  junction = "JUNCTION",
                                  v_call = "V_CALL",
                                  j_call = "J_CALL",
                                  clone = "CLONE",
                                  first = TRUE,
                                  cdr3 = FALSE,
                                  nproc = 1,
                                  progress = FALSE) {

    # make temproary clone column
    db$CLONE_temp <- as.character(db[[clone]])

    # find unique groups of sequences with same vgene, jgene, and junction length
    uniqueGroups <- db %>%
        dplyr::group_by_(.dots = c("VJ_GROUP", "L")) %>%
        # dplyr::summarise(CLONE_temp=paste(unique(CLONE_temp), collapse=",")) %>%
        dplyr::summarise_(CLONE_temp=interp(~paste(unique(x), collapse=","), x=as.name("CLONE_temp")))

    n_groups <- nrow(uniqueGroups)

    # Create cluster of nproc size and export namespaces
    if(nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }

    # check the progressbar
    if (progress) {
        cat("INTER AND INTRA DISTANCES ANALYSIS> ", "\n")
        pb <- progressBar(n_groups)
    }

    k <- NULL
    # open dataframes
    vec_ff <- foreach(k=1:n_groups,
                      .combine="c",
                      #.packages = c("dplyr"),
                      .errorhandling='stop') %dopar% {
                          clones <- strsplit(uniqueGroups$CLONE_temp[k], split=",")[[1]]
                          l <- uniqueGroups$L[k]
                          n_clones <- length(clones)
                          seqs <- db$JUNC_temp[db$CLONE_temp %in% clones]
                          names(seqs) <- db$CLONE_temp[db$CLONE_temp %in% clones]
                          # cat("group=", k, "/", n_groups, ", number of clones=", n_clones, ", number of sequences=", length(seqs), "\n", sep="")
                          seqs_db <- data_frame(value = seqs, name = names(seqs)) %>%
                              dplyr::group_by_(.dots=c("name", "value")) %>% # alternatively: group_by(name) if name value pair is always unique
                              slice(1) %>%
                              ungroup()
                          seqs <- seqs_db$value
                          names(seqs) <- seqs_db$name
                          dist_mtx <- pairwiseDist(seqs, dist_mat=getDNAMatrix(gap=0))
                          # prealoocate a dataframe
                          nrow_f <- n_clones + n_clones*(n_clones-1)/2
                          vec_f <- rep(NA, nrow_f)
                          # minimum and maximum distance in each clone
                          n <- 0
                          if (n_clones == 1) {
                              if (!all(dist_mtx == 0)) {
                                  n <- n+1
                                  vec_f[n] <- max(dist_mtx)/l
                                  names(vec_f)[n] <- paste(clones[1], "NA", "inter", sep="_")
                              }
                          } else {
                              for (i in 1:(n_clones-1)) {
                                  xx <- dist_mtx[rownames(dist_mtx) == clones[i], colnames(dist_mtx) == clones[i]]
                                  if (!all(xx == 0)) {
                                      n <- n+1
                                      vec_f[n] <- max(xx)/l
                                      names(vec_f)[n] <- paste(clones[i], "NA", "inter", sep="_")
                                  }
                                  for (j in (i+1):n_clones) {
                                      xy <- dist_mtx[rownames(dist_mtx) == clones[i], colnames(dist_mtx) == clones[j]]
                                      xy <- xy[xy>0]
                                      if (length(xy) != 0) {
                                          n <- n+1
                                          vec_f[n] <- min(xy)/l
                                          names(vec_f)[n] <- paste(clones[i], clones[j], "intra", sep="_")
                                      }
                                  }
                              }
                              yy <- dist_mtx[rownames(dist_mtx) == clones[j], colnames(dist_mtx) == clones[j]]
                              if (!all(yy == 0)) {
                                  n <- n+1
                                  vec_f[n] <- max(yy)/l
                                  names(vec_f)[n] <- paste(clones[j], "NA", "inter", sep="_")
                              }
                          }
                          # Update progress
                          if (progress) { pb$tick() }
                          # remove all na's
                          vec_f <- vec_f[!is.na(vec_f)]
                          # return result from each proc
                          if (length(vec_f) == 0) return(NULL)
                          return(vec_f)
                      }
    # Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }

    # convert to a data.frame
    db_dff <- data.frame(keyName=names(vec_ff), VALUE=vec_ff, row.names=NULL)
    db_dff$LABEL <- "inter"
    db_dff$LABEL[grepl("intra", db_dff$keyName)] <- "intra"
    clones_xy <- data.frame(matrix(unlist(stri_split_fixed(db_dff$keyName, "_", n=3)),
                                   nrow=nrow(db_dff),
                                   byrow=T),
                            stringsAsFactors=FALSE)
    db_dff <- cbind(clones_xy, db_dff)
    db_dff$keyName <- NULL
    colnames(db_dff)[colnames(db_dff) == "X1"] <- "CLONE_X"
    colnames(db_dff)[colnames(db_dff) == "X2"] <- "CLONE_Y"
    colnames(db_dff)[colnames(db_dff) == "X3"] <- "empty"
    db_dff$empty <- NULL
    db_dff$CLONE_Y[grepl("inter", db_dff$CLONE_Y)] <- NA
    interVsIntra <- list()
    interVsIntra[[length(interVsIntra)+1]] <- db_dff

    # find threshold
    db.summ <- db_dff %>%
        dplyr::group_by_(.dots = "LABEL") %>%
        dplyr::summarise_(MEAN=interp(~mean(x, na.rm = TRUE), x=as.name("VALUE")),
                          SD=interp(~sd(x, na.rm = TRUE), x=as.name("VALUE")))

    func1.1 <- db.summ$MEAN[db.summ$LABEL == "inter"]
    func1.2 <- db.summ$SD[db.summ$LABEL == "inter"]
    func2.1 <- db.summ$MEAN[db.summ$LABEL == "intra"]
    func2.2 <- db.summ$SD[db.summ$LABEL == "intra"]
    minInt <- 0
    maxInt <- 1
    intxn <- uniroot(intersectPoint, interval = c(minInt, maxInt), tol=1e-8, extendInt="yes",
                     func1.0=1, func1.1=func1.1, func1.2=func1.2,
                     func2.0=1, func2.1=func2.1, func2.2=func2.2)
    threshold <- round(intxn$root, 2)

    ClonalAnalysis <- new("ClonalAnalysis",
                                inter_intra=interVsIntra,
                                threshold=threshold)
    return(ClonalAnalysis)
}

# plot inter-clone-distance vs intra-clone-distance
plotInterVsIntra <- function(db, threshold) {

    . <- NULL

    # Check for valid columns
    columns <- c("VALUE", "LABEL")
    columns <- columns[!is.null(columns)]
    check <- checkColumns(db, columns)
    if (check != TRUE) { stop(check) }
    db <- select(db, c("VALUE", "LABEL"))

    pdat <- db %>%
        group_by_(.dots = "LABEL") %>%
        do(data.frame(loc = density(.$VALUE, na.rm = TRUE)$x,
                      dens = density(.$VALUE, na.rm = TRUE)$y
        )
        )

    pdat <- pdat %>%
        dplyr::group_by_(.dots = "LABEL") %>%
        dplyr::mutate_(dens_norm=interp(~normalize(x), x=as.name("dens")))

    # Flip and offset densities for the groups
    pdat$dens_norm <- ifelse(pdat$LABEL == 'inter', pdat$dens_norm*(-1), pdat$dens_norm)
    # fill color
    fill_manual <- c("inter"="grey30",
                     "intra"="grey60")
    # labels
    # mg <- "Mean \u00b1 SD"
    # Encoding(mg) <- "UTF-8"
    # pm <-"\u00b1"
    # Encoding(pm)<-"UTF-8"
    # mg1_txt <- data.frame(X=-0.9, Y=threshold-0.025, LAB=paste(mg))
    # db1_txt <- data.frame(X=-0.9, Y=threshold-0.05, LAB=paste(meanInter, pm, sdInter))
    # mg2_txt <- data.frame(X=0.9, Y=threshold+0.05, LAB=paste(mg))
    # db2_txt <- data.frame(X=0.9,  Y=threshold+0.025, LAB=paste(meanIntra, pm, sdIntra))

    hline_size <- 1.0
    # txt_size <- 5.0
    # plot
    p <- ggplot(pdat) +
        baseTheme() +
        theme(axis.title.x = element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y=element_text(size=12),
              axis.text.y=element_text(size=11)) +
        ggtitle(paste("Effective threshold=", threshold)) +
        ylab("Normalized hamming distance") +
        scale_fill_manual(name="",
                          values=fill_manual,
                          labels = c("inter"="maximum-distance\nwithin clones  ",
                                     "intra"="minimum-distance\nbetween clones  ")) +
        geom_polygon(aes_string(x="dens_norm", y="loc", fill="LABEL")) +
        geom_hline(yintercept=threshold, color="grey30", linetype=2, size=hline_size) #+
    # geom_text(data=mg1_txt, aes(x=X, y=Y, label=LAB), size=txt_size, inherit.aes = FALSE) +
    # geom_text(data=db1_txt, aes(x=X, y=Y, label=LAB), size=txt_size, inherit.aes = FALSE) +
    # geom_text(data=mg2_txt, aes(x=X, y=Y, label=LAB), size=txt_size, inherit.aes = FALSE) +
    # geom_text(data=db2_txt, aes(x=X, y=Y, label=LAB), size=txt_size, inherit.aes = FALSE)
    return(p)
}

# neighborhoods
calculateNeighborhoods <- function(db,
                                   junction = "JUNCTION",
                                   v_call = "V_CALL",
                                   j_call = "J_CALL",
                                   similarity = NULL,
                                   neighborhood = c("assign", "infer"),
                                   first = TRUE,
                                   cdr3 = FALSE,
                                   progress = FALSE) {
    # groups to use
    groupBy <- c("VJ_GROUP", "L")
    # find unique groups of sequences with same vgene, jgene, and junction length
    uniqueGroups <- data.frame(unique(db[, groupBy]))
    colnames(uniqueGroups) <- groupBy
    rownames(uniqueGroups) <- NULL
    n_groups <- nrow(uniqueGroups)

    # scale parameteres
    neighborhood <- match.arg(neighborhood)
    if (neighborhood == "assign") {
        if (is.null(similarity)) stop("similarity needs to be assigned.")
        d <- round((1-similarity)*l)
    } else if (neighborhood == "infer") {
        d <- NULL
    }

    # check the progressbar
    if (progress) {
        cat("\n")
        cat("NEIGHBORHOOD ANALYSIS> ", "\n")
        pb <- progressBar(n_groups)
    }

    # sigma anlysis
    eps <- rep(NA, nrow(db))
    idx <- 1
    for (i in 1:n_groups) {
        # print(paste(n_groups, i, sep=" "))
        db_group <- filter_(db, ~VJ_GROUP == uniqueGroups$VJ_GROUP[i], ~L == uniqueGroups$L[i])

        if (nrow(db_group) == 1) {
            # Update progress
            if (progress) { pb$tick() }
            next
        }
        l <- unique(db_group$L)
        entrySeq <- db_group$JUNC_temp
        mtx <- pairwiseDist(seq=entrySeq, dist_mat=getDNAMatrix(gap=0))
        n <- dim(mtx)[1]
        epsilon <- rep(0, length=n)
        for (j in 1:n) {
            epsilon[j] <- apply(as.matrix(mtx[j,]), 2 , FUN=neighborhood, d=d)
        }
        eps[idx:(idx+n-1)] <- epsilon/l
        idx <- idx+n
        # Update progress
        if (progress) { pb$tick() }
    }
    # remove all na's
    eps <- eps[!is.na(eps)]
    return(eps)
}

# plot neighborhoods
plotNeighborhoods <- function(sigmas, threshold = NULL) {
    sigmas_df <- data.frame(sigmas=sigmas)
    vline_size <- 1.0
    binwidth <- 0.01
    center <- 0.005
    p <- ggplot(sigmas_df) +
        baseTheme() +
        theme(axis.text=element_text(size=11),
              axis.title=element_text(size=12)) +
        ggtitle(paste("Effective threshold=", threshold)) +
        xlab("Normalized hamming distance") +
        ylab("Density") +
        geom_histogram(aes_string(x="sigmas", y="..density.."),
                       binwidth=binwidth, center=center, alpha=1.0, fill="black", color="white", size=0.25)
    if (!is.null(threshold)) p <- p + geom_vline(xintercept=threshold, linetype=2, size=vline_size)
    return(p)
}


#### defineClonesScoper ####
#' Assigning Ig sequences into clonal groups
#'
#' The \code{defineClonesScoper} function provides an unsupervised pipline for assigning Ig sequences into
#' clonal groups sharing same V gene, J gene, and junction length.
#'
#' @param    db              data.frame with Change-O style columns containing sequence data.
#' @param    junction        name of the column containing nucleotide sequences to compare.
#'                           Also used to determine sequence length for grouping.
#' @param    v_call          name of the column containing the V-segment allele calls.
#' @param    j_call          name of the column containing the J-segment allele calls.
#' @param    first           if \code{TRUE} only the first call of the gene assignments
#'                           is used. if \code{FALSE} the union of ambiguous gene
#'                           assignments is used to group all sequences with any
#'                           overlapping gene calls.
#' @param    cdr3            if \code{TRUE} remove 3 nts from both ends of \code{junction}
#'                           (converts IMGT junction to CDR3 region). if \code{TRUE} remove \code{junction}(s)
#'                           with length less than 7 nts.
#' @param    mod3            if \code{TRUE} remove \code{junction}(s) with number of nucleotides not modulus of 3.
#' @param    iter_max	     the maximum number of iterations allowed for kmean clustering step.
#' @param    nstart	         the number of random sets chosen for kmean clustering initialization.
#' @param    nproc           number of cores to distribute the function over.
#' @param    progress        if \code{TRUE} print a progress bar.
#' @param    out_name        if not \code{NULL} save cloned data.frame and a summary of cloning
#'                           performance. \code{out_name} string is used as the prefix of the successfully
#'                           processed output files.
#' @param    out_dir         specify to change the output directory. The input file
#'                           directory is used if this is not specified while \code{out_name} is specified.
#'
#' @return
#' Returns a modified \code{db} data.frame with clone identifiers in the \code{CLONE} column.
#' if \code{out_name} is not \code{NULL}, it will save the modified \code{db} and a summary
#' of cloning performance in the current directory or the specified \code{out_dir}.
#'
#' @details
#' An unsupervised pipeline to identify B cell clones from adaptive immune receptor
#' repertoire sequencing (AIRR-Seq) datasets. This method is based on spectral clustering
#' of the junction sequences of B cell receptors (BCRs, also referred to as Immunoglobulins,
#' (Igs)) that share the same V gene, J gene and junction length. It uses an adaptive
#' threshold that analyzes sequences in a local neighborhood.
#'
#' @seealso
#' To assess the performance of clonal assignment process check \code{analyzeClones}.
#'
#' @examples
#' # clone data using defineClonesScoper function
#' db <- defineClonesScoper(ExampleDb, junction = "JUNCTION", v_call = "V_CALL",
#'                          j_call = "J_CALL", first = TRUE)
#' @export
defineClonesScoper <- function(db,
                               junction = "JUNCTION",
                               v_call = "V_CALL",
                               j_call = "J_CALL",
                               first = FALSE,
                               cdr3 = FALSE,
                               mod3 = FALSE,
                               iter_max = 1000,
                               nstart = 25,
                               nproc = 1,
                               progress = FALSE,
                               out_name = NULL,
                               out_dir = ".") {


    # Initial checks
    if (!is.data.frame(db)) {
        stop("Must submit a data frame")
    }

    # number of enteries
    rec_count <- nrow(db)

    # temp columns
    temp_cols <- c("VJ_GROUP", "L", "JUNC_temp", "NUMBRE_OF_N", "NUMBRE_OF_DOT", "ID","CLONE_temp")

    # check for invalid columns
    invalid_cols <- c("CLONE", temp_cols)
    if (any(invalid_cols %in% colnames(db))) {
        stop("Column(s) '", paste(invalid_cols[invalid_cols %in% colnames(db)], collapse = "', '"), "' already exist.",
             "\n Invalid column names are: '", paste(invalid_cols, collapse = "', '"), "'.")
    }

    # Check for valid columns
    columns <- c(junction, v_call, j_call) #, fields
    columns <- columns[!is.null(columns)]
    check <- checkColumns(db, columns)
    if (check != TRUE) { stop(check) }

    # Check for invalid characters
    valid_chars <- colnames(getDNAMatrix(gap=0))
    .validateSeq <- function(x) { all(unique(strsplit(x, "")[[1]]) %in% valid_chars) }
    valid_seq <- sapply(db[[junction]], .validateSeq)
    not_valid_seq <- which(!valid_seq)
    if (length(not_valid_seq) > 0) {
        stop("Invalid sequence characters in the ", junction, " column. ",
                length(not_valid_seq)," sequence(s) found.", "\n Valid characters are: '",  valid_chars, "'")
    }

    # chaeck outputs
    if (!is.null(out_name)) {
        if (is.null(out_dir)) stop("out_dir must be specified.")
        if (!dir.exists(out_dir)) stop("out_dir", out_dir, "does not exist.")
    }

    # neighborhood <- match.arg(neighborhood)
    neighborhood <- "infer"
    similarity <- NULL
    # if (is.null(neighborhood)) {
    #     stop(" 'neighborhood' must be specified.")
    # }
    # if (neighborhood == "assign" & is.null(similarity)) stop("similarity needs to be assigned.")

    # add temporary junction column
    db$JUNC_temp <- db[[junction]]

    # add junction length column
    db$L <- stri_length(db$JUNC_temp)

    # check for mod3
    # filter mod 3 junction lengths
    if (mod3) {
        nIint <- nrow(db)
        db <- filter_(db, ~L%%3 == 0)
        nFin <- nrow(db)
        cat("MOD3 FILTER> ", "\n")
        cat(nIint-nFin, " invalid junction length(s) not mod3 in the ", junction, " column. ", "\n")
    }

    # check for cdr3
    # filter junctions with length >6
    if (cdr3) {
        nIint <- nrow(db)
        db <- filter_(db, ~L > 6)
        nFin <- nrow(db)
        cat("CDR3 FILTER> ", "\n")
        cat(nIint-nFin, " invalid junction length(s) < 7 in the ", junction, " column. ", "\n")
        # add cdr3 temp column
        db$JUNC_temp <- substr(db$JUNC_temp, 4, db$L-3)
        # update cdr3 length column
        db$L <- stri_length(db$JUNC_temp)
    }

    # Parse V and J columns to get gene
    db <- groupGenes(db,
                     v_call=v_call,
                     j_call=j_call,
                     first=first)

    # groups to use
    groupBy <- c("VJ_GROUP", "L")
    # if (!is.null(fields)) {
    #     groupBy <- append(groupBy, fields)
    # }

    # unique groups
    uniqueGroups <- data.frame(unique(db[, groupBy]))
    colnames(uniqueGroups) <- groupBy
    rownames(uniqueGroups) <- NULL
    n_groups <- nrow(uniqueGroups)
    groupBy_length <- length(groupBy)

    # Create cluster of nproc size and export namespaces
    if(nproc == 1) {
        # If needed to run on a single core/cpu then, register DoSEQ
        # (needed for 'foreach' in non-parallel mode)
        registerDoSEQ()
    } else if( nproc > 1 ) {
        cluster <- parallel::makeCluster(nproc, type="PSOCK")
        registerDoParallel(cluster)
    } else {
        stop('Nproc must be positive.')
    }

    # check the progressbar
    if (progress) {
        cat("CLONING> ", "\n")
        pb <- progressBar(n_groups)
    }

    # set seed for reproducibility
    set.seed(12345)

    # cloning
    i <- NULL
    db_cloned <- foreach(i=1:n_groups,
                         .combine="rbind",
                         .export=c("filterFunction", "spectralClustering", "krnlMtxGenerator",
                                   "makeAffinity", "kMeanClustering", "laplacian_mtx", "floor_dec", neighborhood),
                         #.packages = c("dplyr"),
                         .errorhandling='stop') %dopar% {
                             ft <- uniqueGroups[i, ]
                             filterBy <- c()
                             for (j in 1:groupBy_length) filterBy <- c(filterBy, as.character(ft[[groupBy[j]]]))
                             db_group <- db %>%
                                 dplyr::filter_(filterFunction(groupBy, filterBy))
                             db_group$ID <- db_group %>%
                                 dplyr::group_by_(.dots = "JUNC_temp") %>%
                                 dplyr::group_indices()
                             if (length(unique(db_group$JUNC_temp)) == 1) {
                                 CLONE <- data.frame(CLONE=as.vector(paste(i, rep(1, times=nrow(db_group)), sep="_")))
                             } else {
                                 idCluster <- spectralClustering(entrySeq=db_group$JUNC_temp,
                                                                 id=db_group$ID,
                                                                 similarity=similarity,
                                                                 neighborhood=neighborhood,
                                                                 iter_max=iter_max,
                                                                 nstart=nstart)
                                 CLONE <- data.frame(CLONE=as.vector(paste(i, idCluster, sep="_")))
                             }
                             # Update progress
                             if (progress) { pb$tick() }
                             # return result from each proc
                             return(bind_cols(db_group, CLONE))
                         }

    # Stop the cluster
    if (nproc > 1) { parallel::stopCluster(cluster) }

    db_cloned$CLONE_temp <- db_cloned %>%
        dplyr::group_by(CLONE) %>%
        dplyr::group_indices()
    db_cloned$CLONE <- db_cloned$CLONE_temp
    db_cloned <- db_cloned[order(db_cloned$CLONE), ]
    db_cloned$CLONE <- as.character(db_cloned$CLONE)

    # print into saveSummary
    if (!is.null(out_name)) {
        if (progress) cat("SAVING SUMMARY> ", "\n")
        db_summ <- db_cloned %>%
            dplyr::group_by_at(vars(one_of(groupBy))) %>%
            dplyr::summarise_(NUMBER_OF_SEQUENCE = interp(~length(x), x=as.name("JUNC_temp")),
                              NUMBER_OF_UNIQUE_SEQUENCE = interp(~length(unique(x)), x=as.name("JUNC_temp")),
                              NUMBER_OF_CLONE = interp(~length(unique(x)), x=as.name("CLONE")),
                              CLONE = interp(~paste(unique(x), collapse = ","), x=as.name("CLONE")))
        colnames(db_summ)[colnames(db_summ) == "L"] <- "JUNCTION_LENGTH"
        fileName <- paste(out_name, "summary-pass.tsv", sep="_")
        writeChangeoDb(db_summ, file = file.path(out_dir, fileName))
    }

    # remove extra columns
    db_cloned <- db_cloned[, !(names(db_cloned) %in% temp_cols)]

    # save cloned db
    if (!is.null(out_name)) {
        if (progress) cat("SAVING DB> ", "\n")
        fileName <- paste(out_name, "clone-pass.tsv", sep="_")
        writeChangeoDb(db_cloned, file = file.path(out_dir, fileName))
        }

    clone_count <- length(unique(db_cloned$CLONE))
    pass_count <- nrow(db_cloned)
    fail_count <- rec_count - pass_count

    cat(paste("CLONES= ", clone_count), "\n", sep="")
    cat(paste("RECORDS= ", rec_count), "\n", sep="")
    cat(paste("PASS= ", pass_count), "\n", sep="")
    cat(paste("FAIL= ", fail_count), "\n", sep="")
    return(db_cloned)
}



#### Clonal analysis ####

#' Clonal assignment analysis
#'
#' The \code{analyzeClones} function performs a series of analysis to assess the performance of
#' \code{defineClonesScoper} function.
#'
#' @param    db              data.frame with Change-O style columns containing sequence data.
#' @param    junction        name of the column containing nucleotide sequences to compare.
#'                           Also used to determine sequence length for grouping.
#' @param    v_call          name of the column containing the V-segment allele calls.
#' @param    j_call          name of the column containing the J-segment allele calls.
#' @param    clone	         name of the data column containing clone identifiers.
#' @param    first           if \code{TRUE} only the first call of the gene assignments
#'                           is used. if \code{FALSE} the union of ambiguous gene
#'                           assignments is used to group all sequences with any
#'                           overlapping gene calls.
#' @param    cdr3            if \code{TRUE} remove 3 nts from both ends of \code{junction}
#'                           (converts IMGT junction to CDR3 region).
#' @param    nproc           number of cores to distribute the function over.
#' @param    progress        if \code{TRUE} print a progress bar.
#'
#' @note
#' Arguments \code{first} and \code{cdr3} must match the corresponding arguments
#' used in the \link{defineClonesScoper} function.
#'
#' @return
#' Returns a \link{ClonalAnalysis} object.
#'
#'@examples
#' # Clonal assignment analysis
#' results <- analyzeClones(ClonedExampleDb, junction = "JUNCTION", v_call = "V_CALL",
#'                           j_call = "J_CALL", clone = "CLONE", first = TRUE)
#' # print threshold (a numeric)
#' results@threshold
#'
#' # get inter and intra conal distances (a data.frame)
#' df <- results@inter_intra[[1]]
#'
#' # density plot of inter versus intra clonal distances  (a ggplot).
#' results@plot_inter_intra
#'
#' # get the neighborhoods used in spectral clustering (a numeric vector).
#' ngs <- results@neighborhoods
#'
#' # plot histogram of neighborhoods (a ggplot).
#' results@plot_neighborhoods
#' @export
analyzeClones <- function(db,
                           junction = "JUNCTION",
                           v_call = "V_CALL",
                           j_call = "J_CALL",
                           clone = "CLONE",
                           first = FALSE,
                           cdr3 = FALSE,
                           nproc = 1,
                           progress = FALSE) {

    # Initial checks
    # Check for valid columns
    columns <- c(junction, v_call, j_call, clone)
    columns <- columns[!is.null(columns)]
    check <- checkColumns(db, columns)
    if (check != TRUE) { stop(check) }

    # add junction temp column
    db$JUNC_temp <- db[[junction]]

    # add junction length column
    db$L <- stri_length(db$JUNC_temp)

    # check for cdr3
    if (cdr3) {
        # add cdr3 temp column
        db$JUNC_temp <- substr(db$JUNC_temp, 4, db$L-3)
        # update cdr3 length column
        db$L <- stri_length(db$JUNC_temp)
    }

    # Parse V and J columns to get gene
    db <- groupGenes(db,
                     v_call=v_call,
                     j_call=j_call,
                     first=first)

    # calculate inter and intra distances
    results <- calculateInterVsIntra(db,
                                     junction = junction,
                                     v_call = v_call,
                                     j_call = j_call,
                                     clone = clone,
                                     first = first,
                                     cdr3 = cdr3,
                                     nproc = nproc,
                                     progress = progress)

    # revoke the results
    df <- results@inter_intra[[1]]
    threshold <- results@threshold
    interVsIntra <- list()
    interVsIntra[[length(interVsIntra)+1]] <- df

    # plot inter and intra distances
    p1 <- list()
    p <- plotInterVsIntra(df,
                           threshold=threshold)
    p1[[length(p1)+1]] <- p

    # calculate neighborhoods
    neighborhood <- "infer"
    similarity <- NULL
    neighborhoods <- calculateNeighborhoods(db,
                                            junction = junction,
                                            v_call = v_call,
                                            j_call = j_call,
                                            similarity = similarity,
                                            neighborhood = neighborhood,
                                            first = first,
                                            cdr3 = cdr3,
                                            progress = progress)

    # plot neighborhoods
    p2 <- list()
    p <- plotNeighborhoods(neighborhoods,
                           threshold = threshold)
    p2[[length(p2)+1]] <- p

    # return all
    analysis <- new("ClonalAnalysis",
                    threshold=threshold,
                    inter_intra=interVsIntra,
                    plot_inter_intra=p1,
                    neighborhoods=neighborhoods,
                    plot_neighborhoods=p2)

    return(analysis)
}
