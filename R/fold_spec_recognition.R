#' Fold-specific GO terms recognition
#'
#' Recognizes fold-specific GO for each fold change interval versus differentialy expressed genes annotation,
#' returns list contains dataframe with fold-specific GO terms data and not fold-fold-specific GO terms data
#' accessible by fs and nfs keys correspondingly
#'
#' @param p_adjust_method method name for correction on multiple correction, type "p.adjust.methods" in R Console
#'                        to obtain all possiple values
#' @param fdr2step value for fdr step 2
#' @param fdr3step value for fdr step 3
#' @param listoftables list contains GO output tables as dataframes and interval names as leys
#' @param wholeintname key for whole interval table specified in listoftables
#' @param listofcolnames list of column names, must contain the following keys: "GO_id" - GO term identifier
#'                                                                              "name" - GO term name
#'                                                                              "namespace" - namespace of GO term
#'                                                                              "qit" - number of genes annotated to a term in query list
#'                                                                              "qtot" - number of genes in a query list
#'                                                                              "padj" - adjusted p-values, if you don't use any corrections on multiple testing assign the name of raw p-values column to this key
#'                                                                              (column names of internal annotation used as default value)
#' @param fisher_alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". You can specify just the initial letter. Only used in the 2 by 2 case.
#'
#' @return list contains dataframe with fold-specific and not fold-fold-specific GO terms data
#'
#' @export
#'
#' @examples
#'     data(up_annot, package = "fsgor")
#'     recognize_fs_terms(up_annot, "1-6", fsgor::colnameslist, "BY", 1, 0.05)
recognize_fs_terms <-
  function (listoftables,
            wholeintname,
            listofcolnames = fsgor::colnameslist,
            p_adjust_method,
            fdr2step,
            fdr3step,
            fisher_alternative = "greater") {

    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("tidyr package needed for this function to work. Please install it.",
           call. = FALSE)
    }

    # extract dataframe for whole interval annotation
    # and cut it by FDR step 2 threshold
    df <- listoftables[[wholeintname]]
    if (fdr2step < 1) {
      df <- subset(df, df[[listofcolnames$padj]] < fdr2step)
    }

    # leave only GO_id, qit and qtot columns
    # and rename qit and qtot to wqit and wqtot where w means whole interval
    df <-
      df[, c(listofcolnames$GO_id,
             listofcolnames$qit,
             listofcolnames$qtot)]
    names(df)[names(df) == listofcolnames$qit] <- "wqit"
    names(df)[names(df) == listofcolnames$qtot] <- "wqtot"

    # create empty dataframe for merging all dataframes by GO_id
    resdf <- data.frame()

    # extract interval names form list of dataframes names
    # and delete the name of whole interval
    filenames <- names(listoftables)
    filenames <- filenames[!filenames %in% wholeintname]

    # extract dataframe for each interval except whole
    # add column with interval name
    # merge data by GO_id's from whole interval data frame
    # bind all dataframes in resdf
    for (i in 1:length(filenames)) {
      data <- listoftables[[filenames[i]]]
      data$filename <- rep(filenames[i], nrow(data))
      data <- data[, c(
        listofcolnames$GO_id,
        listofcolnames$name,
        listofcolnames$namespace,
        listofcolnames$qit,
        listofcolnames$qtot,
        listofcolnames$qit_genes,
        "filename"
      )]
      resdf <-
        rbind(resdf, merge(data, df, by = listofcolnames$GO_id))
    }

    # create matrix with rows contains qit, qtot, wqit, wqtot for each GO term
    cutdf <-
      resdf[, c(listofcolnames$qit, listofcolnames$qtot, "wqit", "wqtot")]
    mat <- matrix(unlist(cutdf), nrow = nrow(cutdf))

    # calculate Fisher exact test for each GO term
    resdf$pvalues <-
      apply(mat, 1, function(x)
        stats::fisher.test(matrix(c(
          x[1], x[2] - x[1], x[3] - x[1], x[4] - x[2] - x[3] + x[1]
        ), nrow = 2),
        alternative = fisher_alternative)$p.value)


    genes_df <- resdf[, c(
      listofcolnames$GO_id,
      listofcolnames$qit_genes,
      "filename"
    )]

    # leave only GO_id, filename and pvalues rows in resulting dataframe
    # and convert it to wide format
    resdf <-
      resdf[, c(
        listofcolnames$GO_id,
        listofcolnames$name,
        listofcolnames$namespace,
        "filename",
        "pvalues"
      )]

    resdf <- tidyr::spread(resdf, "filename", "pvalues")
    genes_df <- tidyr::spread(genes_df, "filename", "qit_genes")

    # move GO_ids from first column to rownames
    rownames(resdf) <- resdf[, 1]
    resdf[, 1] <- NULL
    go_ids <- rownames(resdf)
    go_names <- resdf[, 1]
    go_namespaces <- resdf[, 2]
    resdf <- resdf[, c(-1, -2)]

    # find minimal p-values for each GO term and corresponding interval name
    minp <-
      apply(resdf, 1, function(x)
        c(min(x, na.rm = TRUE), colnames(resdf)[which.min(x)]))

    # from table with genes delete all intervals that are not in minp
    genes <-
      sapply(c(1:length(genes_df[, 1])), function(x)
        genes_df[x, which(colnames(genes_df) == minp[2, ][x])])

    # apply correction on multiple testing
    # and separate GO terms into fold-specific and not fold-specific groups
    # using fdr step 3 threshold
    padj <- stats::p.adjust(minp[1, ], method = p_adjust_method)
    fsinds <- which(as.numeric(padj) < fdr3step)
    nfsinds <- which(as.numeric(padj) >= fdr3step)

    # create dataframes for fold specific and
    # not fold specific terms containing GO_ids, minimal p adjusted values,
    # real names, interval names
    fs <- data.frame(
      ids = go_ids[fsinds],
      namespace = go_namespaces[fsinds],
      name = go_names[fsinds],
      padj = padj[fsinds],
      interval = minp[2, ][fsinds],
      genes = genes[fsinds],
      stringsAsFactors = FALSE
    )

    if (length(fsinds) != 0) {
      rownames(fs) <- c(1:length(fsinds))
    }

    nfs <- data.frame(
      ids = go_ids[nfsinds],
      namespace = go_namespaces[nfsinds],
      name = go_names[nfsinds],
      padj = padj[nfsinds],
      interval = minp[2, ][nfsinds],
      genes = genes[nfsinds],
      stringsAsFactors = FALSE
    )

    if (length(nfsinds) != 0) {
      rownames(nfs) <- c(1:length(nfsinds))
    }

    # create output list
    output <- list("fs" = fs, "nfs" = nfs)
    return(output)
  }
