#' Functional annotation for multiple gene sets. Parallel version.
#'
#' Conducts functional annotation using Fisher exact test
#' for multiple gene sets, returns list with resulting annotations
#'
#' @param gene_amount_border minimal number of genes annotated to a term (1 by default)
#' @param p_adjust_method method for multiple testing correction (to see all possible methods print: p.adjust.methods)
#'                        Benjamini-Hochberg by default
#' @param background_genes vector contains background set of genes
#' @param threads number of threads
#' @param gene_sets list contains quantile names (e.g. "1","2","1-6") as keys and corresponding gene sets as values
#' @param annot_data list object contains annotation data ready for functional annnotation procedure,
#'           ontology relation data, GO terms full names, GO terms namespaces and metadata for annotation
#'
#' @return list with filenames as keys and annotaton data frames as values
#' @export
#'
#' @examples
#'     data(up_genes, annot_data, package = "fsgor")
#'     gene_sets <- div_genes_to_quantiles(up_genes, 6)
#'     par_multiple_sets_annot(gene_sets, annot_data, threads = 6,  gene_amount_border = 10)
#'
par_multiple_sets_annot <-
  function(gene_sets,
           annot_data,
           threads,
           gene_amount_border = 1,
           background_genes = NULL,
           p_adjust_method = "BH") {

    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("parallel package needed for this function to work. Please install it.",
           call. = FALSE)
    }

    annotlist <- annot_data$annot_data
    realnames <- annot_data$real_names
    namespaces <- annot_data$namespaces

    if (!is.null(background_genes)){
      annotlist <- fsgor::intersect_with_background(annotlist, background_genes)
    }

    cl <- parallel::makeCluster(threads)

    parallel::clusterExport(
      cl = cl,
      varlist = c(
        "func_annot_test_internal",
        "gene_sets",
        "gene_amount_border",
        "p_adjust_method"
      ),
      envir = environment()
    )

    outlist <- parallel::parLapply(cl, gene_sets, function(inputgenes) {

      # run functonal annotation
      return(
        func_annot_test_internal(
          annotlist,
          realnames,
          namespaces,
          inputgenes,
          gene_amount_border = gene_amount_border,
          p_adjust_method = p_adjust_method
        )
      )
    })
    names(outlist) <- names(gene_sets)
    parallel::stopCluster(cl)
    return(outlist)
  }
