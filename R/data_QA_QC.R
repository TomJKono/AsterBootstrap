# Functions related to data import and "cleaning" - raising errors and warnings
# if there are incompatible values in a data set.

#' Read a data sheet in preparation for an Aster analysis.
#' 
#' @param filename
#' Path to the data file to read (character string)
#' @param initial
#' Name of the column in the data file to use as the "Initial" node. This is
#' case-sensitive (character string)
#' @param sep
#' Separator or delimiter character used in the data file. Comma (",") for CSV
#' files and tab ("\t") for TSV files. Default is "," for CSV files.
#' @return
#' A data.frame containing the data in the specified file.
#' @examples
#' # Read a CSV file with an initial column called "Initial"
#' data <- read_aster_data(filename="aster_data.csv", initial="Initial", sep=",")
#' # Read a TSV file with an initial column called "Root"
#' data <- read_aster_data(filename="aster_data.tsv", initial="Root", sep="\t")
read_aster_data <- function(filename, initial, sep=",") {
    # Cast them to character strings
    filename <- as.character(filename)
    initial <- as.character(initial)
    # Read the data in a very general way
    dat <- read.delim(
        filename,
        sep=sep,
        quote=c("\"", "'"),
        comment="#")
    # Check the name of the initial node - it should exist in the data frame
    if(!(initial %in% colnames(dat))) {
        stop(paste(
            "Initial column '",
            initial,
            "' not found in data file! Data file contains: ",
            paste(colnames(dat), sep=", ", collapse=", "),
            sep=""),
        call.=TRUE)
    }
    # If the column for the initial node exists, then we return the data.frame
    return(dat)
}

#' Build a graphical model given input data
#' 
#' @param dat
#' The data.frame that will be analyzed with an Aster model. This should be the
#' output from read_aster_data().
#' @param initial
#' The name of the initial node in the graphical model. This should correspond
#' to one of the column names of `dat`.
#' @param vars
#' A vector of the names of the nodes in the graphical model, excluding the
#' initial node. These should correspond to column names of `dat`. See details.
#' @param pred
#' A vector of the names of the predecessor nodes for `vars`. The initial node
#' should be included in this vector and the names should correspond to the
#' column names of `dat`. See details.
#' @param fams
#' A vector of the exponential family distributions that should be used for the
#' edges of the graphical model. See details.
#' @param quiet
#' A boolean (TRUE/FALSE) for whether summary information about the graphical
#' model should be printed. Default TRUE.
#' @return
#' A list with three elements:
#'   $initial: the name of the "Initial" node of the graphical model
#' 
#'   $vars: the names of the nodes of the graphical model
#' 
#'   $pred: the numerical index of the predecessors for each node named in the
#'          $vars element.
#' 
#'   $fam: the numerical index of the exponential family distribution that will
#'         be used to model the edge from predecessor->successor
#' @details
#' This function is mostly a convenience function for specifying a graphical
#' model with a given data.frame. The `vars` argument should be a vector of the
#' names of the non-initial nodes in the graphical model. The `pred`
#' argument should be a vector of the names of the predecessor node for each
#' node given in `vars`, in the same order. The `fams` argument shoud give the
#' exponential family distribution that will be used to model the
#' predecessor->successor edges in the graphical model.
#' 
#' For example, consider a data.frame that has the following columns:
#' 
#'   Initial, Germinated, Fruits, Seeds
#' 
#' where `Initial` is the number of seeds planted in an experimental unit,
#' `Germinated` is a count of the number of seeds germinated, `Fruits` is a
#' count of the number of fruits set by all plants in the unit, and `Seeds` is
#' a count of the number of seeds collected from all fruits.
#' 
#' A graphical model that relates these variables could be:
#' 
#'         Bernoulli      0-Poisson    Poisson
#'   Initial ---> Germinated ---> Fruits ---> Seeds
#' 
#' Where arrows denote edges, and names above arrows denote exponential
#' families that describe the distribution of the successor node, conditional
#' on the predecessor node.
#' 
#' This function would be called as such:
#' 
#' model_spec <- build_graphical_model(data, initial="Initial",
#'                                     vars=c("Germinated", "Fruits", "Seeds"),
#'                                     pred=c("Initial", "Germinated", "Fruits"),
#'                                     fams=c(fam_bernoulli, fam_0poi, fam_poi))
#' 
#' This function will throw errors if the names of the nodes or the designation
#' of predecessor nodes are incompatible with a valid Aster graphical model.
build_graphical_model <- function(dat, initial, vars, pred, fams, quiet=TRUE) {
    # First, check that vars, pred, and fams are all the same length and
    # are greater than 2 (Aster models need at least two nodes)
    arg_lens <- unique(c(length(vars), length(pred), length(fams)))
    if(length(arg_lens) != 1) {
        stop(
            "The arguments for 'vars', 'pred', and 'fams' must all be the same length!",
            call.=TRUE)
    }
    if(arg_lens < 2) {
        stop("You must specify at least two nodes in your graphical model",
             call.=TRUE)
    }
    # Check that the initial node exists
    if(!(initial %in% colnames(dat))) {
        stop(paste(
            "Initial column '",
            initial,
            "' not found in data file! Data file contains: ",
            paste(colnames(dat), sep=", ", collapse=", "),
            sep=""),
        call.=TRUE)
    }
    # Check that all of the variables listed are found as columns in the
    # data.frame
    mismatch <- vars[!vars %in% colnames(dat)]
    if(length(mismatch) > 0) {
        stop(
            paste(
                "Some of the names in 'vars' were not found in the data: ",
                paste(mismatch, sep=", ", collapse=", "),
                sep=""),
            call.=TRUE)
    }
    # Check, the predecessor names, too, including the initial node, which is
    # not part of the "vars" variable
    mismatch <- pred[!pred %in% c(initial, vars)]
    if(length(mismatch) > 0) {
        stop(
            paste(
                "Some of the names in 'pred' were not found in 'vars': ",
                paste(mismatch, sep=", ", collapse=", "),
                sep=""),
            call.=FALSE)
    }
    # Then, set up the predecessor->successor graph structure. Put the initial
    # node onto the vector of variables for this calcualtion. Subtract 1 to
    # account for the extra element (initial)
    pred_idx <- match(pred, c(initial, vars)) - 1
    # Eventually we will want to check for cycles. This is a bit complicated
    # for me right now, so I will come back to it later.
    # detect_cycles()
    # Print a table of the model to verify the specification, but only if the
    # user asks for it
    if(!quiet) {
        cat('Node order:\n', file=stderr())
        node_order <- c(initial, vars[pred_idx+1])
        cat(
            paste(seq(1, length(vars)+1), node_order, sep=": "),
            file=stderr(),
            sep="\n")
        cat('Edges:\n', file=stderr())
        model_relations <- data.frame(
            Predecessor=c(initial, vars[head(pred_idx, -1)+1]),
            Successor=vars,
            Family=fams)
        cat(capture.output(print(model_relations)), file=stderr(), sep="\n")
    }
    # Return a list with the model specifications
    ret <- list(
        initial=initial,
        vars=vars,
        pred=pred_idx,
        fam=fams)
    return(ret)
}

#' Check the validity of a data set for a given Aster model specification
check_data_validity <- function(dat, model_spec) {
    # Define a sub-function to check for error cases,  a 0 after a non-0. We
    # will also check the numbers for the family: in the case of Bernoulli,
    # Poisson, or 0-Poisson, there should only be positive integer values (or
    # 0)
    check_pred <- function(dat_row, dist_fam) {
        # Define the set of exponential families that take integer counts
        int_fams <- c(fam_bernoulli, fam_poi, fam_0poi, fam_trunc_poi)
        pred <- dat_row[1]
        succ <- dat_row[2]
        err_codes <- c(0)
        # First, check for missing values. This is the biggest problem - Aster
        # was developed the address a "missing data" problem!
        if(is.na(pred) | is.na(succ)) {
            err_codes <- c(err_codes, 1)
        } else {
            # Check that the successor is not >0 when the predecessor is 0
            if(succ > 0 & pred == 0) {
                err_codes <- c(err_codes, 1)
            }
            # If the relationship between these nodes is one of the exponential
            # families that requires integer data, we will check that these are
            # positive integers.
            if(dist_fam %in% int_fams) {
                if(!is.integer(pred) | !is.integer(succ)) {
                    err_codes <- c(err_codes, 2)
                }
                if(pred < 0 & !is.na(pred)) {
                    err_codes <- c(err_codes, 3)
                }
                if(succ < 0 & !is.na(succ)) {
                    err_codes <- c(err_codes, 4)
                }
            }
            # Now check the 0-truncated Poisson family. If an edge has a 0-poi
            # family, then it can only be 0 if and only if its predecessor is
            # 0.
            if(dist_fam == fam_0poi) {
                if(pred > 0 & succ == 0) {
                    err_codes <- c(err_codes, 5)
                }
            }
        }
        # Put the error codes into a string with semicolons separating them
        err_str <- paste(err_codes, sep=";", collapse=";")
        return(err_str)
    }
    # Get the order of the nodes in the graphical model
    nodenames <- model_spec$vars
    nodes_ordered <- c(model_spec$initial, nodenames[model_spec$pred+1])
    # For each node (successor) and its predecessor...
    #   This is going to be *rather inefficient*
    err_dat <- data.frame(
        DataRow=c(),
        Predecessor.Name=c(),
        Successor.Name=c(),
        Predecessor.Value=c(),
        Successor.Value=c(),
        Error.Codes=c())
    for(x in seq(2, length(nodes_ordered))) {
        pred_idx <- x-1
        succ_idx <- x
        fam <- model_spec$fam[pred_idx]
        test_dat <- cbind(
            dat[,nodes_ordered[pred_idx]],
            dat[,nodes_ordered[succ_idx]])
        errs <- apply(test_dat, 1, check_pred, fam)
        # Then, for each row tested, if there is an error, we will append a
        # row to the data.frame that is holding our error data
        for(ec_idx in seq_along(errs)) {
            if(errs[ec_idx] != "0") {
                # If there is an error then we will append a row to the
                # data.frame of errors
                ec <- gsub("0;", "", errs[ec_idx])
                err_row <- data.frame(
                    DataRow=ec_idx,
                    Predecessor.Name=nodes_ordered[pred_idx],
                    Successor.Name=nodes_ordered[succ_idx],
                    Predecessor.Value=dat[ec_idx, nodes_ordered[pred_idx]],
                    Successor.Value=dat[ec_idx, nodes_ordered[succ_idx]],
                    Error.Codes=ec)
                err_dat <- rbind(err_dat, err_row)
            }
        }
    }
    return(err_dat)
}

#' Reshape "wide" data to "long" data for Aster
#' 
make_long_data <- function(dat, model_spec, fitness_var) {
    # Check that the fitness variable exists as a column in the data
    if(!fitness_var %in% colnames(dat)) {
        stop(paste(
            "Fitness variable '",
            fitness_var,
            "' not found in data column names!",
            sep=""),
        call.=TRUE)
    }
    # Note that we want to use the "stats" reshape() rather than the function
    # in Hadley Wickham package.
    long_data <- stats::reshape(dat,
        varying=list(model_spec$vars),
        direction="long",
        timevar="varb",
        times=as.factor(model_spec$vars),
        v.names="resp")
    # And add the fitness surrogate variable into the long data
    long_data$fit <- as.numeric(long_data$varb == fitness_var)
    return(long_data)
}
