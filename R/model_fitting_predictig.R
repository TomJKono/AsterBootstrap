# Functions related to fitting Aster models and generating point predictions
# from Aster models.

#' Fit a fixed-effects Aster model to a dataset
fixed_effects_aster <- function(dat, model_spec, blocking=NULL, formula=NULL,
    quiet=TRUE) {
    # If there are blocking factors specified, then coerce them to factor type
    if(!is.null(blocking)) {
        for(bf in blocking) {
            if(!bf %in% colnames(dat)) {
                stop(paste(
                    "Blocking factor '",
                    bf,
                    "' not found in the dataset!",
                    sep=""),
                call.=TRUE)
            }
            dat[[bf]] <- as.factor(dat[[bf]])
        }
    }

    # Build the formula. If the user has specified a custom formula, then use
    # that instead.
    if(is.null(formula)) {
        # The default formula will be a super basic one. It will use just the
        # response (resp), the graphical model variables (varb), the fitness
        # surrogate (fit), and whichever blocking factors were specified.
        if(is.null(blocking)) {
            mod <- resp ~ varb + fit
        } else {
            # We will only support additive "main effects" of blocking
            # factors. Interactions or more complicated linear combinations of
            # blocking factors will need a custom formula.
            block_string <- paste(blocking, collapse="+")
            mod_string <- paste("resp~varb+fit:(", block_string, ")", sep="")
            mod <- as.formula(mod_string)
        }
    } else {
        mod <- formula
    }

    # If the user asks for non-quiet output, we will print the formula
    if(!quiet) {
        cat("Formula:\n", file=stderr())
        cat(capture.output(print(mod)), file=stderr())
        cat("\n", file=stderr())
    }

    # Fit the model and return the model object
    aout <- aster::aster(
        formula=mod,
        pred=model_spec$pred,
        fam=model_spec$fam,
        varvar=varb,
        idvar=id,
        root=dat[[model_spec$initial]],
        data=dat)
    return(aout)
}


#' Fit a random-effects Aster model to a dataset
random_effects_aster <- function(dat, model_spec, reff=NULL, blocking=NULL,
    formula=NULL, quiet=TRUE) {
    # If there are blocking factors specified, then coerce them to factor type
    if(!is.null(blocking)) {
        for(bf in blocking) {
            if(!bf %in% colnames(dat)) {
                stop(paste(
                    "Blocking factor '",
                    bf,
                    "' not found in the dataset!",
                    sep=""),
                call.=TRUE)
            }
            dat[[bf]] <- as.factor(dat[[bf]])
        }
    }

    # Check the random effects specification. If it is NULL (none specified),
    # then we will throw an error. Otherwise, it should be a list of random
    # effects to add.
    if(is.null(reff) | !is.list(reff)) {
        stop("Random effects must be specified as a list! See help page.",
             call.=TRUE)
    }
    # This is a bit ugly, but hopefully the list of random effects is not so
    # large that this for() loop becomes a problem
    rand_effs <- list()
    for(re in names(reff)) {
        # Check that the random effect variable is a single character string.
        if(!is.character(reff[[re]]) & length(reff[[re]]) != 1) {
            stop(
                paste("Random effect '",
                      re,
                      "' is not a single character string.",
                      sep=""),
                call.=TRUE)
        }
        # Check that the random effects are part of the data to which the
        # model will be fit
        if(!reff[[re]] %in% colnames(dat)) {
            stop(paste(
                "Random effect '",
                re,
                "' not found in the dataset!",
                sep=""),
            call.=TRUE)
        }
        # Remove data rows where the random effects columns have missing data
        if(any(is.na(dat[[reff[[re]]]]))) {
            miss_rows <- which(is.na(dat[[reff[[re]]]]))
            warning(paste(
                "Input data has missing values for random effect '",
                re,
                "' on rows ",
                paste(miss_rows, collapse=", "),
                ". Removing these rows.",
                sep=""),
            immediate.=TRUE,
            call.=TRUE
            )
            dat <- dat[-miss_rows,]
        }
        # Finally if we get here, then make the random effects formulae
        revar <- reff[[re]]
        re_f_string <- paste("~0+fit:", reff[[re]], sep="")
        rand_effs[[re]] <- as.formula(re_f_string)
    }

    # Build the formula. If the user has specified a custom formula, then use
    # that instead.
    if(is.null(formula)) {
        # The default formula will be a super basic one. It will use just the
        # response (resp), the graphical model variables (varb), the fitness
        # surrogate (fit), and whichever blocking factors were specified.
        if(is.null(blocking)) {
            mod <- resp ~ varb + fit
        } else {
            # We will only support additive "main effects" of blocking
            # factors. Interactions or more complicated linear combinations of
            # blocking factors will need a custom formula.
            block_string <- paste(blocking, collapse="+")
            mod_string <- paste("resp~varb+fit:(", block_string, ")", sep="")
            mod <- as.formula(mod_string)
        }
    } else {
        mod <- formula
    }

    # If the user asked for non-quiet output, print what we will be fitting
    # to the terminal.
    if(!quiet) {
        cat("Random effects:\n", file=stderr())
        cat(capture.output(print(rand_effs)), file=stderr())
        cat("\n\n", file=stderr())
        cat("Formula:\n", file=stderr())
        cat(capture.output(print(mod)), file=stderr())
        cat("\n", file=stderr())
    }

    # Then, fit the model and return the model obj
    rout <- aster::reaster(
        fixed=mod,
        random=rand_effs,
        pred=model_spec$pred,
        fam=model_spec$fam,
        varvar=varb,
        idvar=id,
        root=dat[[model_spec$initial]],
        data=dat)
    return(rout)
}

# TODO:
#   See TR658 for Pearson residuals. We will likely want to include some
#   way to plot these.
