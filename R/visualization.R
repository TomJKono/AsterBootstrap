# Functions related to visualizing data and models.

#' Visualize an Aster graphical model
#' 
#' @param mod_spec
#' The model specification from the build_graphical_model() function.
#' @param linear.layout
#' Boolean flag for whether graphical model should be plotted in a straight
#' line or not. Default is TRUE: do plot the graph in a straight line.
#' @return
#' A 'dagitty' object cooresponding to the directed acyclic graph (DAG) that
#' describes the Aster graphical model. This can be plotted with custom
#' plotting functions in the 'dagitty' or 'ggdag' packages.
#' @details
#' This function takes the data in the '$initial', '$vars' and '$pred' elements
#' of the model specification list to generate a 'dagitty' object that
#' represents the DAG for the Aster graphical model. The function then
#' calculates coordinates for plotting and draws the plot.
#' 
#' Currently, only "unbranched" graphical models are supported.
#' 
#' @examples
#' # Build a model specification
#' model_spec <- build_graphical_model(data, initial="Initial",
#'                                     vars=c("Germinated", "Fruits", "Seeds"),
#'                                     pred=c("Initial", "Germinated", "Fruits"),
#'                                     fams=c(fam_bernoulli, fam_0poi, fam_poi))
#' # Draw a plot to interactive R window
#' dag_layout <- draw_graphical_model(model_spec)
draw_graphical_model <- function(mod_spec, linear.layout=TRUE) {
    # Get the nodes from the model specification
    nodes <- mod_spec$vars
    # Put them in order based on the predecessor->successor notation
    nodes_ordered <- c(mod_spec$initial, nodes[mod_spec$pred+1])
    # Build the string-like specification of the DAG, for dagitty
    dag_nodes <- c('dag {')
    for(x in seq(2, length(nodes_ordered))) {
        prev <- x - 1
        curr <- x
        node_desc <- paste(nodes_ordered[prev], "->", nodes_ordered[curr], sep="")
        dag_nodes <- c(dag_nodes, node_desc)
    }
    dag_nodes <- c(dag_nodes, '}')
    dag_spec <- paste(dag_nodes, sep=" ", collapse=" ")
    dag <- dagitty::dagitty(dag_spec)
    # Set the coordinates, if 'linear.layout' is TRUE
    if(linear.layout) {
        x_coord <- seq(0, length(nodes_ordered))
        names(x_coord) <- nodes_ordered
        y_coord <- rep(0, length(nodes_ordered))
        names(y_coord) <- nodes_ordered
        dagitty::coordinates(dag) <- list(x=x_coord, y=y_coord)
    } else {
        dag <- dagitty::graphLayout(dag)
    }
    # Draw the plot!
    plot(dag)
    return(dag)
}
