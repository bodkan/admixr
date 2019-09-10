library(DiagrammeR)
library(magrittr)

extract_feature <- function(lines, what) {
    lines %>%
        .[grepl(what, .)] %>%
        gsub(paste0("^.*", what, " +"), "", .) %>%
        gsub(" +", " ", .) %>%
        gsub(" $", "", .) %>%
        strsplit(" ")
}

parse_qpgraph <- function(path, integers = FALSE) {
    graph <- DiagrammeR::create_graph(attr_theme = NULL)
    
    lines <- readLines(path) %>% .[!grepl("^#", .)]

    ## for whatever reason, there are two possible versions of
    ## a qpGraph graph file - this concerns only the "labels"
    ## portion of the file (which has either two or three columns)
    ## or the fact that edges in the output qpGraph file contain
    ## drift values
    init_file <- !any(grepl("vertex", lines))

    nodes <- extract_feature(lines, "label") %>%
        unlist %>%
        matrix(ncol = ifelse(init_file, 2, 3), byrow = TRUE)
    nodes <- nodes[, if (init_file) 2:1 else 1:2] %>%
        as.data.frame %>%
        setNames(c("type", "label"))

    branches <- extract_feature(lines, "edge")
    admixtures <- extract_feature(lines, "admix")

    for (i in 1:nrow(nodes))
        graph <- DiagrammeR::add_node(graph,
                                      label = nodes[i, "label"],
                                      type = nodes[i, "type"])

    # first add missing inner nodes from the edge definition
    for (branch in branches) {
        node1 <- branch[2]
        node2 <- branch[3]
        if (!node1 %in% get_node_df(graph)$type)
            graph <- DiagrammeR::add_node(graph, node1, node1)
        if (!node2 %in% get_node_df(graph)$type)
            graph <- DiagrammeR::add_node(graph, node2, node2)
    }

    # collect information about all leaves and inner nodes
    all_nodes <- DiagrammeR::get_node_df(graph)

    # then walk through edges again and connect the nodes
    for (branch in branches) {
        from <- all_nodes[all_nodes$type == branch[2], "label"]
        to <- all_nodes[all_nodes$type == branch[3], "label"]

        drift <- DiagrammeR::edge_data(value = as.numeric(branch[4]))
        if (integers) drift$value <- round(1000 * drift$value)

        style <- DiagrammeR::edge_aes(style = "solid", color = "black")

        graph <- DiagrammeR::add_edge(graph, from, to,
                                      edge_aes = style, edge_data = drift)
    }

    # add admixture edges
    for (admixture in admixtures) {
        source1 <- all_nodes[all_nodes$type == admixture[2], "label"]
        source2 <- all_nodes[all_nodes$type == admixture[3], "label"]
        target <- all_nodes[all_nodes$type == admixture[1], "label"]

        prop1 <- as.numeric(admixture[4])
        prop2 <- as.numeric(admixture[5])
        if (integers) {
            prop1 <- round(100 * prop1)
            prop2 <- round(100 * prop2)
        }

        prop1 <- DiagrammeR::edge_data(value = paste0(prop1, "%"))
        prop2 <- DiagrammeR::edge_data(value = paste0(prop2, "%"))

        style <- DiagrammeR::edge_aes(style = "dotted", color = "black")

        graph <- graph %>%
            DiagrammeR::add_edge(source1, target,
                                 edge_aes = style, edge_data = prop1) %>%
            DiagrammeR::add_edge(source2, target,
                                 edge_aes = style, edge_data = prop2)

    }

    graph %>% DiagrammeR::set_edge_attr_to_display(attr = value)
}

g <- parse_qpgraph("~/local/AdmixTools-5.1/examples.qpGraph/gr1x")
render_graph(g)

g <- parse_qpgraph("~/local/AdmixTools-5.1/examples.qpGraph/sim1:gr1x.ggg")
render_graph(g)

g <- parse_qpgraph("~/local/AdmixTools-5.1/examples.qpGraph/sim1:gr1x.ggg", integers = T)
render_graph(g)

grViz("~/local/AdmixTools-5.1/examples.qpGraph/sim1:gr1x.dot")



admixture_graph <- function(snps) {
    data <- list(eigenstrat = snps, graph = DiagrammeR::create_graph())
    class(data) <- "admixture_graph"
    data
}

