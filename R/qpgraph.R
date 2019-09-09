library(DiagrammeR)
library(magrittr)

extract_feature <- function(lines, what) {
    lines %>%
        .[grepl(what, .)] %>%
        gsub(paste0(what, " +"), "", .) %>%
        gsub(" +", " ", .) %>%
        gsub(" $", "", .) %>%
        strsplit(" ")
}

parse_qpgraph <- function(path) {
    graph <- DiagrammeR::create_graph(attr_theme = NULL)

    lines <- readLines(path) %>% .[!grepl("^#", .)]
    
    ## for whatever reason, there are two possible versions of
    ## a qpGraph graph file
    nodes <- extract_feature(lines, "label")
    branches <- extract_feature(lines, "edge")
    admixtures <- extract_feature(lines, "admix")

    for (node in nodes)
        graph <- DiagrammeR::add_node(graph,
                                      label = node[1],
                                      type = node[2])

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

        style <- DiagrammeR::edge_aes(color = "black")

        graph <- DiagrammeR::add_edge(graph, from, to, edge_aes = style)
    }

    # add admixture edges
    for (admixture in admixtures) {
        source1 <- all_nodes[all_nodes$type == admixture[2], "label"]
        source2 <- all_nodes[all_nodes$type == admixture[3], "label"]
        target <- all_nodes[all_nodes$type == admixture[1], "label"]

        prop1 <- DiagrammeR::edge_data(value = paste0(admixture[4], "%"))
        prop2 <- DiagrammeR::edge_data(value = paste0(admixture[5], "%"))

        style <- DiagrammeR::edge_aes(style = "dotted")

        graph <- DiagrammeR::add_edge(graph, source1, target, edge_aes = style, edge_data = prop1)
        graph <- DiagrammeR::add_edge(graph, source2, target, edge_aes = style, edge_data = prop2)
    }

    graph %>%
        set_edge_attr_to_display(attr = value)
}

g <- parse_qpgraph("~/local/AdmixTools-5.1/examples.qpGraph/gr1x")

render_graph(g)


g %>% generate_dot %>% writeLines(con = "/tmp/test.txt")

grViz("/tmp/test.txt")

grViz("~/local/AdmixTools-5.1/examples.qpGraph/sim1:gr1x.dot")



admixture_graph <- function(snps) {
    data <- list(eigenstrat = snps, graph = DiagrammeR::create_graph())
    class(data) <- "admixture_graph"
    data
}

        
