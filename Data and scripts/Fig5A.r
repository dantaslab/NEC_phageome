library(igraph)
library(readxl)
library(dplyr)
library(ggplot2)

##https://kateto.net/netscix2016.html
phage_host = read_excel("250421_Imp_phage_host.xlsx")
phage_host_early = filter(phage_host, NEC_onset=="Early")
links = as.matrix(with(phage_host_early, table(Contig, Phage_Host)))

phage_info = read_excel("250421_Imp_phage_host.xlsx", sheet = "621_ImpVir_info")
phage_info_early = filter(phage_info, NEC_onset=="Early")
nodes = phage_info_early
head(nodes)
head(links)

contig_case = phage_info$Contig[phage_info$NEC_status=="case"]
contig_control = phage_info$Contig[phage_info$NEC_status=="control"]

net2 = graph_from_incidence_matrix(links)
# V(net2)$color = c("orange", "steel blue")[V(net2)$type+1]

V(net2)$who = ifelse(V(net2)$name %in% contig_case, 'Case',
                      ifelse(V(net2)$name %in% contig_control, 'Control', 
                             'phage'))
V(net2)$color = c('#f28a5f', '#3291a4', '#d5e3ed')[factor(V(net2)$who)]
V(net2)$shape = c("circle", "rectangle")[V(net2)$type+1] # square
V(net2)$label = ""
V(net2)$label[V(net2)$type==T] = colnames(links) #nodes2$media[V(net2)$type==F] #colnames(links)
V(net2)$label.cex=0.5
V(net2)$label.font=2
# E(net2)$weight = 1

set.seed(127)
# pdf("250421_net_ImpVir_phage_bacteria_v1.pdf", 10, 10)
plot(net2, vertex.label.color="black", vertex.size=(V(net2)$type+0.12)*20, edge.color="light gray", 
     edge.curved=0, edge.width=1, vertex.size2 = 5) 
dev.off()

ceb = cluster_edge_betweenness(net2) 
clp = cluster_label_prop(net2)
V(net2)$shape = c("square", "circle", "rectangle")[factor(V(net2)$who)]
# library(RColorBrewer)
# cluster_colors = brewer.pal(length(ceb), "Set3")
# V(net2)$color = cluster_colors[membership(ceb)]
set.seed(127)
pdf("250421_net_ImpVir_phage_bacteria_v3.1.pdf", 10, 10)
plot(ceb, net2, vertex.label.color="black", vertex.size=(V(net2)$type+0.12)*20, edge.color="light gray", 
     edge.curved=0, edge.width=1, vertex.size2 = 7)
dev.off()

set.seed(138)
pdf("250426_net_ImpVir_phage_bacteria_early_v3.1.pdf", 10, 10)
plot(ceb, net2, vertex.label.color="black", vertex.size=(V(net2)$type+0.075)*20, edge.color="light gray", 
     edge.curved=0, edge.width=1, vertex.size2 = 7)
dev.off()











