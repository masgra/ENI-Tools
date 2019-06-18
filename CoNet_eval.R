library(igraph)


# shoul biom files be extracted again  
refresh.files = F


nc.path <- readLines( file("localpath.txt","r") ,n=1)

n.nodes = 999

# load node degree distribution from CoNet export (node table export)
d <- list()
d$RC$Node.degree <- read.csv(paste(nc.path,"Results/Network-Analysis/RC-002_node.csv", sep=""), header = T, sep = ",")$Degree
d$RD$Node.degree <- read.csv(paste(nc.path,"Results/Network-Analysis/RD-002_node.csv", sep=""), header = T, sep = ",")$Degree

# creade edge distribution of random graph model
d$RC.model$Node.degree <- igraph::degree(igraph::sample_gnm(n=n.nodes,m = sum(d$RC$Node.degree)/2, directed = F, loops = F))
d$RD.model$Node.degree <- igraph::degree(igraph::sample_gnm(n=n.nodes,m = sum(d$RC$Node.degree)/2, directed = F, loops = F))

# calculate degree distribution, omit zeros
data <- lapply(d, function(x) {res <- data.frame(Nodes=as.numeric(names(table(x))), Degree=c(table(x)) )
                              res[res==0] <- NA
                              return(res[complete.cases(res),])})


# fit the power law distribution P(x) = C* x^-alpha
fit_power_law = function(data) {
  reg = lm(log(data$Degree) ~ log(data$Node))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  return(list(fun=power.law.fit, "alpha"=alpha,"R_square"=R.square)) 
}

# calculate fit for each distribution, including alpha and R-square
plf <- lapply(data, function(x) fit_power_law(x))



library("ggplot2")
# plot desity distributions 
ggplot()  +  
  geom_point(data = data$RC.model, aes(x=Nodes, y=Degree, color="red"),size=2) +
  geom_point(data = data$RC, aes(x=Nodes, y=Degree, color="green"), size=2 ) +
  geom_line(data = data.frame(Nodes=seq(1,max(data$RC$Nodes),by=.01), 
                              Degree=plf$RC$fun(seq(1,max(data$RC$Nodes),by=.01))), aes(x=Nodes,y=Degree, color="black"))+
  labs(title =bquote('Node degree Distribution of'~R[C]), x="Number of nodes", y="Degree" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5), legend.position=c(.8,.9)) +
  scale_y_continuous(trans='log2', limits = c(1,512))+ 
  scale_x_continuous(trans='log2')+ scale_colour_manual(name = "", 
                                                        labels=c(bquote(R[C]~'power-law'), bquote(~R[C]), 'Random Network'),
                                                        values=c("#FF3366", "#FF3366", "#333333"),
                                                        guide = guide_legend(override.aes = list(
                                                            linetype = c("solid", rep("blank", 2)),
                                                            shape = c(NA, rep(16, 2)))))


ggplot(data = data$RD.model, aes(x=Nodes, y=Degree, color="red"))  +  
  geom_point(size=2) +
  geom_point(data = data$RD, aes(x=Nodes, y=Degree, color="green"), size=2 )+ 
  geom_line(data = data.frame(Nodes=seq(1,max(data$RD$Nodes),by=.01), 
                              Degree=plf$RD$fun(seq(1,max(data$RD$Nodes),by=.01))), aes(x=Nodes,y=Degree,color="black"))+
  labs(title =bquote('Node degree Distribution of'~R[D]), x="Number of nodes", y="Degree" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5), legend.position=c(.8,.9)) +
  scale_y_continuous(trans='log2', limits = c(1,512))+ 
  scale_x_continuous(trans='log2')+ scale_colour_manual(name = "", 
                                                     labels=c(bquote(R[D]~'power-law'), bquote(~R[D]), 'Random Network'),
                                                     values=c("#0066CC", "#0066CC", "#333333"),
                                                     guide = guide_legend(override.aes = list(
                                                       linetype = c("solid", rep("blank", 2)),
                                                       shape = c(NA, rep(16, 2)))))



#library(poweRlaw)
## power law fiting for cumulative density plot
#data_pl <-list(); 
#data_pl$model <- displ$new(data$model$Degree[data$model$Degree>0])
#data_pl$RC <- displ$new(data$RC$Degree[data$RC$Degree>0])
#data_pl$RD <- displ$new(data$RD$Degree[data$RD$Degree>0])
## calculate characteristic parameters
#est_pl <- lapply(data_pl, function(x) estimate_xmin(x))
## p-value calculation (should be < 0.05)
#bs <- lapply(data_pl, function(x) bootstrap_p(x,no_of_sims=1000, threads=2, seed = 42))
#sapply(bs, function(x) (x$p)) #print p-values  
#
## fit power law model to data
#mapply(function(x,y)(x$setXmin(y$xmin)),data_pl,est_pl) 
#mapply(function(x,y)(x$setPars(y$pars)),data_pl,est_pl)
#plot.data <- lapply(data_pl, function(x) plot(data_pl$RC, draw = F))
#fit.data <- lapply(data_pl, function(x) lines(data_pl$RC, draw = F))
#
#
##cummulative density plot with power law fit
#ggplot() + 
#  geom_line(data=fit.data$RC, aes(x=x, y=y), colour="red") + 
#  geom_line(data=fit.data$RD, aes(x=x, y=y), colour="red") + 
#  geom_line(data=fit.data$model, aes(x=x, y=y), colour="red") + 
#  geom_point(data=plot.data$model, aes(x=x, y=y), color="blue") + 
#  labs(x="k", y="CDF") + theme_bw() 



## ------------- Wilcoxon rank sum test of Edge Betweenness between intersect and union (Mann-Whitney tes) --------------

# load node degree distribution from CoNet export (node table export)
ds <- list()
ds$RD$EB <- read.csv(paste(nc.path,"Results/Network-Analysis/RD-002_edge.csv", sep=""), header = T, sep = ",")$EdgeBetweenness
ds$RC$EB <- read.csv(paste(nc.path,"Results/Network-Analysis/RC-002_edge.csv", sep=""), header = T, sep = ",")$EdgeBetweenness
ds$RD$EB.sub <- read.csv(paste(nc.path,"Results/Network-Analysis/RD-002_intersect-edge.csv", sep=""), header = T, sep = ",")$EdgeBetweenness
ds$RC$EB.sub <- read.csv(paste(nc.path,"Results/Network-Analysis/RC-002_intersect-edge.csv", sep=""), header = T, sep = ",")$EdgeBetweenness

# H0: distributions are not shifted? --> accept because p>0.05 (0.1215) 
wilcox.test(ds$RC$EB, ds$RC$EB.sub, correct = F)
# H0: EB.sub distributions is not greater than EB? -> accepted. p > 0.05 (0.9393)
wilcox.test(ds$RC$EB.sub, ds$RC$EB, correct = F,alternative = c("greater"))

# H0: distributions are not shifted? --> rejected because p<0.05 (~0)
wilcox.test(ds$RD$EB, ds$RD$EB.sub, correct = F)
# H0: EB.sub distributions is not greater than EB? -> accepted. p > 0.05 (1.0)
wilcox.test(ds$RD$EB.sub, ds$RD$EB, correct = F, alternative = c("greater"))

#RC: distribution of Edge Betweenness from sub set is similar to overall distribution. 
#RD: distribution of Edge Betweenness from sub set is significantly different to overall distribution. 
#    However, the shift is not significant upward (it is significant towards lower betweenness).
#    --> some metrics might generate edges with very high betweenness which are not supported by other metrices. 
#        (the correllation metrices)





