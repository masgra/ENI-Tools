library(igraph)


# shoul biom files be extracted again  
refresh.files = F


nc.path <- readLines( file("localpath.txt","r") ,n=1)

n.nodes = 999

# load node degree export
data <- list()
data$RC <- read.csv(paste(nc.path,"Results/Network-Analysis/RC-002_node.csv", sep=""), header = T, sep = ",")
data$RD <- read.csv(paste(nc.path,"Results/Network-Analysis/RD-002_node.csv", sep=""), header = T, sep = ",")
# get degree distribution
data$RC.dist <- data.frame(Nodes=as.numeric(names(table(data$RC$Degree))), Degree= c(table(data$RC$Degree)))
data$RC.dist <- data$RC.dist[data$RC.dist$Degree>0,]
data$RD.dist <- data.frame(Nodes=as.numeric(names(table(data$RD$Degree))), Degree= c(table(data$RD$Degree)))
data$RD.dist <- data$RD.dist[data$RD.dist$Degree>0]

# creade edge distribution of random graph model
data$model$graph <- igraph::sample_gnm(n=n.nodes,m = 1117, directed = F, loops = F)
data$model$Degree <- igraph::degree(data$model$graph)
data$model.dist <- data.frame(Nodes=0:max(data$model$Degree), Degree=degree_distribution(data$model$graph)*n.nodes)
data$model.dist <- data$model.dist[data$model.dist$Degree>0,]



# fit the power law distribution P(x) = C* x^-alpha
fit_power_law = function(data) {
  # delete blank values
  data <- data[data$Degree>0,]
  data <- data[data$Nodes>0,]
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
plf <- list()
plf$RC <- fit_power_law(data$RC.dist)
plf$RD <- fit_power_law(data$RD.dist)
plf$model <- fit_power_law(data$model.dist)



library("ggplot2")
# plot desity distributions 
ggplot()  +  
  geom_point(data = data$model.dist, aes(x=Nodes, y=Degree),size=2) +
  #geom_line(data = data.frame(Nodes=seq(1,max(data$model.dist$Nodes),by=.01), 
  #                            Degree=plf$model$fun(seq(1,max(data$model.dist$Nodes),by=.01))), aes(x=Nodes,y=Degree))+
  geom_point(data = data$RC.dist, aes(x=Nodes, y=Degree) , color="blue", size=2 ) +
  geom_line(data = data.frame(Nodes=seq(1,max(data$RC.dist$Nodes),by=.01), 
                              Degree=plf$RC$fun(seq(1,max(data$RC.dist$Nodes),by=.01))), aes(x=Nodes,y=Degree), color="blue")+
  labs(title ='Node degree Distribution', x="Number of nodes", y="Degree" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(trans='log2')+scale_x_continuous(trans='log2') 




ggplot()  +  
  geom_point(data = data$model.dist, aes(x=Nodes, y=Degree),size=2) +
  #geom_line(data = data.frame(Nodes=seq(1,max(data$model.dist$Nodes),by=.01), 
  #                            Degree=plf$model$fun(seq(1,max(data$model.dist$Nodes),by=.01))), aes(x=Nodes,y=Degree))+
  geom_point(data = data$RD.dist, aes(x=Nodes, y=Degree) , color="green", size=2 )+ 
  geom_line(data = data.frame(Nodes=seq(1,max(data$RD.dist$Nodes),by=.01), 
                              Degree=plf$RD$fun(seq(1,max(data$RD.dist$Nodes),by=.01))), aes(x=Nodes,y=Degree), color="green")+
  labs(title ='Node degree Distribution', x="Number of nodes", y="Degree" )+
  theme(text = element_text(size=12), axis.text = element_text(size=10), plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(trans='log2')+scale_x_continuous(trans='log2') 
















library(poweRlaw)
# power law fiting for cumulative density plot
data_pl <-list(); 
data_pl$model <- displ$new(data$model$Degree[data$model$Degree>0])
data_pl$RC <- displ$new(data$RC$Degree[data$RC$Degree>0])
data_pl$RD <- displ$new(data$RD$Degree[data$RD$Degree>0])
# calculate characteristic parameters
est_pl <- lapply(data_pl, function(x) estimate_xmin(x))
# p-value calculation (should be < 0.05)
bs <- lapply(data_pl, function(x) bootstrap_p(x,no_of_sims=1000, threads=2, seed = 42))
sapply(bs, function(x) (x$p)) #print p-values  

# fit power law model to data
mapply(function(x,y)(x$setXmin(y$xmin)),data_pl,est_pl) 
mapply(function(x,y)(x$setPars(y$pars)),data_pl,est_pl)
plot.data <- lapply(data_pl, function(x) plot(data_pl$RC, draw = F))
fit.data <- lapply(data_pl, function(x) lines(data_pl$RC, draw = F))


#cummulative density plot with power law fit
ggplot() + 
  geom_line(data=fit.data$RC, aes(x=x, y=y), colour="red") + 
  geom_line(data=fit.data$RD, aes(x=x, y=y), colour="red") + 
  geom_line(data=fit.data$model, aes(x=x, y=y), colour="red") + 
  geom_point(data=plot.data$model, aes(x=x, y=y), color="blue") + 
  labs(x="k", y="CDF") + theme_bw() 














