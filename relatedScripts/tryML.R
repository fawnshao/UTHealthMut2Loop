# devtools::install_github('thomasp85/ggraph')
# install.packages(c("dplyr", "ggraph", "igraph"))
library(randomForest)
library(dplyr)
library(ggraph)
library(igraph)
setwd("/Volumes/GoogleDrive/My Drive/housekeeping_genes/ML_test")
data <- read.csv("GTEx.GTRD.forWeka.csv", row.names = 1)
tree_func <- function(final_model, 
                      tree_num) {
    
    # get tree by index
    tree <- randomForest::getTree(final_model, 
                                  k = tree_num, 
                                  labelVar = TRUE) %>%
        tibble::rownames_to_column() %>%
        # make leaf split points to NA, so the 0s won't get plotted
        mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
    
    # prepare data frame for graph
    graph_frame <- data.frame(from = rep(tree$rowname, 2),
                              to = c(tree$`left daughter`, tree$`right daughter`))
    
    # convert to graph and delete the last node that we don't want to plot
    graph <- graph_from_data_frame(graph_frame) %>%
        delete_vertices("0")
    
    # set node labels
    V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
    V(graph)$leaf_label <- as.character(tree$prediction)
    V(graph)$split <- as.character(round(tree$`split point`, digits = 2))
    
    # plot
    plot <- ggraph(graph, 'dendrogram') + 
        theme_bw() +
        geom_edge_link() +
        geom_node_point() +
        geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
        geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
        geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                        repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.background = element_blank(),
              plot.background = element_rect(fill = "white"),
              panel.border = element_blank(),
              axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(size = 18))
    
    print(plot)
}


rf <- randomForest(x = data[,-567], y = data[, 567], 
                          ntree = 5000, mtry = 100, importance = T, 
                   keep.inbag = T, keep.forest = T)
rf.unsupervised <- randomForest(x = data[,-567], 
                   ntree = 5000, mtry = 100, importance = T, 
                   keep.inbag = T, keep.forest = T)

rf <- randomForest(type ~ . , data = data,  ntree = 400)
ff < forestFloor(rf, X, binary_reg = T, calc_np = T)
Col < fcol(ff,cols=1,outlier.lim = 2.5)
plot(ff, col=Col, plot_GOF = T)

tree_func(rf, 5)
tree_num <- which(rf$finalModel$forest$ndbigtree == min(rf$finalModel$forest$ndbigtree))
tree_func(final_model = rf$finalModel, tree_num)











library(MASS)
attach(Boston)
set.seed(101)
train <- sample(1:nrow(Boston), 300)
Boston.rf <- randomForest(medv ~ . , data = Boston , subset = train)
Boston.rf
oob.err <- double(13)
test.err <- double(13)

# mtry is no of Variables randomly chosen at each split
for(mtry in 1:13) 
{
    rf <- randomForest(medv ~ . , data = Boston, subset = train, mtry = mtry, ntree = 400) 
    oob.err[mtry] <- rf$mse[400] #Error of all Trees fitted
    
    pred <- predict(rf, Boston[-train,]) #Predictions on Test Set for each Tree
    test.err[mtry] <- with(Boston[-train,], mean( (medv - pred)^2)) #Mean Squared Test Error
    
    cat(mtry," ") #printing the output to the console
}


data(iris)
rf1 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
rf2 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
rf3 <- randomForest(Species ~ ., iris, ntree=50, norm.votes=FALSE)
rf.all <- combine(rf1, rf2, rf3)
print(rf.all)

set.seed(1)
data(iris)
iris.rf <- randomForest(Species ~ ., iris, keep.forest=FALSE)
plot(margin(iris.rf))




# Generate scaled 4*5 matrix with random std normal samples
set.seed(101)
mat <- scale(matrix(rnorm(20), 4, 5))
dimnames(mat) <- list(paste("Sample", 1:4), paste("Var", 1:5))

# Perform PCA
myPCA <- prcomp(mat, scale. = F, center = F)
myPCA$rotation # loadings
myPCA$x # scores
biplot(myPCA)
plot(myPCA$x[,1:2], col = )
# Perform SVD
mySVD <- svd(mat)
mySVD # the diagonal of Sigma mySVD$d is given as a vector
sigma <- matrix(0,4,4) # we have 4 PCs, no need for a 5th column
diag(sigma) <- mySVD$d # sigma is now our true sigma matrix



install.packages(c("mlbench","forestFloor", "AUC"))
rm(list=ls())
set.seed(1)
library(mlbench)
library(randomForest)
library(forestFloor)
library(AUC)
data(PimaIndiansDiabetes)
y = PimaIndiansDiabetes$diabetes
X = PimaIndiansDiabetes
X = X[,!names(X)=="diabetes"]

#train default model and the most regularized model with same predictive performance
rf.default = randomForest(X,y,ntree=5000)
rf.robust = randomForest(X,y,sampsize=25, ntree=5000, mtry=4,
                         keep.inbag = T,keep.forest = T)
#verify similar performance
plot(roc(rf.default$votes[,2],y), main="ROC: default black, robust is red")
plot(roc(rf.robust$votes[,2],y), col = 2, add = T)
auc(roc(rf.default$votes[,2],y))
auc(roc(rf.robust$votes[,2],y))

#compute feature contributions
ff = forestFloor(rf.robust,X,binary_reg = T,calc_np=T)
Col = fcol(ff,cols=1,outlier.lim = 2.5)

#the plot in this blog
plot(ff,col=Col,plot_GOF = T)
#some 3d plots
show3d(ff,c(1,5),5,col=Col,plot_GOF = T)
library(rgl); rgl.snapshot("3dPressure.png")


tree_num <- which(rf.default$finalModel$forest$ndbigtree == min(rf.default$finalModel$forest$ndbigtree))

tree_func(final_model = rf.default$finalModel, 5)
