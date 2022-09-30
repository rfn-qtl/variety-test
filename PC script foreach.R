#####################################
# Pre-Comercial test analysis - 2022
# Roberto Fritsche-Neto
# rfneto@agcenter.lsu.edu
# Last update: Sep 30 2022
#####################################

### loading files
#loading files
data <- read.csv("PC_2022_phenotype.csv")
metadata <- read.csv("2022 PC Entry List.csv")

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

entries <- intersect(unique(data$germplasmName), unique(metadata[,4]))
fillers <- outersect(unique(data$germplasmName), unique(metadata[,4]))

# coerce col to factors or numeric
library(magrittr)
library(dplyr)
cols <- colnames(data)[1:30]
data %<>%
  mutate_each_(funs(factor(.)), cols)
cols <- colnames(data)[length(colnames(data))-1]
data %<>%
  mutate_each_(funs(as.numeric(.)), cols)

data$row <- as.numeric(data$rowNumber)
data$col <- as.numeric(data$colNumber)

#################################
# fitting a model for each trial
#################################
require(foreach)
require(doParallel)
require(doMC)
library(SpATS)
library(lme4)
library(car)
library(data.table)
library(ggplot2)

# setting the number of cores that will be used
detectCores()
registerDoParallel(cores = detectCores()) # type the number of cores you want to use
getDoParWorkers()

trials <- unique(data$locationName)
st <- 1:length((trials))

system.time(
results.st <- foreach(i = st, 
                          .packages = c("SpATS", "car", "data.table"), 
                          .combine = "rbind", 
                          .multicombine = TRUE, 
                          .errorhandling = "remove",
                          .verbose = TRUE
 ) %dopar% {

      # subset the data  
    sample <- droplevels.data.frame(data[data$locationName == trials[i],])

      outlier <- boxplot.stats(sample$Yield.LSU_01.0000138)$out
      log.outliers <- data.frame(sample[sample$Yield.LSU_01.0000138 %in% outlier,])
      sample$Yield.LSU_01.0000138[sample$Yield.LSU_01.0000138 %in% outlier] <- NA
      
      nrow <- max(sample$row)
      ncol <- max(sample$col)
      nseg.row <- nrow
      nseg.col <- ncol
      
      fitF <- SpATS(response = "Yield.LSU_01.0000138", 
                                fixed = ~ 1, 
                                random = ~ replicate + rowNumber + colNumber, 
                                spatial = ~ PSANOVA(col, row, nseg = c(nseg.col, nseg.row)), 
                                genotype = "germplasmName", 
                                genotype.as.random = FALSE, 
                                data = sample)
      
      # obtaining the spatial trends - from raw data to BLUES and BLUPS
      #plot.SpATS(fitF) 
    
      # Estimate BLUEs
      blues <- predict.SpATS(fitF, which = "germplasmName")
      blues <- blues[,c(1,7,8)]
      colnames(blues)[2:3] <- c("BLUE", "sep_BLUE")
      
      # Now, run as random
      fitR <- SpATS(response = "Yield.LSU_01.0000138", 
                    fixed = ~ 1, 
                    random = ~ replicate + rowNumber + colNumber, 
                    spatial = ~ PSANOVA(col, row, nseg = c(nseg.col, nseg.row)), 
                    genotype = "germplasmName", 
                    genotype.as.random = TRUE, 
                    data = sample)

      # to obtain the heritability via the package function we can use 
      h2g <- getHeritability(fitR)
      # Broad-sense heritability based on Cullis method
      Vg <- fitR$var.comp["germplasmName"]
      ng <- length(unique(sample$germplasmName))
      C11_g <- fitR$vcov$C11_inv
      trC11_g <-sum(diag(C11_g))
      av2 <- 2/ng * (trC11_g - (sum(C11_g) - trC11_g) / ng-1) # mean var of a difference between genotypic BLUPS
      H2.Cullis <- 1 - av2 / (2 * Vg)
      # Estimate BLUPs for Grain Yield
      blups <- predict.SpATS(fitR, which = "germplasmName")
      blups <- blups[,c(7,8)]
      colnames(blups)[1:2] <- c("BLUP", "sep_BLUP")
      # Reliability
      rel <- mean(1 - blups$sep_BLUP^2 / fitR$var.comp["germplasmName"])
      # weights for ID's - adjust residual for further analysis
      vcov.mme <- fitR$vcov$C11_inv
      w <- diag(vcov.mme)
      
      output <- data.frame(blues,
                           w = w,
                           blups,
                           Location = as.character(unique(sample$locationName)), 
                           h2g = h2g,
                           r = rel,
                           H.cullis = H2.Cullis,
                           trait = fitR$model$response[1]
                          )
        }
)

output <- results.st
head(output)

# estimate how many trials were eliminated by singularities
length(unique(output$Location))
# eliminate trials with low Cullis heritability
unique(output$H.cullis)
output <- output[output$H.cullis > 0.50,]
# the final dimensions after QC
length(unique(output$Location))
# saving the file
saveRDS(output, "results.st")

##################################################
# second step - joint analysis
##################################################

# merge metadata and remove fillers
output <- merge(output, metadata)
length(unique(output$germplasmName))
all(unique(output$germplasmName) %in% entries)
all(entries %in% unique(output$germplasmName))
head(output)
dim(output)

# Fitting genotype by environment model - joint analysis
library(sommer)
fitMET <- mmer(BLUE ~ 1,
               random= ~ Location + germplasmName + germplasmName:Location,
               weights = w,
               rcov= ~ units,
               data=output, 
               verbose = FALSE)

# Broad-sense heritability
nloc <- length(unique(output$Location))
(h2g.MET <- vpredict(fitMET, h2 ~ V2 / ( V2 + V3/nloc))) # MET level

# predicting the BLUP - main effect
BLUPs.main <- predict.mmer(object = fitMET, classify = "germplasmName")
BLUPs.main$pvals
# reliability
mean(1 - BLUPs.main$pvals$standard.error^2 / summary(fitMET)$varcomp[2,1])

# predicting the BLUP per enviroment
BLUPs.MET <- predict.mmer(object = fitMET, classify = c("germplasmName","Location"))
BLUPs.MET$pvals

# this dataset will be used to request the year of cross or pedigree
saveRDS(BLUPs.MET, "BLUPs.MET") 
saveRDS(BLUPs.main, "BLUPs.main") 

rice.handbook.table <- round(reshape2::acast(BLUPs.MET$pvals, germplasmName ~  Location, value.var = "predicted.value"))
write.csv(rice.handbook.table, "rice.handbook.table.csv")

# barplot graph with confidence interval using main
data.plot <- BLUPs.main$pvals
data.plot <- merge(data.plot, metadata)
dim(data.plot)

limits <- aes(ymax = data.plot$predicted.value + data.plot$standard.error*1.96,
              ymin = data.plot$predicted.value - data.plot$standard.error*1.96)

p <- ggplot(data = data.plot, aes(x = reorder(germplasmName, -predicted.value), y = predicted.value, color = reorder(Gen_type, -predicted.value), 
                                  fill=reorder(BP, -predicted.value))) + 
  geom_bar(position = position_dodge(), 
         stat="identity") +
  scale_fill_brewer(palette="Set1") +
  geom_errorbar(limits, position = position_dodge(),
                width = 0.5) +
  labs(x = "Variety", y = "Yield") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(legend.position = "bottom", legend.box = "horizontal") +
  annotate("text", x=length(unique(data.plot$germplasmName))/2, 
           y=max(data.plot$predicted.value), 
           label= paste("Accuracy = ", round(sqrt(h2g.MET[1]), 2)*100, "%")) +
  labs(x = "", fill = "Breeding Program") +
  labs(color = "")

ggsave(filename = './Overall_performances.tiff',
       plot = p,
       device = 'tiff',
       width = 280,
       height = 140,
       units = 'mm',
       dpi = 300)


#######################################
# third step - MET analysis
#######################################
library(statgenGxE)

data.MET <- BLUPs.MET$pvals
head(data.MET)
colnames(data.MET)[c(2,4)] <- c("Variety", "Yield") 

## Create a TD object
dropsTD <- statgenSTA::createTD(data = data.MET, genotype = "Variety", trial = "Location")

## Fit a model where trials are nested within scenarios.
dropsVarComp <- gxeVarComp(TD = dropsTD, trait = "Yield")
summary(dropsVarComp)
vc(dropsVarComp)
## Compute heritability
herit(dropsVarComp)


# Finlay-Wilkinson Analysis
## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "Yield")
summary(dropsFW)

# reorganizing the data
test <- matrix(dropsFW$TD$`Lake Arthur`$beta) %*% matrix(sort(dropsFW$envEffs[,2]), nrow = 1) + dropsFW$estimates[,2]
locations <- c(trials) 
locations <- as.character(locations[match(sort(dropsFW$envEffs[,2]), dropsFW$envEffs[,2])])
colnames(test) <- sort(dropsFW$envEffs[,2])
rownames(test) <- dropsFW$estimates[, 1]
test <- data.frame(Variety = rownames(test), test)
colnames(test)[2:ncol(test)] <- sort(dropsFW$envEffs[,2])
test <- melt(test)
test$variable <- as.numeric(as.character(test$variable))


## Line plot for MET
q <- ggplot(data = BLUPs.MET$pvals, aes(x = reorder(Location, BLUPs.MET$pvals$predicted.value),
                                        y = predicted.value, 
                                        group = germplasmName, 
                                        colour = germplasmName)) + 
  geom_line() + geom_point() + 
  labs(x = "Location", y = "Yield") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename = './Actual_performance_across_locations.tiff',
       plot = q,
       device = 'tiff',
       width = 300,
       height = 400,
       units = 'mm',
       dpi = 300)


## Create line plot for Finlay Wilkinson analysis.
q <- ggplot(data = test, aes(x = variable, 
                               y = value, 
                               group = Variety, 
                               colour = Variety)) + 
  geom_line() + geom_point() + 
  labs(x = "Location", y = "Yield") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_vline(xintercept = meanEnV <- mean(dropsFW$envEffs[,2]), colour = "grey") +
  scale_x_continuous(breaks = round(dropsFW$envEffs[,2]), sec.axis = dup_axis(labels = locations))


ggsave(filename = './Stability_adpatability_across_locations.tiff',
       plot = q,
       device = 'tiff',
       width = 300,
       height = 400,
       units = 'mm',
       dpi = 300)


# GGE Biplot
library(metan)
head(data.MET)
model <- gge(data.MET, Location, Variety, Yield, svp = "symmetrical")
a <- plot(model, type = 1)
b <- plot(model, type = 2)
c <- plot(model, type = 3)
d <- arrange_ggplot(a, b, c, tag_levels = "a")

ggsave(filename = './GGE_Biplots.tiff',
       plot = d,
       device = 'tiff',
       width = 400,
       height = 250,
       units = 'mm',
       dpi = 300)

###### the end ###########