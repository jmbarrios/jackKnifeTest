library(rgdal)
library(raster)
library(ENMeval)
library(dplyr)
library(magrittr)
library(iterpc)

covarFileList <- tools::list_files_with_exts('data/', "tif")
env <- raster::stack(covarFileList)

occs <- read.csv('data/data_wo_duplicates.csv')

occs <- dplyr::select(occs, Dec_Long, Dec_Lat)
bg.df <- as.data.frame(dismo::randomPoints(env[[1]], n = 10000))

sp.models <- ENMeval::ENMevaluate(occs, env, bg.df, RMvalues = c(1),
                                  fc = c("L"),
                                  method = "jackknife", bin.output = TRUE,
                                  updateProgress = TRUE)

resultsData <- sp.models@results
resultsData %<>% select(settings, contains('ORmin_'), contains('ProbMin'))
probs <- resultsData %>% select(contains('ProbMin')) %>% unlist(., use.names=TRUE)

zProb <- function(state, probs) {
  success <- probs^state
  fails <- (1-probs)^(1-state)
  
  return(prod(success, fails))
}

zStat <- function(state, probs) {
  d <- sum(state*(1-probs))
  pD <- zProb(state, probs)
  
  r <- c(d,pD)
  names(r) <- c('D', 'pD')
  
  return(r)
}

I <- iterpc(table(c(0, 1)), 3, ordered = TRUE, replace = TRUE)
no_perm <- getlength(I)
results <- data.frame(D = double(), pD = double())

for (i in 1:no_perm) {
  s <- getnext(I)
  z_statistic <- zStat(s, u)
  
  results <- rbind(results, z_statistic)
}
names(results) <- c('D', 'pD')

results %>% filter(D >= 1.5) %>% summarise_at('pD', sum)
