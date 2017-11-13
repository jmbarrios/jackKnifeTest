# Basado en los datos del material suplementario de 
# Pearson, R. G., Raxworthy, C. J., Nakamura, M. and Townsend Peterson, A. 
# (2007), ORIGINAL ARTICLE: Predicting species distributions from small 
# numbers of occurrence records: a test case using cryptic geckos in Madagascar.
# Journal of Biogeography, 34: 102–117. doi:10.1111/j.1365-2699.2006.01594.x

library(iterpc)

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

dataSpecies <- readr::read_csv('data/Example_data.txt')

sp1 <- dataSpecies %>% filter(specie == 'sp1') %>% select(state, prob)
sp2 <- dataSpecies %>% filter(specie == 'sp2') %>% select(state, prob)
sp3 <- dataSpecies %>% filter(specie == 'sp3') %>% select(state, prob)

## Analisis specie 1
I <- iterpc::iterpc(table(c(0, 1)), nrow(sp1), ordered = TRUE, replace = TRUE)
no_perm <- getlength(I)
results <- data.frame(D = double(), pD = double())
probs <- sp1 %>% select(prob) %>% unlist(., use.names = TRUE)

for (i in 1:no_perm) {
  s <- getnext(I)
  z_statistic <- zStat(s, probs)
  
  results <- rbind(results, z_statistic)
}
names(results) <- c('D', 'pD')

obsState <- sp1 %>% select(state) %>% unlist(., use.names = TRUE)
obsStat <- zStat(obsState, probs)

p_value <- results %>% filter(D >= obsStat[1]) %>% summarise_at('pD', sum)

## Analisis specie 2
I <- iterpc::iterpc(table(c(0, 1)), nrow(sp2), ordered = TRUE, replace = TRUE)
no_perm <- getlength(I)
results <- data.frame(D = double(), pD = double())
probs <- sp2 %>% select(prob) %>% unlist(., use.names = TRUE)

for (i in 1:no_perm) {
  s <- getnext(I)
  z_statistic <- zStat(s, probs)
  
  results <- rbind(results, z_statistic)
}
names(results) <- c('D', 'pD')

obsState <- sp2 %>% select(state) %>% unlist(., use.names = TRUE)
obsStat <- zStat(obsState, probs)

p_value <- results %>% filter(D >= obsStat[1]) %>% summarise_at('pD', sum)

## Analisis specie 3
I <- iterpc::iterpc(table(c(0, 1)), nrow(sp3), ordered = TRUE, replace = TRUE)
no_perm <- getlength(I)
results <- data.frame(D = double(), pD = double())
probs <- sp3 %>% select(prob) %>% unlist(., use.names = TRUE)

for (i in 1:no_perm) {
  s <- getnext(I)
  z_statistic <- zStat(s, probs)
  
  results <- rbind(results, z_statistic)
}
names(results) <- c('D', 'pD')

obsState <- sp3 %>% select(state) %>% unlist(., use.names = TRUE)
obsStat <- zStat(obsState, probs)

p_value <- results %>% filter(D >= obsStat[1]) %>% summarise_at('pD', sum)