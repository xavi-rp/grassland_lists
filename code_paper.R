
###############################################################
#######            Code accompanying the paper          #######      
#######     "A machine learning approach to identify    #######
#######  effective parsimonious sets of indicator taxa  #######
#######       for grassland diversity monitoring"       #######      
###############################################################

# Code written by G. Strona


### libraries ####
library(randomForest)
library(mgcv)
library(pdp)
library(akima)
library(fields)

set.seed(123456)



### Clean data sets ####

# Species richness per site
data_div <- read.csv("data_div.csv", header = TRUE)
head(data_div)
nrow(data_div)

# relevé
releve <- read.csv("releve_clean.csv", header = TRUE)
View(releve)


### Evaluate 20 species list accuracy ####

r2_e <- round(cor.test(data_div$releve_diversity, data_div$experts_diversity_20)[[4]]**2, 2)
r2_e

pdf('fig2_old.pdf', width = 8, height = 8)
#png('fig2_old.png', width = 16, height = 8, units = "cm", res = 300, pointsize = 10)
par(mfrow = c(1, 1))
plot(data_div$experts_diversity_20, data_div$releve_diversity, 
     #main = bquote(paste('(a)' ~ {R^2} ~ '=' ~ .(r2_e))), 
     main = bquote(paste('' ~ {R^2} ~ '=' ~ .(r2_e))), 
     cex.main = 1.5, 
     las = 1,  cex.axis = 1, cex.lab = 1, pch = 16, 
     xlab = 'Indicator Richness\n(Expert Surveyors)', ylab = 'Total Richness')

dev.off()


### Genetic algorithm ####

### Explore which combination of species would provide the best correlation with total diversity
### Use a Genetic algorithm
### first expand the releve matrix with presence of families and genera,  so to be able
### to identify indicator groups of varying taxonomic resolution (as it was in the actual surveys)
View(releve)
releve[is.na(releve)]  <-  0
rel_mat <- as.matrix(sapply(releve[, 3:ncol(releve)],  as.numeric))
rel_mat <- 1*(rel_mat>0)
rel_div <- colSums(rel_mat)

### now expand the matrix
spp <- releve$species
genera <- sapply(strsplit(spp, '\\.'), "[[", 1)
genera_u <- sort(unique(genera))
fam_dict <- as.list(read.csv('fam_dict.csv', header=F))
fam_dict  <-  setNames(as.character(fam_dict$V2),  fam_dict$V1)
fams <- fam_dict[genera]
fams_u <- sort(unique(fams))

to_add_gen <- c()
for (gen in genera_u){
  sel_rows <- which(genera == gen)
  if (length(sel_rows) > 1){
    to_add_gen <- rbind(to_add_gen, 1 * (colSums(rel_mat[sel_rows, ]) > 0))}
  else {
    to_add_gen <- rbind(to_add_gen, rel_mat[sel_rows, ])
  }
}

to_add_fam <- c()
for (fam in fams_u){
  sel_rows <- which(fams == fam)
  if (length(sel_rows) > 1){
    to_add_fam <- rbind(to_add_fam, 1 * (colSums(rel_mat[sel_rows, ]) > 0))}
  else {
    to_add_fam <- rbind(to_add_fam, rel_mat[sel_rows, ])
  }
}

rel_mat <- rbind(rel_mat, to_add_gen, to_add_fam)
all_names <- c(spp, genera_u, fams_u)

n_ind_sp <- 20
pop_size <- 100
get_fit <- function(rc, rel_mat, rel_div){
  return(cor.test(colSums(rel_mat[unique(rc), ]), rel_div, na.rm = TRUE)[[4]]**2)
}


spp_n <- length(all_names)
mutate <- function(pop, recomb_rate = 0.001, mut_rate = 0.005, spp_n = 3320){
  N <- length(pop)
  C <- ncol(pop)
  rand_cells_n <- as.integer(N * recomb_rate)
  rand_cells_1 <- sample(1:N, rand_cells_n)
  rand_cells_2 <- sample(1:N, rand_cells_n)
  val1 <- pop[rand_cells_1]
  val2 <- pop[rand_cells_2]
  pop[rand_cells_2] <- val1
  pop[rand_cells_1] <- val2  
  pop[sample(1:N, N * mut_rate)] <- sample(1:spp_n, N * mut_rate)
  return (pop)
}


get_l <- function(x){return(length(unique(x)))}


### Run the following lines to obtain the best sets of indicator species 
### using the genetic algorithm; it takes a few hours. Alternatively, 
### run the following to load the data from simulations I run:

run_GA <-  "yes"
run_GA <-  "no"

if(run_GA == "yes"){
  best_pops <- c()
  for (ga_rep in 1:100){
    pop <- c()
    for (i in 1:pop_size){
      pop <- rbind(pop, sample(1:ncol(rel_mat), n_ind_sp))
    }
    for (rep in 1:5000){
      pop <- mutate(pop)
      fitness <- apply(pop,1,get_fit,rel_mat = rel_mat,rel_div = rel_div)
      pop <- rbind(pop[order(fitness, decreasing = TRUE)[1:10], ],
                   pop[sample(1:pop_size, pop_size, replace = TRUE, prob = fitness), ][1:(pop_size - 10), ])
    }
    best_pops <- rbind(best_pops, pop[which(fitness == max(fitness))[1], ])
    print(c(ga_rep, max(fitness)))
  }
  
  
  ### now run 10000 generations on the best_pops to see if there is room to improve
  pop_best <- best_pops
  for (rep in 1:10000){
    pop_best <- mutate(pop_best)
    fitness <- apply(pop_best, 1, get_fit, rel_mat = rel_mat, rel_div = rel_div)
    pop_best <- rbind(pop_best[order(fitness, decreasing = TRUE)[1:10], ],
                      pop_best[sample(1:pop_size, pop_size, replace = TRUE, prob = fitness), ][1:(pop_size - 10), ])
    print(c(rep, max(fitness)))
  }
  
  
  
  #### replicate for 12 categories
  n_ind_sp <- 12
  
  best_pops_12 <- c()
  for (ga_rep in 1:100){
    pop <- c()
    for (i in 1:pop_size){
      pop <- rbind(pop, sample(1:ncol(rel_mat), n_ind_sp))
    }
    for (rep in 1:5000){
      pop <- mutate(pop)
      fitness <- apply(pop, 1, get_fit, rel_mat = rel_mat, rel_div = rel_div)
      pop <- rbind(pop[order(fitness, decreasing = TRUE)[1:10], ],
                   pop[sample(1:pop_size, pop_size, replace = TRUE, prob = fitness), ][1:(pop_size - 10), ])
    }
    best_pops_12 <- rbind(best_pops_12, pop[which(fitness == max(fitness))[1], ])
    print(c(ga_rep, max(fitness)))
  }
  
  
  
  pop_best_12 <- best_pops_12
  for (rep in 1:10000){
    pop_best_12 <- mutate(pop_best_12)
    fitness <- apply(pop_best_12, 1, get_fit, rel_mat = rel_mat, rel_div = rel_div)
    pop_best_12 <- rbind(pop_best_12[order(fitness, decreasing = TRUE)[1:10], ],
                         pop_best_12[sample(1:pop_size,pop_size,replace=T,prob=fitness),][1:(pop_size-10),])
    print(c(rep, max(fitness)))
  }
  
}else{
  load('best_pops.Rdata', verbose = TRUE)
}

### Could be confounding,  however:
# pop_best != best_pops; the first is the ideal best optimization
# obtained by using as starting population the best 100 populations
# in the 100 independent experiments (i.e. best_pops); hence
# best_pops has more variation then pop_best; if we want to focus just on
# theoretical optimum, we should refer to pop_best; if we want to focus
# on different potential sets we should refer to best_pops; same reasoning
# for pop_best_12 vs best_pops_12



### examine the identity of best models
#1) tot number of categories

get_tax_res <- function(x, sp, gen, fam, names){
  n <- length(x)
  x <- names[x]
  p_sp <- length(intersect(x, sp)) / n
  p_gen <- length(intersect(x, gen)) / n
  p_fam <- length(intersect(x, fam)) / n
  return(c(p_sp, p_gen, p_fam))
}



### do plots for the 20 species
best_fitness <- apply(pop_best, 1, get_fit, rel_mat = rel_mat, rel_div = rel_div)
summary(best_fitness)
sd(best_fitness)


best_set <- pop_best[which(best_fitness == max(best_fitness))[1], ]

pdf('fig1.pdf', width = 9, height = 3)
#png('fig1.png', width = 18, height = 6, units = "cm", res = 300, pointsize = 10)
par(mfrow = c(1, 3), font.main = 1)

r2 <- round(cor.test(rel_div, colSums(rel_mat[unique(best_set), ]), na.rm = TRUE)[[4]]**2, 2)
plot(colSums(rel_mat[unique(best_set), ]), rel_div,  
     #main = paste('(a) R2 =', r2), 
     main = bquote(paste('(a)' ~ {R^2} ~ '=' ~ .(r2))), 
     cex.main = 1.5,
     las = 1,  cex.axis = 1, cex.lab = 1, pch = 16, 
     xlab = 'Indicator Richness (Best Set)', ylab = 'Total Richness'
)

taxa_frac <- t(apply(best_pops, 1, get_tax_res, sp = spp, gen = genera_u, fam = fams_u, names = all_names))
boxplot(taxa_frac, names = c('Species', 'Genera', 'Families'),
        las = 1,  cex.axis = 1, cex.lab = 1, 
        xlab = 'Taxonomic Level', ylab = 'Prevalence In Best Sets', main = '(b)', cex.main = 1.5)

## former Fig2 (see above)
plot(data_div$experts_diversity_20, data_div$releve_diversity, 
     main = bquote(paste('(c)' ~ {R^2} ~ '=' ~ .(r2_e))), 
     cex.main = 1.5, 
     las = 1,  cex.axis = 1, cex.lab = 1, pch = 16, 
     xlab = 'Indicator Richness\n(Expert Surveyors)', ylab = 'Total Richness')


dev.off()


### do plots for the 12 groups
best_fitness <-  apply(pop_best_12, 1, get_fit, rel_mat = rel_mat, rel_div = rel_div)
best_set_12 <- pop_best_12[which(best_fitness==max(best_fitness))[1], ]

summary(best_fitness)
sd(best_fitness)


pdf('fig2.pdf', width = 6, height = 3)
#png('fig2.png', width = 12, height = 6, units = "cm", res = 300, pointsize = 6)
par(mfrow = c(1, 2), font.main = 1)

r2 <- round(cor.test(rel_div, colSums(rel_mat[unique(best_set_12), ]), na.rm = TRUE)[[4]]**2, 2)
plot(colSums(rel_mat[unique(best_set_12), ]), rel_div,  
     #main=paste('(a) R2 =', r2), 
     main = bquote(paste('(a)' ~ {R^2} ~ '=' ~ .(r2))), 
     cex.main = 1.5,
     las = 1,  cex.axis = 1, cex.lab = 1, pch = 16, 
     xlab = 'Indicator Richness\n (Best Set, 12 Indicators)', ylab = 'Total Richness'
)

best_taxa_12 <- all_names[unique(as.vector(best_pops_12))]

taxa_frac <- t(apply(best_pops_12, 1, get_tax_res, sp=spp, gen=genera_u, fam=fams_u, names=all_names))
boxplot(taxa_frac, names = c('Species', 'Genera', 'Families'),
        las = 1,  cex.axis = 1, cex.lab = 1, 
        xlab = 'Taxonomic Level', ylab = 'Prevalence In Best Sets', main = '(b)', cex.main = 1.5)

dev.off()


### Prevalence ####
### look at species prevalence in the selected and not selected groups
# 'rel_mat' is the relevé matrix, including genera and families
prev_ind_taxa <- rowSums(rel_mat) / ncol(rel_mat)

### for the 20 species list
prev_data <- c()
for (rep in 1:100){
  n_sp <- length(intersect(all_names[best_pops[rep, ]], spp))
  n_gen <- length(intersect(all_names[best_pops[rep, ]], genera_u))
  n_fam <- length(intersect(all_names[best_pops[rep, ]], fams_u))
  for (sub_rep in range(100)){
    rand_set <- match(c(sample(fams_u, n_fam), 
                        sample(genera_u, n_gen), 
                        sample(spp, n_sp)), all_names)
    prev_data <- rbind(prev_data, cbind(prev_ind_taxa[rand_set], 
                                        prev_ind_taxa[best_pops[rep, ]]))
  }
}

### do the same for the 12 species list
prev_data_12 <- c()
for (rep in 1:100){
  n_sp <- length(intersect(all_names[best_pops_12[rep, ]], spp))
  n_gen <- length(intersect(all_names[best_pops_12[rep, ]], genera_u))
  n_fam <- length(intersect(all_names[best_pops_12[rep, ]], fams_u))
  for (sub_rep in range(100)){
    rand_set <- match(c(sample(fams_u, n_fam), 
                        sample(genera_u, n_gen), 
                        sample(spp, n_sp)), all_names)
    prev_data_12 <- rbind(prev_data_12, cbind(prev_ind_taxa[rand_set], 
                                              prev_ind_taxa[best_pops_12[rep, ]]))
  }
}

# correspondence between key-taxa and relevé species
ind_sp_dict <- read.csv('matched_names_ok1.csv', header = TRUE) 

# figure 3a
boxplot(prev_ind_taxa[ind_sp_dict[, 2]], 
        prev_data[, 1], prev_data[, 2], 
        prev_data_12[, 1], prev_data_12[, 2], 
        names = c('LUCAS Indicators 20',
                  'Random Sets 20', 'Best Sets GA 20', 
                  'Random Sets 12', 'Best Sets GA 12'),
        outline = F,  las = 1, 
        cex.axis = 0.8, cex.lab = 1.2, 
        xlab = 'Prevalence Across Sites', cex = 0.7, horizontal = TRUE,
        ylim = c(0,1))



### Detectability ####
plant_net <- read.csv('plantnet_v2.csv', header = TRUE)

plant_net_sp <- aggregate(cbind(count, score)~searched.name, FUN = 'mean', data = plant_net)
plant_net_gen <- aggregate(cbind(count, score)~genus, FUN = 'mean', data = plant_net)
plant_net_fam <- aggregate(cbind(count, score)~family, FUN = 'mean', data = plant_net)

colnames(plant_net_sp)[1] <- 'taxon'
colnames(plant_net_gen)[1] <- 'taxon'
colnames(plant_net_fam)[1] <- 'taxon'

plant_net_all <- rbind(plant_net_sp, plant_net_gen, plant_net_fam)


### create a dataset with species in the same order as all_names for plant_net data; non matched taxa get scores 0 0
plant_net_scores <- data.frame(cbind(all_names, 0, 0))
colnames(plant_net_scores) <- colnames(plant_net_all)
matched_plant_net_data <- match(plant_net_all[, 1], all_names)
plant_net_all <- plant_net_all[!is.na(matched_plant_net_data), ]
matched_plant_net_data <- matched_plant_net_data[!is.na(matched_plant_net_data)]

plant_net_scores[matched_plant_net_data, 2:3] <- plant_net_all[, 2:3]


### do plots for the 20 species

pn_data <- c()
for (rep in 1:100){
  n_sp <- length(intersect(all_names[best_pops[rep, ]], spp))
  n_gen <- length(intersect(all_names[best_pops[rep, ]], genera_u))
  n_fam <- length(intersect(all_names[best_pops[rep, ]], fams_u))
  for (sub_rep in range(100)){                                      
    rand_set <- match(c(sample(fams_u, n_fam), 
                        sample(genera_u, n_gen), 
                        sample(spp, n_sp)), all_names)
    pn_data <- rbind(pn_data, cbind(plant_net_scores[rand_set, 2:3], 
                                    plant_net_scores[best_pops[rep, ], 2:3]))
  }
}


### do the same for the 12 species list
pn_data_12 <- c()
for (rep in 1:100){
  n_sp <- length(intersect(all_names[best_pops_12[rep, ]], spp))
  n_gen <- length(intersect(all_names[best_pops_12[rep, ]], genera_u))
  n_fam <- length(intersect(all_names[best_pops_12[rep, ]], fams_u))
  for (sub_rep in range(100)){
    rand_set <- match(c(sample(fams_u, n_fam), 
                        sample(genera_u, n_gen), 
                        sample(spp, n_sp)), all_names)
    pn_data_12 <- rbind(pn_data_12, cbind(plant_net_scores[rand_set, 2:3], 
                                          plant_net_scores[best_pops_12[rep, ], 2:3]))
  }
}





pdf('fig3.pdf', width = 9, height = 3)
#png('fig3.png', width = 17, height = 6, units = "cm", res = 300, pointsize = 10)
par(mfrow = c(1, 3), font.main = 1) 
par(oma = c(0, 3, 0, 0))
#par(mai = c(1, 1.0, 0.5, 0.2))
par(mai = c(0.7, 0.7, 0.5, 0.2))

# fig 3a (see above)
boxplot(prev_ind_taxa[ind_sp_dict[, 2]], 
        prev_data[, 1], prev_data[, 2], 
        prev_data_12[, 1], prev_data_12[, 2], 
        names = c('LUCAS Indicators 20',
                  'Random Sets 20', 'Best Sets GA 20', 
                  'Random Sets 12', 'Best Sets GA 12'),
        outline = F,  las = 1, 
        cex.axis = 0.8, cex.lab = 1.2, 
        xlab = 'Prevalence Across Sites', cex = 0.7, horizontal = TRUE,
        main = '(a)',
        ylim = c(0,1))

boxplot(as.numeric(plant_net_scores[ind_sp_dict[, 2], 2]) / 1000, 
        #as.numeric(plant_net_scores[ind_sp_dict[which(ind_sp_dict[, 1]%in%subset_12), 2], 2])/1000, 
        as.numeric(pn_data[, 1]) / 1000, as.numeric(pn_data[, 3]) / 1000, 
        as.numeric(pn_data_12[, 1]) / 1000, as.numeric(pn_data_12[, 3]) / 1000, 
        names = c('LUCAS Indicators 20', 
                  #'LUCAS Indicators 12', 
                  'Random Sets 20', 'Best Sets GA 20', 
                  'Random Sets 12', 'Best Sets GA 12'),
        outline = F, las = 1, 
        cex.axis = 0.8, cex.lab = 1.2,
        xlab = 'PlantNet Count/1000', 
        cex = 0.7, horizontal = TRUE, main = '(b)')

boxplot(as.numeric(plant_net_scores[ind_sp_dict[, 2], 3]), 
        #as.numeric(plant_net_scores[ind_sp_dict[which(ind_sp_dict[, 1]%in%subset_12), 2], 3]), 
        as.numeric(pn_data[, 2]), as.numeric(pn_data[, 4]), 
        as.numeric(pn_data_12[, 2]), as.numeric(pn_data_12[, 4]), 
        names = c('LUCAS Indicators 20', 
                  #'LUCAS Indicators 12', 
                  'Random Sets 20', 'Best Sets GA 20', 
                  'Random Sets 12', 'Best Sets GA 12'), 
        outline = F, las = 1, 
        cex.axis = 0.8, cex.lab = 1.2,
        xlab = 'PlantNet Score', cex = 0.7, horizontal = TRUE, main = '(c)')

dev.off()



### in case you run the genetic algorithm and want to save the suggested sets
#save(list=c('best_pops', 'best_pops_12', 'pop_best', 'pop_best_12'), file='best_pops.Rdata')


### assessing significance

## PlantNet Count
summary(as.numeric(plant_net_scores[ind_sp_dict[, 2], 2]) / 1000) # LUCAS Indicators 20
summary(as.numeric(pn_data[, 1]) / 1000) # Random Sets 20
summary(as.numeric(pn_data[, 3]) / 1000)  # Best Sets GA 20


ks.test(as.numeric(plant_net_scores[ind_sp_dict[, 2], 2]) / 1000, # LUCAS Indicators 20
        as.numeric(pn_data[, 1]) / 1000)  # Random Sets 20

ks.test(as.numeric(plant_net_scores[ind_sp_dict[, 2], 2]) / 1000, # LUCAS Indicators 20
        as.numeric(pn_data[, 3]) / 1000)  # Best Sets GA 20

ks.test(as.numeric(pn_data[, 3]) / 1000,  # Best Sets GA 20
        as.numeric(pn_data[, 1]) / 1000) # Random Sets 20


ks.test(as.numeric(pn_data_12[, 1]) / 1000, # Random Sets 12
        as.numeric(pn_data_12[, 3]) / 1000) # Best Sets GA 12



## PlantNet Score
summary(as.numeric(plant_net_scores[ind_sp_dict[, 2], 3])) # LUCAS Indicators 20
summary(as.numeric(pn_data[, 2])) # Random Sets 20
summary(as.numeric(pn_data[, 4]))  # Best Sets GA 20


ks.test(as.numeric(plant_net_scores[ind_sp_dict[, 2], 3]), # LUCAS Indicators 20
        as.numeric(pn_data[, 2]))  # Random Sets 20

ks.test(as.numeric(plant_net_scores[ind_sp_dict[, 2], 3]), # LUCAS Indicators 20
        as.numeric(pn_data[, 4]))  # Best Sets GA 20

ks.test(as.numeric(pn_data[, 4]),  # Best Sets GA 20
        as.numeric(pn_data[, 2])) # Random Sets 20


ks.test(as.numeric(pn_data_12[, 2]), # Random Sets 12
        as.numeric(pn_data_12[, 4])) # Best Sets GA 12





