######################################################
# harvesting smammal community data 
# for NEON box trapping
# Bryan S McLean 
# 10 Feb 2022
######################################################

library(data.table)


##############################################
## read downloaded NEON box-trapping files
##############################################

# index records dirs
setwd('/Users/mclean/Box/Projects/FuTRES/traits-v-space/_data_/NEON_count-small-mammals-7Feb2022/')
top.dir <- getwd()
data.dirs <- list.dirs()
data.dirs <- data.dirs[grepl('expanded', data.dirs, ignore.case = T, fixed = F)]
sites <- 'GRSM|ORNL|MLBS|BLAN|SCBI|BART|HARV' # list of just the Appalachians sites
data.dirs <- data.dirs[grep(sites, data.dirs)]
Total_traps <- rep(0, length(strsplit(sites, '|', fixed = T)[[1]])); names(Total_traps) <- strsplit(sites, '|', fixed = T)[[1]] # data frame of the total trap effort, per site

# loop through files to build full record-set
for(i in data.dirs){
  
  # read file-wise retaining desired info
  setwd(i)
  data.file <- list.files(pattern = 'pertrapnight')
  tmp <- fread(data.file)
  caps <- tmp[grep("4 -|5 -", tmp$trapStatus), ] # data frame of all captures
  indivs <- caps[which(caps$tagID != ""), ] # data frame of all captures with a unique identifier
  
  # record the total trapping effort
  Total_traps[unique(tmp$siteID)] <- Total_traps[unique(tmp$siteID)] + nrow(tmp)

  # append all captures to file
  if(i == data.dirs[1]){
    caps.all <- caps
  } else {
    caps.all <- rbind(caps.all, caps)
  }
  
  # append only known individuals to file
  if(i == data.dirs[1]){
    indivs.all <- indivs
  } else {
    indivs.all <- rbind(indivs.all, indivs)
  }
  
  # get an update
  print(paste('completed search for record set ', i, sep = '')); Sys.sleep(0.00001)
  
  # loop back up
  setwd(top.dir)
  
  if(i == data.dirs[1]){
    plot(nrow(caps),nrow(indivs), xlim = c(0,2500), ylim = c(0,2500))
  } else {
    points(nrow(caps),nrow(indivs))
  }
}


##############################################
## create (NEON)site-by-species matrix
##############################################

# get all unique sites and spp in the data set
sites <- unique(caps.all$siteID)
sites <- sites[order(sites)]
spp <- unique(caps.all$scientificName)[grep(" ", unique(caps.all$scientificName))] # removing species with no ID
spp <- spp[order(spp)]

# want total captures and total abundances
cap.counts <- indiv.counts <- matrix(0, nrow = length(spp), ncol = length(sites), dimnames = list(spp, sites))
Total_captures <- rep(NA, length(sites))
cap.counts <- rbind(cap.counts, Total_captures)
indiv.counts <- rbind(indiv.counts, Total_captures)

for(i in sites){
  
  tmp <- caps.all[caps.all$siteID == i, ]
  tmp <- tmp[grep(" ", tmp$scientificName), ]
  cap.counts[names(table(tmp$scientificName)), i] <- table(tmp$scientificName)
  cap.counts['Total_captures', i] <- nrow(tmp)
  rm(tmp)
  
  tmp <- indivs.all[indivs.all$siteID == i,]
  tmp <- tmp[grep(" ", tmp$scientificName), ]
  indiv.counts[names(table(tmp$scientificName)), i] <- table(tmp$scientificName)
  indiv.counts['Total_captures', i] <- nrow(tmp)
  rm(tmp)

  }

cap.counts <- rbind(cap.counts, Total_traps = Total_traps[colnames(cap.counts)])
indiv.counts <- rbind(indiv.counts, Total_traps = Total_traps[colnames(indiv.counts)])

write.csv(cap.counts, '/Users/mclean/Box/Projects/GI_smammal-community-traits/_data_/NEON_Appalachians-smammals_total-captures.csv', row.names = T, quote = F)
write.csv(indiv.counts, '/Users/mclean/Box/Projects/GI_smammal-community-traits/_data_/NEON_Appalachians-smammals_total-individuals.csv', row.names = T, quote = F)








































