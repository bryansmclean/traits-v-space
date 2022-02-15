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
data.dirs <- data.dirs[grep('GRSM|ORNL|MLBS|BLAN|SCBI|BART|HARV', data.dirs)] # just the Appalachians sites

# loop through files to build full record-set
for(i in data.dirs){
  
  # read file-wise retaining smammal captures
  setwd(i)
  data.file <- list.files(pattern = 'pertrapnight')
  tmp <- fread(data.file)
  caps <- tmp[grep("4 -|5 -", tmp$trapStatus), ] # all captures
  indivs <- caps[which(caps$tagID != ""), ] # all captures with a unique identifier
  
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

write.csv(cap.counts, '/Users/mclean/Box/Projects/GI_smammal_community-traits/_data_/NEON_Appalachians-smammals_total-captures.csv', row.names = T, quote = F)
write.csv(indiv.counts, '/Users/mclean/Box/Projects/GI_smammal_community-traits/_data_/NEON_Appalachians-smammals_total-individuals.csv', row.names = T, quote = F)








































