######################################################
# harvesting smammal community data 
# for NEON box trapping
# Bryan S McLean 
# 10 Feb 2022, last updated 28 Feb 2023
######################################################

library(data.table)


##############################################
## read NEON box-trapping files
## retain captures and record effort
##############################################

# index the occurrence record dirs
setwd("~/Library/CloudStorage/Box-Box/Projects/Chapman_GI_smammal-community-traits/_data_/NEON_count-small-mammals_2023-02-24")
top.dir <- getwd()
data.dirs <- list.dirs()
data.dirs <- data.dirs[grepl('expanded', data.dirs, ignore.case = T, fixed = F)]
sites <- 'GRSM|ORNL|MLBS|BLAN|SCBI|BART|HARV' # list of just Appalachians sites (this should be all that are present - was a targeted download)
data.dirs <- data.dirs[grep(sites, data.dirs)]

# data frame of the total trap effort, per site
Total_effort <- rep(0, length(strsplit(sites, '|', fixed = T)[[1]]))
# data frames of the timing of trap effort, per site
Total_effort_earliestStartDate <- rep(365, length(strsplit(sites, '|', fixed = T)[[1]]))
Total_effort_latestEndDate <- rep(1, length(strsplit(sites, '|', fixed = T)[[1]]))
names(Total_effort) <- names(Total_effort_earliestStartDate) <- names(Total_effort_latestEndDate) <- strsplit(sites, '|', fixed = T)[[1]] 

# loop through files to build full record set
for(i in data.dirs){
  
	# read file-wise
 	setwd(i)
 	data.file <- list.files(pattern = 'pertrapnight')
 	tmp <- fread(data.file)
 	if(length(grep("1 -", tmp$trapStatus)) == nrow(tmp)) {setwd(top.dir); next} # skip files that have no trap effort
	if(length(grep("1 -", tmp$trapStatus)) > 0) {tmp <- tmp[-grep("1 -", tmp$trapStatus), ]}  # limit to actual trapnights (i.e., excluding traps not set)
	
	# add Julian days to the record set
	julianday <- rep(NA, nrow(tmp))
	for(j in 1:nrow(tmp)){julianday[j] <- strptime(tmp$collectDate[j],"%Y-%m-%d")$yday} #convert to julian day
	tmp <- cbind(tmp, julianday)

	# record trap effort, timing, captures, individuals
	Total_effort[unique(tmp$siteID)] <- Total_effort[unique(tmp$siteID)] + nrow(tmp) # add the trap effort to site total
	eff.start <- min(tmp$julianday, na.rm = F) # the earliest date of trapping at this site
	eff.end <- max(tmp$julianday, na.rm = F) # the latest date of trapping at this site
	caps <- tmp[grep("4 -|5 -", tmp$trapStatus), ] # data frame of captures (includes multiple captures in single trap)
	inds <- tmp[which(tmp$tagID != ""), ] # data frame of captures that have a unique identifier
  
  	# track the magnitude and timing of effort (to control for effort across sites)
	if(i == data.dirs[1]){
		Total_effort_earliestStartDate[unique(tmp$siteID)] <- eff.start
		Total_effort_latestEndDate[unique(tmp$siteID)] <- eff.end
  		} else {
		if(eff.start < Total_effort_earliestStartDate[unique(tmp$siteID)]) {Total_effort_earliestStartDate[unique(tmp$siteID)] <- eff.start} # update with earlier dates
		if(eff.end > Total_effort_latestEndDate[unique(tmp$siteID)]) {Total_effort_latestEndDate[unique(tmp$siteID)] <- eff.end} # update with later dates
		}		

	# append the captures to file
	if(i == data.dirs[1]){
		caps.all <- caps
  		} else {
  		caps.all <- rbind(caps.all, caps)
		}
  
	# append just the known individuals to file
	if(i == data.dirs[1]){
    	inds.all <- inds
  		} else {
    	inds.all <- rbind(inds.all, inds)
  		}
  
	# update on progress
	print(paste('read record set ', i, sep = '')); Sys.sleep(0.0000001)
  
	# loop up to next record set/dir
	setwd(top.dir)

	# summary plot
	if(i == data.dirs[1]){
    	plot(nrow(inds), nrow(caps), 
    		xlim = c(0,1000), 
    		ylim = c(0,1000), 
    		xlab = "Total Inds.", 
    		ylab = "Total Captures", 
    		col = 'darkblue', 
    		main = "Unique Individuals and Captures per Month"
    		)
    	abline(a = 1, b = 1, lty = 2, col = 'darkblue')
  		} else {
    	points(nrow(inds), nrow(caps), col = 'darkblue')
    	}
    
    rm(tmp);rm(caps);rm(inds);rm(eff.start);rm(eff.end);
    
    }


##############################################
## create community (site x species) matrices
##############################################

# get unique sites and spp from the capture data
sites <- unique(caps.all$siteID)
sites <- sites[order(sites)]
spp <- unique(caps.all$scientificName)[grep(" ", unique(caps.all$scientificName))] # keep only species with some taxonomy data
spp <- spp[order(spp)]

# create matrices to hold total captures and total individuals
cap.counts <-  ind.counts <- matrix(0, nrow = length(spp), ncol = length(sites), dimnames = list(spp, sites))
Totals <- rep(NA, length(sites)) # row to hold totals
cap.counts <- rbind(cap.counts, Totals, Total_effort, Total_effort_earliestStartDate, Total_effort_latestEndDate)
ind.counts <- rbind(ind.counts, Totals, Total_effort, Total_effort_earliestStartDate, Total_effort_latestEndDate)

# populate community matrices
for(i in sites){

	# community matrix based on captures
	tmp <- caps.all[caps.all$siteID == i, ]
	tmp <- tmp[grep(" ", tmp$scientificName), ]
 	cap.counts[names(table(tmp$scientificName)), i] <- table(tmp$scientificName)
	cap.counts['Totals', i] <- nrow(tmp)
	rm(tmp)
  
	# community matrix based on individuals
	tmp <- inds.all[inds.all$siteID == i, ]
	tmp <- tmp[grep(" ", tmp$scientificName), ]
	ind.counts[names(table(tmp$scientificName)), i] <- table(tmp$scientificName)
	ind.counts['Totals', i] <- nrow(tmp)
	rm(tmp)

  }

# write to files
write.csv(cap.counts, '~/Library/CloudStorage/Box-Box/Projects/Chapman_GI_smammal-community-traits/_data_/NEON_Appalachians-smammals_total-captures_2023-02-28.csv', row.names = T, quote = F)
write.csv(ind.counts, '~/Library/CloudStorage/Box-Box/Projects/Chapman_GI_smammal-community-traits/_data_/NEON_Appalachians-smammals_total-individuals_2023-02-28.csv', row.names = T, quote = F)


##############################################
## create community (site x species) matrices,
## controlled for temporal differences in 
## effort among sites and years
##############################################

# get unique sites and spp from the capture data
sites <- unique(caps.all$siteID)
sites <- sites[order(sites)]
spp <- unique(caps.all$scientificName)[grep(" ", unique(caps.all$scientificName))] # keep only species with some taxonomy data
spp <- spp[order(spp)]

# create matrices to hold total captures and total individuals
cap.counts.timeControl <-  ind.counts.timeControl <- matrix(0, nrow = length(spp), ncol = length(sites), dimnames = list(spp, sites))
Totals <- rep(NA, length(sites)) # row to hold totals
Total_effort_earliestStartDateControlled <- Total_effort_latestEndDateControlled <- rep(NA, length(Total_effort_latestEndDate))
cap.counts.timeControl <- rbind(cap.counts.timeControl, Totals, Total_effort, Total_effort_earliestStartDateControlled, Total_effort_latestEndDateControlled)
ind.counts.timeControl <- rbind(ind.counts.timeControl, Totals, Total_effort, Total_effort_earliestStartDateControlled, Total_effort_latestEndDateControlled)

# calculate inner temporal window of sampling across sites and years
latestStartDate <- max(Total_effort_earliestStartDate) # the latest start date across years at each site
earliestEndDate <- min(Total_effort_latestEndDate) # the earliest end date across years at each site

# populate community matrices controlled for duration of effort across sites
for(i in sites){

	# community matrix based on captures
	tmp <- caps.all[caps.all$siteID == i, ]
	tmp <- tmp[grep(" ", tmp$scientificName), ]
	tmp <- tmp[tmp$julianday >= latestStartDate & tmp$julianday <= earliestEndDate, ]
 	cap.counts.timeControl[names(table(tmp$scientificName)), i] <- table(tmp$scientificName)
	cap.counts.timeControl['Totals', i] <- nrow(tmp)
	cap.counts.timeControl['Total_effort_earliestStartDateControlled', i] <- min(tmp$julianday) # the controlled start date
	cap.counts.timeControl['Total_effort_latestEndDateControlled', i] <- max(tmp$julianday) # the controlled end date
	rm(tmp)
  
	# community matrix based on individuals
	tmp <- inds.all[inds.all$siteID == i, ]
	tmp <- tmp[grep(" ", tmp$scientificName), ]
	tmp <- tmp[tmp$julianday >= latestStartDate & tmp$julianday <= earliestEndDate, ]
	ind.counts.timeControl[names(table(tmp$scientificName)), i] <- table(tmp$scientificName)
	ind.counts.timeControl['Totals', i] <- nrow(tmp)  
	ind.counts.timeControl['Total_effort_earliestStartDateControlled', i] <- min(tmp$julianday) # the controlled start date
	ind.counts.timeControl['Total_effort_latestEndDateControlled', i] <- max(tmp$julianday) # the controlled end date
	rm(tmp)

  }

# last step is to modify total trapnights (Total_effort) based on the controlled temporal window
setwd(top.dir)
for(i in sites){
	for(j in data.dirs[grep(i,data.dirs)]){

		# read file-wise
		# looping through sites and site dirs
 		setwd(j)
 		data.file <- list.files(pattern = 'pertrapnight')
 		tmp <- fread(data.file)
 		if(length(grep("1 -", tmp$trapStatus)) == nrow(tmp)) {setwd(top.dir); next} # skip files that have no trap effort, as above
		if(length(grep("1 -", tmp$trapStatus)) > 0) {tmp <- tmp[-grep("1 -", tmp$trapStatus), ]}  # limit to actual trapnights (i.e., excluding traps not set), as above
	
		# add Julian days to the record set
		julianday <- rep(NA, nrow(tmp))
		for(k in 1:nrow(tmp)){julianday[k] <- strptime(tmp$collectDate[k],"%Y-%m-%d")$yday} #convert to julian day
		tmp <- cbind(tmp, julianday)

		# record trap effort, timing, captures, individuals
		removes <- c(which(tmp$julianday < latestStartDate), which(tmp$julianday > earliestEndDate)) # trap effort falling outside inner sampling window
		cap.counts.timeControl['Total_effort', i] <- cap.counts.timeControl['Total_effort', i] - length(removes) # subtract from previous total
		ind.counts.timeControl['Total_effort', i] <- ind.counts.timeControl['Total_effort', i] - length(removes) # subtract from previous total

		# update on progress
		print(paste('read record set ', i, sep = ''))
		print(paste('removed ', length(removes), ' trap events based on time controlling', sep = ''))
		Sys.sleep(0.0000001)
  		rm(tmp)
 		
		# loop up to next record set/dir
		setwd(top.dir)
		}
	}

# write to files
write.csv(cap.counts.timeControl, '~/Library/CloudStorage/Box-Box/Projects/Chapman_GI_smammal-community-traits/_data_/NEON_Appalachians-smammals_total-captures-timeControl_2023-02-28.csv', row.names = T, quote = F)
write.csv(ind.counts.timeControl, '~/Library/CloudStorage/Box-Box/Projects/Chapman_GI_smammal-community-traits/_data_/NEON_Appalachians-smammals_total-individuals-timeControl_2023-02-28.csv', row.names = T, quote = F)

