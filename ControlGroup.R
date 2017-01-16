
# create control group patients

# read data

# items and labvalue
chart <- read.csv("chartevents.csv", header= T)
# lactate and other items
lab<- read.csv("labevents.csv", header= T)  
# blood culture
microbio <- read.csv("microbiologyevents.csv", header= T)  

# lab events - hl.bc, hl.only, and bc.only
length(table(lab$subject_id))  # 13753 patients with labs taken

# lactate lab values
lab.lactate <- lab[lab$itemid == 50010,]
length(table(lab.lactate$subject_id))  # 10697 total people had lactate tested

# high lactate lab values
lab.hl <- lab.lactate[lab.lactate$meanvalue >= 4,]  
length(table(lab.hl$subject_id))  # 3274 total people with high lactate values
length(table(micro$subject_id))  # 8930 total people with blood cultures taken

# identify individuals with a high lactate
hl.patients <- lab.hl[!duplicated(lab.hl$subject_id) & !is.na(lab.hl$subject_id),]  
# identify individuals with a blood culture taken
bc.patients <- micro[!duplicated(micro$subject_id),]  

# as vectors
hl <- as.vector(hl.patients$subject_id)
bc <- as.vector(bc.patients$subject_id)

# intersection and difference
hl.bc<- intersect(hl, bc)  #  2646 patients with blood culture and high lactate
hl.only <- setdiff(hl, bc)  #  628 patients with high lactate only
bc.only <- setdiff(bc, hl)  #  6278 patients with blood culture only (Will read 6284 because of 6 patients have blood cultures, but no labs)

# rearrange and merge the data sets
micro.1 <- data.frame(micro$subject_id, micro$halfhour, micro$spec_itemid)
micro.1$itemid <- NA
micro.1$meanvalue <- NA
micro.1$stdvalue <- NA
names(micro.1) <- c("SUBJECT_ID", "HALFHOUR", "SPEC_ITEMID", "ITEMID", "MEANVALUE", "STDVALUE")

lab.1 <- data.frame(lab$subject_id, lab$halfhour, lab$itemid, lab$meanvalue, lab$stdvalue)
lab.1$spec_itemid <- NA
names(lab.1) <- c("SUBJECT_ID", "HALFHOUR", "ITEMID", "MEANVALUE", "STDVALUE", "SPEC_ITEMID")

complete <- rbind(micro.1, lab.1)

# remove rows where CHARTTIME is missing
y <- complete[complete$HALFHOUR!="",]

# Ccount the number of people in the data set  
length(table(y$SUBJECT_ID)) # 13759

# control group - hl, no bc w/in 48 hr period; no hl, bc w/in 48 hr period; neither hl or bc w/in 48 hr period
control.group <- complete[!complete$SUBJECT_ID %in% hl.bc, ]

library(plyr)
control.group <- control.group[order(control.group$SUBJECT_ID, control.group$HALFHOUR),]

# Convert HALFHOUR to Date/Time format
library(lubridate)
control.group$HALFHOUR <- ymd_hms(control.group$HALFHOUR)

# start times and end times for each patient
start.time <- ddply(control.group, "SUBJECT_ID", summarise, min(HALFHOUR))
end.time <- ddply(control.group, "SUBJECT_ID", summarise, max(HALFHOUR))

start.time2 <- start.time[,2]+86400 # 86400 seconds = 24 hours
end.time2 <- end.time[,2]-86400

# control.group2 <- control.group[control.group$HALFHOUR >= start.time2 & control.group$HALFHOUR <= end.time2,]

pool <- data.frame(start.time[,1],start.time2,end.time2)
names(pool) <- c("SUBJECT_ID", "start.time", "end.time")
pool <- pool[!is.na(pool$start.time) | !is.na(pool$end.time),]  

pool$diff <- pool$end.time - pool$start.time
set.seed(12345)
random.time <- sapply(1:nrow(pool), function(i) runif(1, 0, max = pool[i,4]))
pool$middle.time <- pool$start.time+as.period(random.time,unit="seconds")
head(pool)

minutes <- floor( as.numeric( substr( pool$middle.time, 15, 16 ) ) / 30 ) * 30
minutes <- ifelse( minutes == 30, "30", "00" )
seconds <- rep ("00", nrow(pool))

pool$middle.time <- as.character( pool$middle.time )
substr( pool$middle.time , 15, 16 ) <- minutes
substr( pool$middle.time , 18, 19 ) <- seconds
pool$middle.time <- as.POSIXct( pool$middle.time,format="%Y-%m-%d %H:%M:%S")
head( pool)

pool$left.time <- pool$middle.time - 86400
pool$right.time <- pool$middle.time + 86400

head(pool)
write.csv(pool, "pool.csv")

pool2 <- pool[,c(1,5)]
head(pool2)
write.csv(pool2, "pool2.csv")
