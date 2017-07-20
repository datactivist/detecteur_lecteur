# Cleaning function
formatlog <- function(file) {
  newlog <- read_log(file)
  colnames(newlog) <- c("clientip" ,"V2", "V3", "timestamp", "doc", "status", "size", "refer", "agent")
  newlog$doc <- gsub(".org/.*",".org",newlog$doc) 
  counts <- newlog %>% 
    group_by(timestamp, doc)%>% 
    summarise(countconsult = n()) 
  counts$timestamp <- gsub(" \\+0100", "", counts$timestamp)
  counts$timestamp <- gsub("(\\d\\d:\\d\\d:\\d\\d):\\d\\d", "\\1:00", counts$timestamp)
  counts$timestamp <- dmy_hms(counts$timestamp)
  counts <- counts %>% filter(!grepl("PURGE", doc)) %>%
    group_by(timestamp)%>% 
    summarise(countconsult = sum(countconsult))
  return(counts)
}
# Joining several days
j0 <- formatlog("02272017-TOTAL-JOURNALS.log")
j1 <- formatlog("02282017-TOTAL-JOURNALS.log")
j2 <- formatlog("03012017-TOTAL-JOURNALS.log")

counts_fin <- rbind(j0,j1,j2)

#Plotting data 
res <- AnomalyDetectionTs(counts_fin, max_anoms=0.02, direction='both', only_last="day", plot=TRUE)
res$plot
res$anoms

save(res, file="output.RData")
save(counts_fin, file="input.RData")
