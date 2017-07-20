#Function to clean the logs as little as possible before getting the journal name

cleanlog <- function(namelog) {
  logs <- read_log(namelog)
  colnames(logs) <- c("clientip", "V2", "V3", "timestamp", "doc", "status", "size", "refer", "agent")
  logs <- tbl_df(logs)
  logs$doc <- gsub("org/.*", "org", logs$doc) #Only keep the .org DNS
  logs <- logs %>% filter(!grepl("PURGE", doc)) #We remove the "PURGE" command (that do not imply any consultation)
  return(logs)
}

#Function to display Anomaly Detection function for each revue in Open Edition logs (with some cleaning implied)

get_revue_anomalies <- function(full_consultation, revues) {
  df_anoms <- data.frame(matrix(ncol = 3, nrow = 0)) #We initialize an empty dataframe to contain the anoms
  colnames(df_anoms) <- c("timestamp", "freq", "revue")
  for (revue in revues) {
    counts_revue <- full_consultation %>% filter(doc == revue) #Get only a subsample of the logs for the revue (essential to get quickly in order to decrease computation)
    counts_revue$timestamp <- gsub(" \\+0100", "", counts_revue$timestamp) #Remove time zone
    counts_revue$timestamp <- gsub("(\\d\\d:\\d\\d:\\d\\d):\\d\\d", "\\1:00", counts_revue$timestamp) #Remove second precision (with regex)
    counts_revue$timestamp <- dmy_hms(counts_revue$timestamp) #Get date format
    counts_revue <- counts_revue %>% group_by(timestamp) %>% summarise(count_consult=n()) #merge to each minute to get the freq counts
    revue_title = gsub("GET http://", "", revue) #clean revue name
    graph_title <- paste("Anomalies de consultation sur ", revue_title, sep=" ") #set the graph title with the proper revue name
    res_revue = AnomalyDetectionTs(counts_revue, max_anoms=0.02, direction='both', plot=TRUE, title=graph_title) #Get the Anomaly Detection
    show(res_revue$plot) #Show the graphic (but could be an if to save some computation)
    revue_anoms <- res_revue$anoms #initiate a new df with only anoms
    revue_anoms$revue <- revue_title #add the revue name
    df_anoms <- rbind(df_anoms, revue_anoms) #bind the anom df of this revue to the full df of each revues
    Sys.sleep(0) #add a sys sleep (in order to properly display the graphs and force update)
  }
  return(df_anoms)
}

get_one_revue_anomalies <- function(full_consultation, revue) {
    counts_revue <- full_consultation %>% filter(doc == revue)
    counts_revue$timestamp <- gsub(" \\+0100", "", counts_revue$timestamp)
    counts_revue$timestamp <- gsub("(\\d\\d:\\d\\d:\\d\\d):\\d\\d", "\\1:00", counts_revue$timestamp)
    counts_revue$timestamp <- dmy_hms(counts_revue$timestamp)
    counts_revue <- counts_revue %>% group_by(timestamp) %>% summarise(count_consult=n())
    revue_title = gsub("GET http://", "", revue)
    revue_title <- paste("Anomalies de consultation sur ", revue_title, sep=" ")
    #    res_revue = AnomalyDetectionVec(counts_revue[,2], max_anoms=0.49, period=length(counts_lectures), direction='both', plot=TRUE, title=revue_title)
    res_revue = AnomalyDetectionTs(counts_revue, max_anoms=0.02, direction='both', plot=TRUE, title=revue_title)
    show(res_revue$plot)
    return(res_revue$anoms)
}


revue_anom <- get_one_revue_anomalies(j_full, "GET http://polis.revues.org")

j3 <- cleanlog("02282017-TOTAL-JOURNALS.log")
j2 <- cleanlog("02272017-TOTAL-JOURNALS.log")
j1 <- cleanlog("02262017-TOTAL-JOURNALS.log")
j_full <- rbind(j1,j2,j3)
