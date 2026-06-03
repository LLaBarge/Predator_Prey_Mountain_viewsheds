
# For each encounter, check whether an alarm call was recorded in the
# ad lib data near the time of closest approach. Alarm subtypes:
#   A  = generic alarm
#   CH = chirp (major predator response)
#   RA = rahow (adult male alarm, major predator, but also given by Vera)
#   CT = chutter (snake/human alarm)
#   RP = rraup (predator bird alarm)

ALARM_WINDOW_MIN <- 15
# Parse dates and times for ALL ad lib records (not just baboon obs)
adlib$date_clean_all <- sapply(adlib[[date_col]], parse_date_gps_informed,
                               gps_dates_ref = gps_dates)
adlib$date_clean_all <- as.Date(adlib$date_clean_all, origin = "1970-01-01")

adlib$time_clean <- sapply(adlib[[names(adlib)[grepl("TIME", names(adlib), ignore.case = TRUE)][1]]],
                           parse_scan_time)
adlib$datetime <- as.POSIXct(paste(adlib$date_clean_all, adlib$time_clean),
                             format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Filter to alarm vocalisations
alarm_subtypes <- c("A", "CH", "RA", "CT", "RP")
subtype_col <- names(adlib)[grepl("SUBTYPE", names(adlib), ignore.case = TRUE)][1]

alarms <- adlib %>%
  filter(toupper(trimws(.data[[subtype_col]])) %in% alarm_subtypes,
         !is.na(datetime))

print(table(trimws(toupper(alarms[[subtype_col]]))))


# For each encounter, find closest approach time and check for alarms
encounter_results$closest_approach_time <- as.POSIXct(NA)
encounter_results$alarm_present <- FALSE
encounter_results$alarm_subtypes_found <- ""
encounter_results$alarm_observations <- ""
encounter_results$n_alarms_in_window <- 0

for (i in 1:nrow(encounter_results)) {
  enc_id <- encounter_results$encounter_id[i]
  
  # Get all proximity obs for this encounter
  enc_prox <- proximity %>%
    filter(encounter_id == enc_id) %>%
    arrange(distance_m)
  
  if (nrow(enc_prox) == 0) next
  
  # Time of closest approach
  closest_time <- enc_prox$timestamp[1]
  encounter_results$closest_approach_time[i] <- closest_time
  
  # Search for alarms within the window
  window_start <- closest_time - ALARM_WINDOW_MIN * 60
  window_end   <- closest_time + ALARM_WINDOW_MIN * 60
  
  matching_alarms <- alarms %>%
    filter(datetime >= window_start & datetime <= window_end)
  
  if (nrow(matching_alarms) > 0) {
    encounter_results$alarm_present[i] <- TRUE
    encounter_results$n_alarms_in_window[i] <- nrow(matching_alarms)
    encounter_results$alarm_subtypes_found[i] <- paste(
      unique(trimws(toupper(matching_alarms[[subtype_col]]))),
      collapse = ", "
    )
    encounter_results$alarm_observations[i] <- paste(
      unique(trimws(matching_alarms[[obs_col]])),
      collapse = "; "
    )
  }
}


# total encounters with 350m
n_total <- nrow(encounter_results)
n_alarm <- sum(encounter_results$alarm_present)
n_silent <- n_total - n_alarm

n_silent/n_total
# Save updated encounter results
encounter_results_final <- encounter_results
write.csv(encounter_results_final, "encounter_behavioral_analysis_with_alarms.csv",
          row.names = FALSE)
saveRDS(encounter_results_final, "encounter_results_with_alarms.rds")
