data <- read.csv("Data/june_26_sim_data.csv")

data <- data[data$num_studies < 11, ]

reml_data <- data[data$method == "REML", ]

# Looking at REML only
# Get data set for all ts
reml1 <- reml_data[reml_data$true_ts == 0, ]
reml2 <- reml_data[reml_data$true_ts == 0.008, ]
reml3 <- reml_data[reml_data$true_ts == 0.03, ]
reml4 <- reml_data[reml_data$true_ts == 0.16, ]

######################################################
# REML - How often we go from 0 to nonzero
######################################################
counter <- 0
for (i in 1:5000) {
  loop_data <- reml1[reml1$sim == i, ]
  for (j in 1:(nrow(loop_data) - 1)) {
    if (loop_data$est_ts[j] == 0 && loop_data$est_ts[j+1] != 0) {
      counter <- counter + 1
      break
    }
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml2[reml2$sim == i, ]
  for (j in 1:(nrow(loop_data) - 1)) {
    if (loop_data$est_ts[j] == 0 && loop_data$est_ts[j+1] != 0) {
      counter <- counter + 1
      break
    }
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml3[reml3$sim == i, ]
  for (j in 1:(nrow(loop_data) - 1)) {
    if (loop_data$est_ts[j] == 0 && loop_data$est_ts[j+1] != 0) {
      counter <- counter + 1
      break
    }
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml4[reml4$sim == i, ]
  for (j in 1:(nrow(loop_data) - 1)) {
    if (loop_data$est_ts[j] == 0 && loop_data$est_ts[j+1] != 0) {
      counter <- counter + 1
      break
    }
  }
}

######################################################
# REML - How often we go from nonzero to 0
######################################################

counter <- 0
for (i in 1:5000) {
  loop_data <- reml1[reml1$sim == i, ]
  for (j in 1:(nrow(loop_data) - 1)) {
    if (loop_data$est_ts[j] != 0 && loop_data$est_ts[j+1] == 0) {
      counter <- counter + 1
      break
    }
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml2[reml2$sim == i, ]
  for (j in 1:(nrow(loop_data) - 1)) {
    if (loop_data$est_ts[j] != 0 && loop_data$est_ts[j+1] == 0) {
      counter <- counter + 1
      break
    }
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml3[reml3$sim == i, ]
  for (j in 1:(nrow(loop_data) - 1)) {
    if (loop_data$est_ts[j] != 0 && loop_data$est_ts[j+1] == 0) {
      counter <- counter + 1
      break
    }
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml4[reml4$sim == i, ]
  for (j in 1:(nrow(loop_data) - 1)) {
    if (loop_data$est_ts[j] != 0 && loop_data$est_ts[j+1] == 0) {
      counter <- counter + 1
      break
    }
  }
}

######################################################
# REML - Percent of sim reps that never include a 0 estimate
######################################################

counter <- 0
for (i in 1:5000) {
  loop_data <- reml1[reml1$sim == i, ]
  if (0 %in% loop_data$est_ts) {
    counter <- counter + 0
  } else {
    counter <- counter + 1
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml2[reml2$sim == i, ]
  if (0 %in% loop_data$est_ts) {
    counter <- counter + 0
  } else {
    counter <- counter + 1
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml3[reml3$sim == i, ]
  if (0 %in% loop_data$est_ts) {
    counter <- counter + 0
  } else {
    counter <- counter + 1
  }
}

counter <- 0
for (i in 1:5000) {
  loop_data <- reml4[reml4$sim == i, ]
  if (0 %in% loop_data$est_ts) {
    counter <- counter + 0
  } else {
    counter <- counter + 1
  }
}






