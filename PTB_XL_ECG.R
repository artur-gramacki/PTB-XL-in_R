# PTB_XL_ECG.R
# copyright (C) 2026 Artur Gramacki and Jarosław Gramacki
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#
# What is this script?
# --------------------
# This script is a port of the Python script example_physionet.py to R, 
# available at https://physionet.org/content/ptb-xl/1.0.3/
#

pkgs <- c("EGM", "jsonlite")

to_install = !pkgs %in% installed.packages()
if(any(to_install)) {
	install.packages(pkgs[to_install])
}

library("EGM")
library("jsonlite")

sampling_rate <- 500
number_of_recors <- 21799
ecg_length <- 10
lead_12 <- c("I", "II", "III", "AVR", "AVL", "AVF", "V1", "V2", "V3", "V4", "V5", "V6" )

# Range for displaying a sample ECG (in sec.)
start <- 0
stop <- 10
# Which ECG to display
ecg_id <- 146

# //////////////////////////////////////////////////////////////////////////////
# Load and convert annotation data ----
# //////////////////////////////////////////////////////////////////////////////
Y <- read.csv(file = "ptbxl_database.csv", header = TRUE, sep = ",")

# //////////////////////////////////////////////////////////////////////////////
# Load scp_statements.csv for diagnostic aggregation ----
# //////////////////////////////////////////////////////////////////////////////
agg_df <- read.csv(file = "scp_statements.csv", header = TRUE, sep = ",")
colnames(agg_df)[1] <- "index"
agg_df <- agg_df[which(agg_df$diagnostic == 1),]

# //////////////////////////////////////////////////////////////////////////////
# Apply diagnostic superclass ----
# //////////////////////////////////////////////////////////////////////////////
for (i in 1:nrow(Y)) {
	if (i %% 1000 == 0) cat(i, "/", nrow(Y), "\n")
	y <- Y$scp_codes[i]
	y_json <- gsub("'", "\"", y)
	obj <- fromJSON(y_json)
	tmp <- c()
	for (j in 1:length(obj)) {
		obj_name <- names(obj[j])
		idx <- which(agg_df$index == obj_name)
		if (length(idx) > 0) {
			dc <- agg_df$diagnostic_class[idx]
			tmp <- c(tmp, dc)
		}
	}
	tmp <- unique(tmp)
	tmp <- sort(tmp)
	tmp <- paste(tmp, collapse = " | ", sep = "")
	Y$diagnostic_superclass[i] <- tmp
}

# //////////////////////////////////////////////////////////////////////////////
# Load raw signal data ----
# //////////////////////////////////////////////////////////////////////////////
# dim = c(rows, columns, layers, ...)
# data is stored colum-wise (column-major order)
X <- array(
	data = rep(NA, nrow(Y)), 
	dim = c(sampling_rate * ecg_length, 12, nrow(Y))
)

dimnames(X) <- list(
	c(formatC(seq(0, ecg_length, length.out = sampling_rate * ecg_length), format = "f", digits = 3)),  
	seq(1:12), 
	seq(1:number_of_recors)
)

for (i in 1:nrow(Y)) {
	if (i %% 100 == 0) cat(i, "/", nrow(Y), "\n")
	if (sampling_rate == 500) {
		file <- Y$filename_hr[i]
	} 	
	if (sampling_rate == 100) {
		file <- Y$filename_lr[i]
	} 	
	out <- read_wfdb(
		record = basename(file),
		record_dir = dirname(file), 
		units = "physical"
	)
	X[, , i]	<- as.matrix(out$signal[,2:13])
	dimnames(X)[[2]] <- out$header$lead
	dimnames(X)[[3]][[i]] <- basename(file)
}

if (sampling_rate == 500) {
	save(X, file = "ECG_records500.RData")
} 	
if (sampling_rate == 100) {
	save(X, file = "ECG_records100.RData")
} 		

# load("ECG_records100.RData")
# load("ECG_records500.RData")

# //////////////////////////////////////////////////////////////////////////////
# Sample prints ----
# //////////////////////////////////////////////////////////////////////////////
head(X[, , 1], 10)
tail(X[, , 1], 10)

head(X[, , number_of_recors], 10)
tail(X[, , number_of_recors], 10)

# //////////////////////////////////////////////////////////////////////////////
# Split data into train and test ----
# //////////////////////////////////////////////////////////////////////////////
test_fold = 10

# Train
rows_train <- which(Y$strat_fold != test_fold)
X_train <- X[, , rows_train]
y_train = Y[rows_train, ]$diagnostic_superclass

# Test
rows_test <- which(Y$strat_fold == test_fold)
X_test <- X[, , rows_test]
y_test = Y[rows_test, ]$diagnostic_superclass

# //////////////////////////////////////////////////////////////////////////////
# ECG sample plots ----
# //////////////////////////////////////////////////////////////////////////////
from <- start * sampling_rate
to <- stop * sampling_rate
range <- seq(from, to)
row <- which(Y$ecg_id == ecg_id)

par(mfrow = c(13, 1), pty = "m")
# mai: c(bottom, left, top, right)
par(mai = c(0.05, 0.7, 0.1, 0.2))
for (j in 1:12) {
	plot(X[range, j, row], type = "l", ylab = lead_12[j], las = 1, xlab = "", xaxt = "n")
}

lab <- seq(from = dimnames(X)[[1]][from + 1], to = dimnames(X)[[1]][to], length.out = 11)
axis(
	side = 1, las = 1, cex.axis = 0.9,
	at = seq(0, length(range), length.out = 11),
	labels = c(formatC(lab, format = "f", digits = 2))
)

par(mai = c(0.05, 0.7, 0.2, 0.2))
plot(1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
txt <- paste("ECG ID: ", ecg_id, ", sampling rate: ", sampling_rate, "Hz", 
    				 ", diag: ", Y$diagnostic_superclass[row], sep = "")
text(1,1, txt, cex = 1.5, col = "blue")

