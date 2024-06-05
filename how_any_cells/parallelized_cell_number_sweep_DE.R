suppressMessages(suppressWarnings(library(parallel)))


no_cores <- detectCores() - 1  # Leave one core free for system processes
cl       <- makeCluster(no_cores)

cell_type

results <- parLapply(cl, 1:10000, function(i) {
	
	
		

})

stopCluster(cl)
