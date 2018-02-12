setwd("/Users/kgosik/Documents/Projects/WebApps/GraphQLPractice/GraphQL-js/WDL_Tasks")


m <- readLines("methods.txt")

m[1:20]

methods_data <- data.frame(matrix(NA, ncol = 5, nrow = 271))
colnames(methods_data) <- m[1:5]

i <- 6
j <- 1
k <- 1
beginning <- TRUE
repeat{
  
  if( j == 1 & !{beginning} ) i <- i + 1 #next
  
  if(j == 1) {
    if(beginning) {
      methods_data[k, j] <- paste0(c(m[i], m[(i+1)]), collapse = ":")
      beginning <- FALSE
    }else{
      j <- j + 1
    }
  }
  
  if(j == 2){
    if(grepl("@", m[i])) {
      methods_data[k, j] <- ""
      j <- j + 1
    }else{
      if(is.na(methods_data[k, j])) {
        methods_data[k, j] <- m[i]
      }else{
        j <- j + 1
      }
    }
  }
  
  if(j == 3){
    if(is.na(methods_data[k, j])) {
      methods_data[k, j] <- m[i]
    }else{
      j <- j + 1
    }
  }
  
  if(j == 4){
    if(is.na(methods_data[k, j])) {
      methods_data[k, j] <- m[i]
    }else{
      j <- j + 1
    }
  }
  
  if(j == 5){
      methods_data[k, j] <- m[i]
      j <- 1
      k <- k + 1
      beginning <- TRUE
  } 
  i <- i + 1
  i
  j
  beginning
  head(methods_data)
  
  if( i == 1553 ) break
}
 
 
 
 ##https://api.firecloud.org/ga4gh/v1/tools/broadinstitute_cga:WXS_hg19_MutationCalling_QC_v1-3_BETA/versions/1/plain-WDL/descriptor
 
 
 out <- paste0("https://api.firecloud.org/ga4gh/v1/tools/",
        methods_data$Method,
        "/versions/",
        methods_data$Snapshots,
        "/plain-WDL/descriptor")
 
 
 writeLines(out, "Methods_APIs.txt")
 