task gmql {
  # A WDL task representing the GMQL query
  
  String expPath	# reference to the experiment dataset
  String refPath	# reference to the reference dataset

  command <<<
  
  # In order to use the GMQL application you need to provide to it a file containing the query in a textual format.
  # Therefore what is done in the next lines of code consists in:
  # 1) Using a variable 'query' to hold the string of the query. The string has some placeholders %s for making it customizable with the paths of the datasets.
  # 2) We use the printf function to assign to each placeholder the dataset path
  # 3) we save the resulting string in a predefined location (/usr/src/myapp/query.gmql)
  # 4) Finally we launch the GMQL application with java -jar /usr/src/myapp/uber-GMQL-Cli-2.0.jar -scriptpath /usr/src/myapp/query.gmql
  # 5) The GMQL application will create a folder with the results (/usr/src/myapp/results/)
  # 6) We zip the results and output them
  
  # The query is defined as a string with some parameters (%s)
  # Follows the query explanation:
  # Put in variable R the data in the first dataset
  # Put in variable E the data in the second dataset
  # Perform a JOIN operation between R and E taking only the regions with distance less than 1000
  # Save the results in a staging area for later processing
  query="R = SELECT() %s; \
         E = SELECT() %s; \
         O = JOIN( DLE(1000) ) R E;\
         MATERIALIZE O INTO file:///usr/src/myapp/results/;"

  printf "$query" ${refPath} ${expPath} > /usr/src/myapp/query.gmql	# We use the printf function to associate the variables %s to each dataset path
  
  # We launch the execution using the GMQL jar and the specified options
  java -jar /usr/src/myapp/uber-GMQL-Cli-2.0.jar -scriptpath /usr/src/myapp/query.gmql
  
  # Finally we zip the result folder 
  tar -zcvf out.tar.gz /usr/src/myapp/results/exp/
  
  >>>
  output {
    File response = stdout()
    File out = "out.tar.gz"		# The output is the location of the result folder
  }
  
  runtime {
        docker: "pp86/gmql"
        memory: "24 GB"
        disks: "local-disk 100 SSD"
    }
}



workflow GMQL_Tutorial_1_A_simple_query {
  String expPath
  String refPath
  
  call gmql{
     input:
        expPath=expPath,
        refPath=refPath
  }

}