task gmql {
  String expPath
  String refPath

  command <<<

  query="R = SELECT() %s; \
         E = SELECT() %s ; \
         O = JOIN( DLE(1000) ) R E;\
         MATERIALIZE O INTO file:///usr/src/myapp/results/;"

  printf "$query" ${refPath} ${expPath} > /usr/src/myapp/query.gmql
  
  cat  /usr/src/myapp/query.gmql

  java -jar /usr/src/myapp/uber-GMQL-Cli-2.0.jar -scriptpath /usr/src/myapp/query.gmql
  
  tar -zcvf out.tar.gz /usr/src/myapp/results/exp/
  
  >>>
  output {
    File response = stdout()
    File out = "out.tar.gz"
  }
  
  runtime {
        docker: "pp86/gmql"
        memory: "24 GB"
        disks: "local-disk 100 SSD"
    }
}



workflow ContEstMuTectOncotatorWorkflow {
  String expPath
  String refPath
  
  call gmql{
     input:
        expPath=expPath,
        refPath=refPath
  }

}