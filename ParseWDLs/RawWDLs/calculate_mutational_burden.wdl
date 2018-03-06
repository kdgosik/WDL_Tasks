workflow calculateMutationalBurden {
    String patientId
    File snvHandle
    File? indelHandle
    File coverageHandle

    call calculateMutationalBurdenTask {
        input: patientId=patientId,
            snvHandle=snvHandle,
            indelHandle=indelHandle,
            coverageHandle=coverageHandle
    }
    
    output  {
        calculateMutationalBurdenTask.mutational_burden
    }
}

task calculateMutationalBurdenTask {
    String patientId
    File snvHandle
    File? indelHandle
    File coverageHandle

    command {
        python /home/calc_mutburden.py -patient_id ${patientId} -snv ${snvHandle} ${"-indel " + indelHandle} -coverage ${coverageHandle}
    }

    output  {
        File mutational_burden="${patientId}.mutational_burden.txt"
    }

    runtime {
        docker: "breardon/calc_mutational_burden:1.1"
        memory: "1 GB"
    }
}
