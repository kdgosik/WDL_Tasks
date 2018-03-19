workflow cellranger {
    String sampleId
    Array[File] fastqs
    String referenceName
    File transcriptomeTarGz
    Boolean? secondary
    Int? expectCells
    String diskSpace

    call CellRanger {
        input:
        sampleId = sampleId,
        fastqs = fastqs,
        reference = referenceName,
        transcriptomeTarGz = transcriptomeTarGz,
        secondary = secondary,
        expectCells = expectCells,
        diskSpace = diskSpace
   }
}

task CellRanger {
    String sampleId
    Array[File] fastqs
    String reference
    File transcriptomeTarGz
    Int? expectCells
    Boolean? secondary
    String diskSpace

    command {
        set -e
        mkdir transcriptome_dir
        tar xf ${transcriptomeTarGz} -C transcriptome_dir --strip-components 1
        ln -s /usr/bin/python3 python
        export PATH=$PATH:.
        python <<CODE
        import os
        from subprocess import call
        dirs = dict()
        for f in ["${sep='","' fastqs}"]:
            dirs.setdefault(os.path.dirname(f), True)
        expect_cells = '${expectCells}'
        secondary = '${secondary}'
        call_args = list()
        call_args.append('cellranger')
        call_args.append('count')
        call_args.append('--jobmode=local')
        call_args.append('--transcriptome=transcriptome_dir')
        call_args.append('--id=' + '${sampleId}')
        call_args.append('--fastqs=' + ','.join(list(dirs.keys())))
        if secondary is not 'true':
            call_args.append('--nosecondary')
        if expect_cells is not '':
            call_args.append('--expect-cells=' + str(expect_cells))
        call(call_args)
        CODE
        }
    output {
        File barcodes = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/barcodes.tsv"
        File genes = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/genes.tsv"
        File matrix = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/matrix.mtx"
        File qc = "${sampleId}/outs/metrics_summary.csv"
        File report = "${sampleId}/outs/web_summary.html"
    }

    runtime {
        docker: "regevlab/cellranger-2.0.2"
        memory: "416 GB"
        bootDiskSizeGb: 12
        disks: "local-disk ${diskSpace} HDD"
        cpu: 64
        preemptible: 2
    }
}
