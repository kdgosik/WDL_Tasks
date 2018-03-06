task bundle_plots {
    Array[File] CNA_plots
    String sample_set_name

    command {
        python /bundle_helper.py --prefix ${sample_set_name} ${sep=' ' CNA_plots}
        pdflatex ${sample_set_name}.tex
    }
    runtime {
        docker: "jnktsj/pdflatex:1"
        memory: "2 GB"
        cpu: 1
    }
    output {
        File bundledPDF = "${sample_set_name}.pdf"
    }
}

workflow ichorSummaryBundle {
    Array[File] CNA_plots
    String sample_set_name

    call bundle_plots {
        input: CNA_plots = CNA_plots,
               sample_set_name = sample_set_name
    }
    output {
        bundle_plots.*
    }
}