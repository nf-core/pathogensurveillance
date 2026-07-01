#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/pathogensurveillance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/pathogensurveillance
    Website: https://nf-co.re/pathogensurveillance
    Slack  : https://nfcore.slack.com/channels/pathogensurveillance
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NFCORE_PATHOGENSURVEILLANCE  } from './workflows/pathogensurveillance'
include { PIPELINE_INITIALISATION      } from './subworkflows/local/utils_nfcore_pathogensurveillance_pipeline'
include { PIPELINE_COMPLETION          } from './subworkflows/local/utils_nfcore_pathogensurveillance_pipeline'
include { getGenomeAttribute           } from './subworkflows/local/utils_nfcore_pathogensurveillance_pipeline'
include { completionEmail              } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { completionSummary            } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { imNotification               } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap             } from 'plugin/nf-schema'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        params.reference_data
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PATHOGENSURVEILLANCE (
        PIPELINE_INITIALISATION.out.sample_data_tsv,
        PIPELINE_INITIALISATION.out.reference_data_tsv
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PATHOGENSURVEILLANCE.out.multiqc_report
    )

    ch_summary_params = PIPELINE_COMPLETION.out.summary_params
    ch_multiqc_report_list = PIPELINE_COMPLETION.out.multiqc_report_list

    onComplete:
    if (params.email || params.email_on_fail) {
        completionEmail(
            ch_summary_params.val,
            params.email,
            params.email_on_fail,
            params.plaintext_email,
            params.outdir,
            params.monochrome_logs,
            ch_multiqc_report_list.getVal(),
        )
    }
    completionSummary(params.monochrome_logs)
    if (params.hook_url) {
        imNotification(ch_summary_params.val, params.hook_url)
    }

    onError:
    log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
