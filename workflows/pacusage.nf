include { VALIDATE_INPUTS } from '../modules/local/validate_inputs'
include { RECORD_INPUT_CHECKSUMS } from '../modules/local/record_input_checksums'
include { RECORD_SOFTWARE_VERSIONS } from '../modules/local/record_software_versions'
include { BUILD_REPORT } from '../modules/local/build_report'
include { PREPARATION } from '../subworkflows/local/preparation'
include { DISCOVERY } from '../subworkflows/local/discovery'
include { QUANTIFICATION } from '../subworkflows/local/quantification'
include { STATISTICS } from '../subworkflows/local/statistics'

workflow PACUSAGE {
    main:
    if (!params.input || !params.fasta || !params.gtf || !params.assembly) {
        error "Required parameters: --input, --assembly, --fasta, and --gtf"
    }

    // File parameters are resolved against the launch directory here, because
    // VALIDATE_INPUTS runs in its own task directory.
    def pathParameters = [
        'input', 'fasta', 'gtf', 'known_pacs', 'chromosome_aliases', 'pas_motif_catalog'
    ]
    def parameterMap = [:]
    params.each { key, value ->
        parameterMap[key] = key in pathParameters && value ? file(value).toString() : value
    }
    def encodedParams = groovy.json.JsonOutput.toJson(parameterMap)
        .bytes
        .encodeBase64()
        .toString()
    // Staged as an input so that -resume reruns the statistics whenever the
    // script changes.
    def statistics_script = file("${projectDir}/scripts/fit_usage_model.R", checkIfExists: true)

    VALIDATE_INPUTS(
        file(params.input, checkIfExists: true),
        encodedParams,
        file("${projectDir}/nextflow_schema.json", checkIfExists: true)
    )
    RECORD_SOFTWARE_VERSIONS(VALIDATE_INPUTS.out.software_versions, statistics_script)
    normalized = VALIDATE_INPUTS.out.normalized_samples
    resolved = VALIDATE_INPUTS.out.resolved_params
    annotation = Channel.value(file(params.gtf, checkIfExists: true))

    PREPARATION(normalized, annotation)
    RECORD_INPUT_CHECKSUMS(
        VALIDATE_INPUTS.out.input_checksums,
        normalized,
        PREPARATION.out.alignment_qc
    )
    DISCOVERY(
        PREPARATION.out.samples,
        normalized,
        PREPARATION.out.reference,
        annotation,
        resolved
    )
    QUANTIFICATION(
        DISCOVERY.out.evidence,
        DISCOVERY.out.atlas,
        DISCOVERY.out.run_resolution,
        DISCOVERY.out.kernel,
        resolved
    )
    STATISTICS(
        normalized,
        QUANTIFICATION.out.wide,
        QUANTIFICATION.out.pau,
        DISCOVERY.out.atlas,
        resolved,
        statistics_script
    )

    report_artifacts = Channel.empty()
        .mix(normalized)
        .mix(VALIDATE_INPUTS.out.input_validation)
        .mix(VALIDATE_INPUTS.out.control_mapping)
        .mix(PREPARATION.out.reference_qc)
        .mix(PREPARATION.out.alignment_qc.flatten())
        .mix(PREPARATION.out.strandedness_qc.flatten())
        .mix(DISCOVERY.out.calibration)
        .mix(DISCOVERY.out.kernel_diagnostics)
        .mix(DISCOVERY.out.filtering_qc.flatten())
        .mix(DISCOVERY.out.discovery_qc)
        .mix(DISCOVERY.out.atlas)
        .mix(QUANTIFICATION.out.pau)
        .mix(QUANTIFICATION.out.quantification_qc.flatten())
        .mix(STATISTICS.out.motif_scores)
        .mix(STATISTICS.out.motif_sensitivity)
        .mix(STATISTICS.out.results)
        .mix(STATISTICS.out.kmer_results)
        .collect()
    BUILD_REPORT(report_artifacts, resolved)
}
