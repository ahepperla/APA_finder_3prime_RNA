include { VALIDATE_INPUTS } from '../modules/local/validate_inputs'
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

    def parameterMap = [:]
    params.each { key, value -> parameterMap[key] = value }
    def encodedParams = groovy.json.JsonOutput.toJson(parameterMap)
        .bytes
        .encodeBase64()
        .toString()

    VALIDATE_INPUTS(
        Channel.fromPath(params.input, checkIfExists: true),
        encodedParams
    )
    normalized = VALIDATE_INPUTS.out.normalized_samples.first()
    resolved = VALIDATE_INPUTS.out.resolved_params.first()
    annotation = Channel.value(file(params.gtf, checkIfExists: true))

    PREPARATION(normalized, resolved, annotation)
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
        QUANTIFICATION.out.testable,
        QUANTIFICATION.out.pau,
        DISCOVERY.out.atlas,
        resolved
    )

    report_artifacts = Channel.empty()
        .mix(VALIDATE_INPUTS.out.input_validation)
        .mix(VALIDATE_INPUTS.out.control_mapping)
        .mix(PREPARATION.out.reference_qc)
        .mix(PREPARATION.out.alignment_qc.flatten())
        .mix(PREPARATION.out.strandedness_qc.flatten())
        .mix(DISCOVERY.out.calibration)
        .mix(DISCOVERY.out.filtering_qc.flatten())
        .mix(DISCOVERY.out.discovery_qc)
        .mix(DISCOVERY.out.atlas)
        .mix(QUANTIFICATION.out.pau)
        .mix(QUANTIFICATION.out.filtering)
        .mix(QUANTIFICATION.out.quantification_qc.flatten())
        .mix(STATISTICS.out.motif_scores)
        .mix(STATISTICS.out.motif_sensitivity)
        .mix(STATISTICS.out.results)
        .mix(STATISTICS.out.kmer_results)
        .collect()
    BUILD_REPORT(report_artifacts)
}
