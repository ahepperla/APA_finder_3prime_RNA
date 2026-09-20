nextflow.enable.dsl = 2

params.batch_root = null

process ASSERT_SINGLE_BOOTSTRAP_BATCH {
    tag "${family}:${batch.baseName}"

    input:
    tuple val(family), path(batch)

    output:
    path '*.seen'

    script:
    """
    touch '${batch.baseName}.seen'
    """
}

workflow {
    if (!params.batch_root) {
        error "Set --batch_root to a directory containing bootstrap batch files."
    }

    batched_output = Channel
        .fromPath("${params.batch_root}/*.rds", checkIfExists: true)
        .collect()
        .map { batch_files -> tuple('test-family', batch_files) }

    flattened_batches = batched_output
        .flatMap { family, batch_files ->
            def files = batch_files instanceof List ? batch_files : [batch_files]
            files.collect { batch -> tuple(family, batch) }
        }

    ASSERT_SINGLE_BOOTSTRAP_BATCH(flattened_batches)
}
