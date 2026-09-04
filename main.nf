#!/usr/bin/env nextflow

include { PACUSAGE } from './workflows/pacusage'

def helpMessage() {
    def schema = new groovy.json.JsonSlurper()
        .parse(file("${projectDir}/nextflow_schema.json"))
    def lines = [
        'PACusage: polyadenylation site discovery and differential usage',
        '',
        'Usage:',
        '  nextflow run /path/to/pacusage -profile local,conda -params-file analysis.yaml -resume',
        ''
    ]
    schema.definitions.each { sectionName, group ->
        lines << "${group.title}:"
        group.properties.each { name, definition ->
            def defaultText = definition.containsKey('default') ?
                " [default: ${definition.default}]" : ''
            lines << "  --${name.padRight(38)} ${definition.description ?: ''}${defaultText}"
        }
        lines << ''
    }
    lines.join('\n')
}

workflow {
    if (params.help) {
        log.info helpMessage()
    } else {
        PACUSAGE()
    }
}
