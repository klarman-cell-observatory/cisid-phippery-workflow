/*
 * This Source Code Form is subject to the terms of the GNU GENERAL PUBLIC LICENCE
 * License, v. 3.0. 
 */


/* 
 * 'PhIP-Flow' - A Nextflow pipeline for running common phip-seq analysis workflows
 * 
 * Fred Hutchinson Cancer Research Center, Seattle WA.
 * 
 * Jared Galloway
 * Kevin Sung
 * Sam Minot
 * Erick Matsen 
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define the default parameters - example data get's run by default
 */ 
params.sample_table     = params.sample_table ?: "/phipflow_T/data/sample_table.csv"
params.peptide_table    = params.peptide_table ?: "/phipflow_T/data/peptide_table.csv"
params.reads_prefix     = params.reads_prefix ?: "/phipflow_T"
params.results          = "$PWD/results/"

log.info """\
P H I P - F L O W!
Matsen, Overbaugh, and Minot Labs
Fred Hutchinson CRC, Seattle WA
================================
sample_table    : $params.sample_table
peptide_table   : $params.peptide_table
results         : $params.results

"""

/* 
 * Import modules 
 */
nextflow.enable.dsl=2

include { ALIGN } from './workflows/alignment.nf'
include { STATS } from './workflows/statistics.nf'
include { DSOUT } from './workflows/output.nf'
include { AGG } from './workflows/aggregate.nf'

workflow {
    ALIGN | STATS | DSOUT | AGG
}
