#!/home/johannes/bin nextflow

//compares the GC content of a bunch of input fastas and a reference fasta
//and plots it

//+++++++++++++++++ SETUP++++++++++++++++++++++++
params.workdir = "/home/johannes/rdrive/PPG_SEQ_DATA-LICHTJ-SE00182/johannes/notebook/bcinerea_test"

//Folder containing fasta files of genomes
params.input = "${params.workdir}/input/*.fasta"

//Reference genome fasta file
params.reference = "${params.workdir}/reference/*.fasta"

params.outdir = "${params.workdir}/output"
//+++++++++++++++++++++++++++++++++++++++++++++++

//if you are working in a folder structure that has a separate
//folder for a reference genome leave everything as it is
input = Channel
.fromPath(params.input)
.map{[it.getParent(), it]}
.set{ genomes }

//if you just want to compare a bunch of sequences comment
//out the next four lines and change 'from sequences' to
//'from input' in the 'OcculterCut' process
sequences = Channel
.fromPath(params.reference)
.map{[it.getBaseName(), it]}
.concat(genomes)

process OcculterCut {
  tag { id }

  input:
  set id, "genome.fasta" from sequences

  output:
  file "composition" into allCompositions

  """
/opt/occultercut/current/OcculterCut -f genome.fasta
paste compositionGC.txt <(yes "$id" | head -n \$(cat compositionGC.txt | wc -l)) > composition
"""
}

process MakeFigure {
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file "input.*.txt" from allCompositions.toList()

  output:
  file "gc_composition.svg"

  """
#!/usr/bin/env Rscript

library(magrittr)
library(ggplot2)
library(dplyr)

Sys.setenv("DISPLAY"=":0.0")
list.files(".", pattern="input.*.txt") %>%
  lapply(function(filename) {read.table(filename, col.names=c("gcperc", "genomeperc", "strain"))}) %>%
  do.call(what=rbind) %>%
  ggplot(aes(x=gcperc, y=genomeperc * 100)) +
  geom_line() +
  facet_grid(strain ~ ., scale="free") +
  theme_minimal() +
  xlab("GC (%)") +
  ylab("Genome representation (%)") -> plot
ggsave("gc_composition.svg", plot=plot)
  """
}  