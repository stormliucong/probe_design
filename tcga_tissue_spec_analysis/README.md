Select tissue specific probes from TCGA 450K data.

betaMatrixGenerator.pl generates sequence of beta Matrix to reduce RAM burden.
betaMatrixGeneratorBugCorrection.pl corrects a bug for the sequnce of betaMatrix by adding a 99 file.
wholeBloodMethylGEOPreprocessing.R collects 5 blood 450K methylation samples and create a big beta Matrix.
tissueSpecCpG.R conducts M-value and t-test for each CpG in a group v.s. others manner.
pickOutCpGAndGene.R pick out tissue-specific genes based on certain criteria.
