A new paradigm for personalized medicine is emerging in the biomedical sciences: quantified wellness. Rather than focus on treatment of disease by dissecting mechanisms and etiologies, this approach aims to monitor and quantify patient "wellness" by longitudinal tracking with a variety of physiological and molecular data [1]. This exciting new approach to medicine, however, raises a number of basic questions related to data integration that have yet to be addressed, primarily because the enabling technologies are newly developed, incompletely understood, and poorly integrated.  My proposal aims to develop a basic science understanding of how multiple omic-scale molecular readouts of cellular physiology should be combined with genetic information to establish molecular "fingerprints" predictive of large-scale phenotypes, such as health, using yeast as a tractable model system. Specifically, I will study a multifaceted process in yeast that approximates wellness, aging. The project will be conducted within the context of a tightly controlled genetic and experimental system so that I can directly evaluate relationships predicted between genetics, longitudinal molecular fingerprints, the environment, and quantitative phenotypes using targeted genetic interventions. I will combine multiple high-throughput molecular technologies with computational and statistical methods to develop predictive models of cellular function.

A growing number of technological developments are changing the scale of phenotyping in human populations [2]. Apple, for instance, recently unveiled a ResearchKit that will allow people using its mobile technologies to enroll in human research studies. There is, therefore, a growing need to understand how to harness large collections of heterogeneous data and integrate them with genetic information to model and ultimately predict biological processes in individuals. This is a problem best addressed by a basic science approach. As an HFSP fellow, I will contribute to this emerging demand by developing data-driven approaches to answer two fundamental questions in systems genetics: (1) how does the influence of genetic variation change across molecular scales, from immediate effects on gene expression to downstream consequences for metabolic pathways and ultimately larger-scale phenotypes like aging? and (2) how are genetic influences and their molecular fingerprints shaped by the environment? I will conduct my experiments in a yeast population that is both genetically and experimentally tractable. The experimental methods will combine omic-scale measurements with a microfluidic system that can image multiple features of single cells and precisely control their environments. Together, these datasets will allow me to quantify genetic, molecular, and environmental factors that influence aging in a genetically diverse population using state-of-the-art machine learning methods.

# Experimental Design and Methodology

The design of this project integrates experimental approaches from several disciplines with computation. I have an interdisciplinary background with expertise in both experimental and computational approaches.

Microfluidics will be used to measure multiple aging phenotypes in mother cells. High-throughput omic technologies will be applied to quantify the molecular composition of cells at several stages of the aging process. Statistical learning algorithms will be applied to these large, heterogeneous datasets to derive multiomic fingerprints predictive of aging state. Finally, targeted genetic approaches will be used confirm whether genetic variation predicted by the model have a quantitative effect on fingerprints or aging. Each of these methods will be performed in a genetically diverse yeast population of 140 haploid segregants in two growth conditions to characterize the environmental (in)dependence of the genetic markers and their resulting age-specific molecular fingerprints.

## Yeast as a model system for aging

Yeast is a well established model of aging in eukaryotes. Longevity in yeast can be quantified as replicative life span (RLS), or the number of daughters produced before irreversible cell cycle arrest. Recently, a microfluidic platform was developed to measure RLS [3]. This system isolates mother cells in a tightly controlled environment where physical features of the cell (e.g., cell size and shape) can be quantified in addition to RLS. The precise environmental control exerted by this system will be used later in the experimental methodology to characterize how the environment modifies aging phenotypes. As a first step in the project, I will establish this experimental system at EMBL to measure RLS and related phenotypes in a genetically diverse yeast population described below. Several groups at EMBL have expertise in microfluidics. I will work closely with these groups to establish this technology in our lab.

## Identifying genetic determinants of aging in a segregating yeast population

A critical feature of this project is that it will be conducted within a model system that is genetically tractable. This means that phenotypic variations can be mapped to genetic loci without having to account for substantial confounding factors (e.g., population structure) like one would have to do in human populations. The Steinmetz lab maintains haploid yeast strains from a cross between the laboratory strain of S. cerevisiae, S288c, and a pathogenic strain isolated from the lungs of an AIDS patient, YJM789 [4], that has been used previously to identify genetic factors related to several biological processes, including gene expression [5]. The cross is high-resolution, with ~52,000 genetic markers spread across the genome at a median distance of 78 basepairs between markers. 140 of these haploid yeast segregants will be characterized for RLS by imaging using the microfluidic system describe above. Genotypes associated with variation in these phenotypes will be detected by quantitative trait loci (QTL) mapping using an ensemble learning method called random forests. Recent studies have shown that random forest-based QTL approaches have greater statistical power to detect causal variants and can even detect nonlinear interactions between loci [6].

Surprisingly, a recent study of RLS in yeast was able to identify only five significant genetic loci (QTLs) related to aging [7], which is unexpected given the expected complexity of aging. My study will greatly increase the power to detect causal genetic factors because it will include more sensitive computational methods (random forest-based QTL detection), a more finely quantified phenotype (microfluidics), and a population with far greater genomic resolution (~52,000 genetic markers compared to ~3,000).

## Multiomic molecular fingerprints of aging

In parallel, I will measure the molecular composition of each yeast strain at several stages of the aging process, from "young" to "old" cells. All of the methods described below are routinely applied in either the Steinmetz lab or in the labs of our collaborators. We have active collaborations at EMBL, Stanford, and the University of Luxembourg where I will perform each these measurements.

### Isolation of cells of a specific age

To approximate longitudinal profiling in a yeast population, I will isolate and profile cells that are a specific age. Age-specific enrichment can be accomplished by activating an antibiotic encoding plasmid that is retained and replicated only in mother cells. With this system it is possible to both enrich for mother cells by selection and quantify their age by sequencing (i.e., counting plasmid copy number). Since this system is activated once in a single generation, it can be used to ensure that all mothers cells in a population are the same starting age. This technology is currently unpublished and being actively developed within the Steinmetz lab.

Given that a typical mother cell has an RLS of 20-30 divisions, I will collect samples for multiomic profiling across three age classes: "young" (generations 1-5), "middle-age" (generations 5-10) and "old" (generations 10-20) cells.

For each of three age classes, 140 strains with different genetic backgrounds and two environmental conditions, I will measure multiple omic-scale profiles. With biological replication a total of 6,720 genome-scale molecular profiles will be collected using the following technologies:

### Epigenomics

Epigenetic and regulatory DNA state will be profiled by DNase-seq. DNase-seq provides a readout of chromatin accessibility and can even measure small DNA regions bound by regulatory factors, genome-wide. Since chromatin state and gene regulation are known to be highly sensitive to cell state, these measurements will be particularly informative for modeling how the environment conditions genetic factors [8].

### Transcriptomics

The transcriptome of each strain will be measured using two methods developed in the Steinmetz lab, 5'P-seq and 3'Tag-seq. 3'Tag-seq quantifies polyadenylation sites, whereas 5'P-seq determines transcription start sites [9,10]. Both methods not only provide robust estimates of transcript abundance levels, but they also measure heterogeneity in the transcriptome. Since variation in isoform preference has recently been linked to physiological effects [11], it will interesting to measure how control of isoform usage varies with age.

### Proteomics

Comprehensive proteomic data will be collected with mass spectrometry. I will take advantage of a newly described sample preparation method developed at EMBL called Single-Pot Solid-Phase-enhanced Sample Preparation (SP3) [12]. SP3 has been optimized previously for proteomic analysis in yeast. Spike-in standards will be included to normalize across all batches and strains.

### Metabolomics

Quantification of select intra- and extracellular metabolites (absolute abundance) will be performed by targeted LC-MS/MS. The Steinmetz lab has an ongoing international collaboration that has now quantified nearly all metabolites in central carbon metabolism of yeast for more than 500 samples [unpublished].

## Environmental dependence of genetic factors and their molecular fingerprints

The environment modifies the effects of genetic variation [5]. I will determine how the environment alters the determinants and hallmarks of aging by leveraging a documented heat tolerance phenotype of the YJM789 strain. I will re-characterize RLS QTLs and molecular fingerprints at elevated growth temperatures in each segregant to assign an environmental (in)dependence to each genotype and molecular fingerprint.

### Data integration and computational modeling

Numerous studies have demonstrated the utility of integrating diverse data types, particularly to understand the function of complex biological systems. Each of the multiomics datasets collected in this study will be processed using established platform-specific pipelines with error models specific for each data type. These heterogeneous datasets will be combined and integrated with prior knowledge, including eQTL maps from a previous study in the lab [5], using data fusion approaches [13].  Molecular fingerprints that are diagnostic of aging in the context of specific genetic backgrounds will be identified using ensemble learning methods, particularly random forests.

## Validation

The computational models will implicate specific genetic loci that influence aging. A final component of the project will be to validate several of these loci using targeted allele replacement with CRISPR/Cas-9 technology.  Briefly, specific genotypes will be swapped at the single-nucleotide level to test whether particular variants are causally associated with RLS and/or molecular fingerprints diagnostic of a particular age.

# Risk management and feasibility

This project leverages a variety of methods that are either well established or in development. To minimize the risk of such a large and interdisciplinary project, I will rely on an worldwide network of collaborators that has been established by my host supervisor. Within this network I will have support from leading experts for every component of this project.

# Summary and Significance

Aging is a complex, poorly understood biological process that is influenced by both genetic and non-genetic factors. I expect to observe differences in aging not only between strains with different genetic backgrounds but also within isogenic strains at different stages of the aging process. By combining molecular fingerprints with genetics, I will be able to identify single and combinations of loci that predict or directly influence the aging process. The model will make predictions about the molecular scale at which particular genetic variations manifest, from immediate consequences on transcription to downstream effects on metabolism, and will distinguish genetic factors with environmental dependence. Deconstruction of this complex phenotype will help determine whether molecular profiles and genetics are sufficient to diagnose a complex biological process and provide proof-of-principle approaches for integrating large heterogeneous datasets with genetic information to drive longitudinal, data-driven evaluation of complex phenotypes, including health in humans.  

[1] Hood L, Price N. Demystifying Disease, Democratizing Healthcare. (2014) Science Translational Medicine. 6:225
[2] The coming era of human phenotyping. (2015) Nat Biotechnol. 33(6):567.
[3] Mancera E, et al. High-resolution mapping of meiotic crossovers and non-crossovers in yeast. (2008) Nature. 454(7203):479-85.
[4] Jo MC, et al. High-throughput analysis of yeast replicative aging using a microfluidic system. (2015) PNAS. 112(3):9364-9369.
[5] Gagneur J, et al. Genotype-environment interactions reveal causal pathways that mediate genetic effects on phenotype. (2013) PLoS Genet. 9(9):e1003803
[6] Stephan J, et al. A random forest approach to capture genetic effects in the presence of population structure. (2015) Nat Commun. 6:7432.
[7] Stumpferl SW, et al. Natural genetic variation in yeast longevity. (2012) Genome Res. 22(10):1963-73.
[8] Hesselberth JR, et al. Global mapping of protein-DNA interactions in vivo by digital genomic footprinting. (2009) Nat Methods. 6(4):283-9
[9] Wilkening S, et al. An efficient method for genome-wide polyadenylation site mapping and RNA quantification. (2013) Nucleic Acids Res. 1;41(5):e65.
[10] Pelechano V, et al. Widespread Co-translational RNA Decay Reveals Ribosome Dynamics. (2015) Cell. 4;161(6):1400-12.
[11] Pelechano V, et al. Extensive transcriptional heterogeneity revealed by isoform profiling. (2013) Nature. 497(7447):127-31
[12] Hughes CS, et al. Ultrasensitive proteome analysis using paramagnetic bead technology. (2014) Mol Syst Biol. 10:757
[13] Wang B, et al. Similarity network fusion for aggregating data types on a genomic scale. Nat Methods. 11(3):333-7.
