# Multiomic fingerprints of wellness in a genetically tractable yeast population

A new paradigm for personalized medicine is emerging in the biomedical sciences: quantified wellness. Rather than focus on treatment of disease by dissecting mechanisms and etiologies, this approach aims to monitor and quantify patient "wellness" in "real-time", suggesting personalized interventions to sustain health before a patient manifests disease [1]. This exciting new approach to medicine, however, raises a number of basic science questions that have yet to be addressed, primarily because the technologies that would enable such large-scale longitudinal monitoring of a biological system are newly developed, incompletely understood, and poorly integrated.  My proposal aims to develop a basic science understanding of how multiple omic-scale molecular readouts of cellular physiology should be combined and integrated with genetic information to establish molecular "fingerprints" predictive of large-scale phenotypes, such as health, using yeast as a tractable model system. Specifically, I will study a multifaceted process in yeast that approximates wellness, aging. The project will be conducted within the context of a tightly controlled genetic and experimental system so that I can directly evaluate predicted relationships predicted between genetics, longitudinal molecular fingerprints, the environment, and quantitative phenotypes using targeted genetic, molecular, and environmental interventions. I will combine multiple high-troughput molecular technologies with state-of-the-art computational and statistical methods to develop integrated and predictive models of cellular function.

The project is motivated by simultaneous changes in biology and technology that are reshaping the practice of medicine. For decades, biomedical research has focused on understanding and treating disease mechanisms and etiologies. This has been a particularly powerful way to develop drugs that directly interface with disease mechanisms. More recently, next-generation sequencing and related technologies have allowed investigators to achieve deeper insight into the genetics determinants of disease. Genome Wide Association Studies (GWAS), for example, have identified multiple genetic loci associated with several human diseases, ranging from diseases with relatively simple pathophysiologies to much more complex scenarios. Despite many successes of these approaches, there are also notable shortcomings. In particular: (1) many diseases are polygenic with complex etiologies such that it is difficult to discern causal mechanism(s), (2) genetic loci identified by GWAS for complex diseases tend to have small effect sizes and oftentimes their precise relationship to disease is unclear, and (3) many common diseases may be driven by rare variations such that is difficult to identify the mechanism(s) responsible for disease in an individual. Unless addressed, these pitfalls will be magnified as the number and diversity of biological measurements increases. Apple, for example, has recently unveiled its ResearchKit, which will allow people using its mobile technologies to enroll in human research studies. This is just one of a growing number of technological developments that will dramatically change the scale of  phenotyping in human populations [2]. There is, therefore, a growing need to understand how large collections of heterogeneous data should be integrated with genetic information to understand, model, and ultimately predict fundamental biological processes in individuals.

This is a question best addressed by basic science approaches. There are a number of fundamental scientific questions embedded within this problem that cut across disciplines, from mathematics to computer science and biochemistry. As an HFSP fellow, I will attempt to relate several these questions to my field of systems biology, namely: (1) how the influence of genetics varies across molecular scales and (2) how genetic influences and their molecular fingerprints are shaped by the environment. I will conduct my experiments in a yeast population that is both genetically and experimentally tractable. My experiments will combine recent developments in omic-scale measurements with microfluidic systems that can image multiple features of single cells and precisely control their environments. Together, these data sets will allow me to quantify genetic, molecular, and environmental factors that influence aging.

# Experimental Design and Methodology

The experimental design for this project integrates experimental methods from several disciplines with computation. I have an interdisciplinary background with expertise in both experimental and computational approaches. Microfluifdics will be used to measure aging phenotypes in mother cells. Multiple high-throughput methods will be applied to quantify molecular fingerprints. Statistical learning approaches will be applied to these large, heterogeneous data sets to derive multiomic fingerprints indicative of the aging process. Finally, targeted genetic approaches will be used confirm novel associations. All of these will be conducted in a genetically diverse yeast popuation across two growth conditions to elucidate and their environmental (in)dependence.

## Yeast as a model system for aging

Yeast is a well established model of aging in eukaryotes. A yeast’s longevity can be quantified as replicative life span (RLS), or the number of daughters produced before irreversible cell cycle arrest. Recently, a microfluidic platform has been developed to measure RLS [1]. This system can isolate mother cell in tightly controlled environments. Importantly, using this system the mother cell can be profiled for, including cell shape and size, in addition to RLS. The ability to precisely control the environment of cells we be used again below to determine how the environment conditions aging responses. As a first step in my project I will establish an experimental system to measure RLS and related parameters in a genetically diverse yeast population described below. Several groups at EMBL have expertise in microfluidics. I will work closely with these groups to establish this technology in our lab.

## Identifying genetic determinants of aging in segregrating yeast population

Each of 140 haploid yeast segregants from a cross between the laboratory strain of S. cerevisiae, S288c, and a pathogenic strain isolated from the lungs of an AIDS patient, YJM789, will be characterized for RLS.

A recent study was able to identify only five significant genetic loci (QTLs) related to RLS in yeast [2]. Our study will significantly increase the power to detect additional genetic factors by combining more sensitive computational methods (random forest-based QTL detection [3]) with more finely quantified phenotypes and segregants with greater genomic resolution (>15,000 markers compared to ~3,000).

## Multiomic molecular fingerprints of aging

In parallel, I will. Many of these methods are currently being optimized. We have active collaborations at EMBL, Stanford, and the University of Luxembourg to help perform these studies.

## Epigenomics

Epigenetic state and regulatory DNA will be profiled by DNase-seq

## Transcriptomics

The transcriptomes of using two methods developed in the Steinmetz lab, 5'PSeq and 3' Tagseq. 3' Tagseq quantifies polyadenylation sites, while 5'Pseq determines the transcription start sites []. These methods provide robust estimations of transcript abundance levels and also reveal the complexity of the transcriptome. The additional information on transcript variation will help us to pinpoint mitochondrial phenotypes caused by transcript isoform variation.

## Proteomeomics

Comprehensive proteomic data will be collected using data-dependent acquisition on a QExactive Orbitrap mass spectrometer. I will take advantage of the SP3 method optimized for proteomic analysis in yeast, to prepare cellular extracts for MS analysis []. A spike-in standard will be included to normalize across all batches and strains. Subsequent steps including enzymatic digestion, LC-MS/MS measurement and computational analysis will be carried out using well-established protocols available in the Steinmetz lab.

## Metabolomics

Quantification of select intra- and extracellular metabolites (absolute abundance). In collaboration we have been optimizing a multi-step procedure to reproducibly generate biomass, withdraw defined amounts of biological material, measure absolute biovolumes (i.e., total amount of cells surrounded by an intact membrane) by volume displacement (coulter counter principle) [], precisely inactivate metabolism by cold methanol quenching , nonselectively extract and measure the metabolites quantitatively. Targeted analysis of absolute intra- and extracellular metabolite concentrations will be performed using an LC-MS/MS. Currently, using a QExactive Orbitrap LS-MS/MS we are able to quantify nearly all metabolites in central carbon metabolism of yeast.

## Data integration and computational modeling

Each of the multiomics data sets will be processed using established pipelines. Rigorous quality control will be enforced and the resulting data integrated into a single multi-omic dataset using error models specific for each data type. Novel associations between specific genotypes and molecular phenotypes within this dataset, we will employ statistical learning techniques. Specifically, we will use ensemble learning methods, such as random forests, to discover meaningful associations across diverse molecular features []. These predicted associations will be interpreted mechanistically within a systems-level model of metabolism and genetic regulation. Thus, we will be able to predict specific mechanisms through which several characterized genotypes lead to an observed phenotype. We will build this systems-level model by integrating up-to-date, public, genome-scale models of metabolism and genetic regulation in yeast using established methods.

These heterogeneous datasets will be combined and integrated with prior knowledge, including eQTL maps from a previous study in the lab [4], using a machine learning approach that I will build to identify combinations of molecular profiles that are diagnostic of the aging process in the context of particular genetic backgrounds.

## Environmental dependence

I will determine how the environment conditions genetic contributions to aging by leveraging a documented heat tolerance phenotype of YJM789. I will re-characterize QTLs and molecular fingerprints at elevated growth temperatures to assign an environmental (in)dependence to each of the genotypes and molecular reporters identified above.

## Validation

# Risk management and feasibility

# Significance

This project generated will open source framework that can be used to study, including human wellness. The computational component will consist of statistical learning approaches that can be generalized to other types of data, such as those routinely collected from human patients.

More fundamentally, this study will begin to address  

[1] Hood L, Price N. Demystifying Disease, Democratizing Healthcare. (2014) Science Translational Medicine. 6:225
[2] The coming era of human phenotyping. (2015) Nat Biotechnol. 33(6):567.
[3] Jo MC, Liu W, Gu L, Dang W, Qin L. High-throughput analysis of yeast replicative aging using a microfluidic system. (2015) PNAS. 112(3):9364-9369.
[4] Stumpferl SW, Brand SE, Jiang JC, Korona B, Tiwari A, Dai J, Seo JG, Jazwinski SM. Natural genetic variation in yeast longevity. (2012) Genome Res. 22(10):1963-73.
[5] Stephan J, Stegle O, Beyer A. A random forest approach to capture genetic effects in the presence of population structure.
(2015) Nat Commun. 6:7432.
[6] Gagneur J, Stegle O, Zhu C, Jakob P, Tekkedil MM, Aiyar RS, Schuon AK, Pe'er D, Steinmetz LM. Genotype-environment interactions reveal causal pathways that mediate genetic effects on phenotype. (2013) PLoS Genet. 9(9):e1003803
[] Hughes CS, Foehr S, Garfield DA, Furlong EE, Steinmetz LM, Krijgsveld J. Ultrasensitive proteome analysis using paramagnetic bead technology. (2014) Mol Syst Biol. 10:757
[] Sonnleitner B, Locher G, Fiechter A. Biomass determination. (1992) J Biotechnol. 25(1-2):5-22.
[] Wilkening S, Pelechano V, Järvelin AI, Tekkedil MM, Anders S, Benes V, Steinmetz LM. An efficient method for genome-wide polyadenylation site mapping and RNA quantification. (2013) Nucleic Acids Res. 1;41(5):e65.
[] Pelechano V, Wei W, Steinmetz LM. Widespread Co-translational RNA Decay Reveals Ribosome Dynamics. (2015) Cell. 4;161(6):1400-12.
.
