Research plan, including references
(maximum of 15,000 characters including spaces and punctuation)

## Multiomic fingerprints of wellness in a genetically tractable yeast population

A new paradigm for personalized medicine is emerging in the biomedical sciences: quantified wellness. Rather than focus on treatment of disease by understanding disease mechanisms and etiologies, this approach to medicine aims to monitor and quantify the "wellness" of patients in real time, providing customized interventions to sustain health before a patient manifests disease. integration of and social technologies to inform . Many of the technologies that enable this approach, however, are newly developed, incompletely understood, and not well integrated.  My proposal aims to develop experimental methods to measure multiple omic-scale molecular readouts of cellular physiology and computational methods to integrate these diverse data types with genetic information to predict molecular "fingerprints" related to health using yeast as a genetically tractable model system. The study will be conducted within the context of a tightly controlled genetic and experimental model system, enabling direct evaluation of the relationships predicted between genetics, longitudinal molecular fingerprints, and quantitative measures of health using targeted genetic, molecular, and environmental interventions.

My project will have two immediate implications. First, the computational methods generated will constitute a proof-of-principle framework that can be used in the future for understanding human wellness. The computational component will consist of statistical learning methods that can be generalized to identify predictive fingerprints from the types of data routinely collected from human patients. Second, I will use data from this study to address two fundamental biological questions that are of particular interest to me. Specifically, I will address (1) how the influence of genetic variation varies across biological scales and (2) how genetic factors and their molecular consequences are shaped by the environment.

The project itself is motivated by simultaneous changes in biomedicine and technology. For decades, biomedical research has focused on understanding and treating disease mechanisms and etiologies. This has been a particularly powerful way to develop drugs that directly interface with disease mechanisms []. More recently, the advent of next-generation sequencing has allowed identification of genetic determinants for more complex diseases. Genome Wide Association Studies (GWAS), for example, have identified genetic factors associated with disease in human populations, including X and Y []. Despite many successes with these approaches, there are also notable shortcomings. In particular: (1) many diseases are polygenic with complex etiologies such that it is difficult to discern the mechanism(s) to drug, (2) genetic loci identified by GWAS for complex diseases tend to have small effect sizes and oftentimes their precise relationship to disease is unclear such that it is difficult to design treatments, and (3) many common diseases may be driven by rare variations such that is difficult to identify the mechanism(s) responsible for disease in a particular individual.

Personalized or precision medicine are first efforts to use modern technologies to tailor treatments to individuals. A recent variant of these approaches focuses on quantification of "wellness" rather than treatment of disease. This approach asks whether it can be more efficacious to monitor and maintain wellness rather than treat common diseases. Such an approach relies on an ability to monitor markers of health and disease in "real time". Such longitudinal profiling allows each individual to serve as their own control (allowing for personalization). Nearly "real time" monitoring of wellness has been facilitated by technological developments []. Indeed, several companies now offer services to quantify wellness and provide data-driven approaches to health management. These companies typically provide DNA sequencing, longitudinal molecular profiling of blood and other biomarkers, activity tracking through, as well as coaching. to improve health and detect anomalies before they manifest disease. From a technical point of view, however, it is unclear how genetic information can be combined with molecular profiles to predict or monitor disease development.

There are a number of fundamental scientific questions. In particular, there is currentl

Simple gentics. For these reasons, I to a model system, where can be developed within a tractabble. The project will have  As a HFSP fellow I will address these questions by quantifying the relationship between genetic, molecular, and environmental factors in a genetically tractable yeast population for a phenotype related to wellness, aging.

#### Experimental Design and Methodology

The experimental design for this project integrates as well as experimental biology and computation. I have an interdisciplinary background with expertise in both experimental and computational methods. Multiple high-throughput approaches will to quantify. Statistical learning approaches

##### Yeast as a model system for aging

Yeast is a well established model of aging in eukaryotes. A yeast’s longevity can be quantified as replicative life span (RLS), or the number of daughters produced before irreversible cell cycle arrest. High-throughput platforms have been established to measure RLS by isolating a single mother cell using microscopy and microfluidic devices [1]. As a first step in my project, I will establish an experimental system to measure RLS in a genetically diverse yeast population generated by the Steinmetz lab.  

##### Genetic determinants of aging in segregrating yeast population

Each of 140 haploid yeast segregants from a cross between the laboratory strain of S. cerevisiae, S288c, and a pathogenic strain isolated from the lungs of an AIDS patient, YJM789, will be characterized for RLS. A recent study identified SIR2 as a primary regulator of lifespan in yeast [2]. I will improve these measurements by using more sensitive QTL detection methods, a more finely measured phenotype, and segregants with greater genomic resolution to increase QTL detection sensitivity.

##### Multiomic fingerprints of wellness

###### Epigenomics

###### Transcriptomics

To comprehensively profile yeast transcriptomes, we will leverage the well-established 5'PSeq and 3' Tagseq methods and corresponding computational analysis developed in the Steinmetz lab. The 3' Tagseq method allows the quantification and identification of polyadenylation sites, while 5'Pseq accurately determines the transcription start sites (27, 56). These methods provide robust estimations of transcript abundance levels and also reveal the complexity of the transcriptome. The additional information on transcript variation will help us to pinpoint mitochondrial phenotypes caused by transcript isoform variation.

###### Proteomeomics

To comprehensively investigate proteomic data, we will perform data-dependent analysis using a QExactive Orbitrap mass spectrometer. We will take advantage of the modified SP3 method optimized for proteomic analysis in yeast, to prepare cellular extracts for MS analysis (28). A spike-in standard will be included to normalize across all batches and strains. Subsequent steps including enzymatic digestion, LC-MS/MS measurement and computational analysis will be carried out using well-established protocols from the GenPhen project.

###### Metabolomics

The generation of intra- and extracellular metabolic data from yeast (absolute abundance of metabolites) is a multi-step process: we will reproducibly generate biomass, withdraw defined amounts of biological material, measure the absolute biovolume (total amount of cells surrounded by an intact membrane) (57) by volume displacement (coulter counter principle) (58)., precisely inactivate metabolism by cold methanol quenching (59), nonselectively extract (60) and measure the metabolites (quantitatively and qualitatively). Within the proposed project a general protocol developed for the analysis of microorganisms andadapted to yeast within GenPhen will be used. Untargeted analysis of metabolites in different sample types will be carried out via an established GC-MS method (34). Targeted analysis of absolute intra- and extracellular metabolite concentrations will be performed using the LCMS/MS platform established within GenPhen. The measurement methods available on the QExactive Orbitrap LS-MS/MS cover almost all metabolites of the central carbon metabolism of yeast and provide absolute concentration via the use of isotopic labelled mass spectrometry coupled to chromatography.

###### Data integration and computational modeling

We will process the multi-omics data using established pipelines. After rigorous quality control, the data will be integrated into a single multi-omic dataset. To identify novel associations between specific genotypes and molecular phenotypes within this dataset, we will employ statistical learning techniques. Specifically, we will use ensemble learning methods, such as random forests, to discover meaningful associations across diverse molecular features (25). These predicted associations will be interpreted mechanistically within a systems-level model of metabolism and genetic regulation. Thus, we will be able to predict specific mechanisms through which several characterized genotypes lead to an observed phenotype. We will build this systems-level model by integrating up-to-date, public, genome-scale models of metabolism and genetic regulation in yeast using established methods (61).

In parallel, I will measure molecular "fingerprints" in each segregant at several stages of aging using a technique pioneered by another EIPOD in the Steinmetz lab to isolate individuals of a specified age. These fingerprints will include measurements of the epigenome and regulatory DNA by DNase-seq, the transcriptome by RNA-sequencing, the proteome by SP3 MS developed at EMBL, and the metabolome by targeted LC-MS to quantify metabolites of central metabolism. These heterogeneous datasets will be combined and integrated with prior knowledge using a machine learning approach that I will build to identify combinations of molecular profiles that are diagnostic of the aging process in the context of particular genetic backgrounds. Finally, I will determine how the environment conditions genetic contributions to aging by leveraging a documented heat tolerance phenotype of YJM789. I will re-characterize QTLs and molecular fingerprints at elevated growth temperatures to assign an environmental (in)dependence to each of the genotypes and molecular reporters identified above.

Having established quantitative relationships between genetics, multi-omic molecular profiles, and the environment, it would be enticing to return full circle to provide actionable insights to human wellness companies. As a final component of my project, therefore, I would like to partner with a company called Arivale (a company started by the head of my previous Institute, Lee Hood) to adapt my computational approaches to discover molecular fingerprints of human wellness.

##### Environmental dependence

##### Validation

#### Risk management and feasibility

#### Significance
