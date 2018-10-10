var detailsArray = ['20180515_SuSiE_Benchmark']
var detailsDict = {'Setup-of-SuSiE-benchmark-1': '20180515_SuSiE_Benchmark'}
var detailsArrayMap = {'20180515_SuSiE_Benchmark': 'Setup of SuSiE benchmark'}
var analysisArray = ['20180918_Use_Penalty', '20180912_Case_Examples_Li2016', '20180911_BF_Exporation_MS', '20180902_Case_Examples_Li2016', '20180829_LD_Heatmap', '20180712_Enrichment_Workflow', '20180711_A_Hard_Case', '20180709_SplicingQTL_Detailed', '20180704_MolecularQTL_Workflow', '20180630_Runtime', '20180620_Purity_Plot_Lite', '20180615_Power_DAP', '20180606_ROC', '20180606_Identify_Interesting_Dataset', '20180606_Coverage_Check', '20180605_PIP_Calibrated', '20180601_One_Causal_CAVIAR_Exam', '20180531_PIP_L1_Comparison', '20180528_Power_Comparison', '20180527_PIP_Workflow', '20180527_PIP_Comparison', '20180525_CS_Outliers', '20180522_Purity_Summary', 'purity', '20180516_Purity_Plot', '20180515_SusieR_Benchmark', '20180515_Extract_Benchmark_Data', '20180508_ELBO', '20180424_DSC_FMO2', '20180416_SingleCondition_FMO2', '20180415_MNMASH_FMO2', '20171210_Simulation_Multivariate', '20171130_MNM_Toys', '20171103_MNMASH_Data', '20170630_Simulation_Study', '20170624_Simulation_Procedures']
var analysisDict = {'Simulation-of-quantitative-phenotype-given-genotypes-1': '20170624_Simulation_Procedures', 'mr-ash-simulation-ash-paper-scenarios-1': '20170630_Simulation_Study', 'Prepare-toy-example-data-set-for-M&M-ASH-model-1': '20171103_MNMASH_Data', 'Toy-M&M-analysis-on-Thyroid-and-Lung-1': '20171130_MNM_Toys', 'Simulation-of-multiple-phenotypes-given-genotypes-and-covariance-1': '20171210_Simulation_Multivariate', 'M&M-analysis-on-FMO2-data-in-GTEx-1': '20180415_MNMASH_FMO2', 'Univariate-analysis-on-FMO2-data-in-GTEx-1': '20180416_SingleCondition_FMO2', 'Fine-mapping-on-FMO2-data-in-GTEx-1': '20180424_DSC_FMO2', 'M&M-ELBO-1': '20180508_ELBO', 'Extract-per-gene-toy-dataset-1': '20180515_Extract_Benchmark_Data', 'SusieR-benchmark-1': '20180515_SusieR_Benchmark', 'SusieR-benchmark-plot-1': '20180516_Purity_Plot', 'Purity-result-summary-1': '20180522_Purity_Summary', 'CS-outlier-scenarios-1': '20180525_CS_Outliers', 'Direct-comparison-of-PIP-for-SuSiE,-DAP,-CAVIAR-and-FINEMAP-1': '20180527_PIP_Comparison', 'Workflow-to-extract-PIP-and-set-information-for-different-methods-1': '20180527_PIP_Workflow', 'Power-comparison-susie-vs-DAP-1': '20180528_Power_Comparison', 'Compare-PIP-of-L1-susie,-DAP-and-CAVIAR-1': '20180531_PIP_L1_Comparison', 'Looking-into-the-3-CAVIAR-outliers-1': '20180601_One_Causal_CAVIAR_Exam', 'Calibration-of-SNP-level-PIP-1': '20180605_PIP_Calibrated', 'Check-susie-coverage-1': '20180606_Coverage_Check', 'Identify-and-extract-interesting-data-set-for-vignettes-1': '20180606_Identify_Interesting_Dataset', 'ROC-comparisons-1': '20180606_ROC', 'Workflow-to-extract-info-for-power-comparison-with-DAP-for-a-hard-case-1': '20180615_Power_DAP', 'SuSiE-Purity-Plot-1': '20180620_Purity_Plot_Lite', 'Comparing-computational-efficiency-of-methods-1': '20180630_Runtime', 'Molecular-QTL-workflow-1': '20180704_MolecularQTL_Workflow', 'A-detailed-look-at-some-of-the-SuSiE-fits-1': '20180709_SplicingQTL_Detailed', 'A-hard-case-fine-mapping-example-1': '20180711_A_Hard_Case', 'Enrichment-analysis-workflow-for-molecular-QTL-results-1': '20180712_Enrichment_Workflow', 'Make-LD-heatmap-for-demonstration-data-sets-1': '20180829_LD_Heatmap', 'Investigating-sQTL-analysis-results-of-interest:-an-overview-1': '20180902_Case_Examples_Li2016', 'knitr::opts-chunk-set(echo-TRUE)-1': '20180911_BF_Exporation_MS', 'Investigating-sQTL-analysis-results-of-interest:-a-deeper-look-1': '20180912_Case_Examples_Li2016', 'Adding-null-component-to-SuSiE-1': '20180918_Use_Penalty'}
var analysisArrayMap = {'20170624_Simulation_Procedures': 'Simulation of quantitat ... en genotypes', '20170630_Simulation_Study': 'mr-ash simulation ash paper scenarios', '20171103_MNMASH_Data': 'Prepare toy example dat ... &M ASH model', '20171130_MNM_Toys': 'Toy M&M analysis on Thyroid and Lung', '20171210_Simulation_Multivariate': 'Simulation of multiple  ... d covariance', '20180415_MNMASH_FMO2': 'M&M analysis on FMO2 data in GTEx', '20180416_SingleCondition_FMO2': 'Univariate analysis on  ... data in GTEx', '20180424_DSC_FMO2': 'Fine mapping on FMO2 data in GTEx', '20180508_ELBO': 'M&M ELBO', '20180515_Extract_Benchmark_Data': 'Extract per gene toy dataset', '20180515_SusieR_Benchmark': 'SusieR benchmark', '20180516_Purity_Plot': 'SusieR benchmark plot', '20180522_Purity_Summary': 'Purity result summary', '20180525_CS_Outliers': 'CS outlier scenarios', '20180527_PIP_Comparison': 'Direct comparison of PI ... and FINEMAP', '20180527_PIP_Workflow': 'Workflow to extract PIP ... rent methods', '20180528_Power_Comparison': 'Power comparison susie vs DAP', '20180531_PIP_L1_Comparison': 'Compare PIP of L1 susie ... P and CAVIAR', '20180601_One_Causal_CAVIAR_Exam': 'Looking into the 3 CAVIAR outliers', '20180605_PIP_Calibrated': 'Calibration of SNP level PIP', '20180606_Coverage_Check': 'Check susie coverage', '20180606_Identify_Interesting_Dataset': 'Identify and extract in ... or vignettes', '20180606_ROC': 'ROC comparisons', '20180615_Power_DAP': 'Workflow to extract inf ... a hard case', '20180620_Purity_Plot_Lite': 'SuSiE Purity Plot', '20180630_Runtime': 'Comparing computational ... y of methods', '20180704_MolecularQTL_Workflow': 'Molecular QTL workflow', '20180709_SplicingQTL_Detailed': 'A detailed look at some ... e SuSiE fits', '20180711_A_Hard_Case': 'A hard case fine-mapping example', '20180712_Enrichment_Workflow': 'Enrichment analysis wor ... QTL results', '20180829_LD_Heatmap': 'Make LD heatmap for dem ... on data-sets', '20180902_Case_Examples_Li2016': 'Investigating sQTL anal ... an overview', '20180911_BF_Exporation_MS': 'knitr::opts_chunk$set(echo = TRUE)', '20180912_Case_Examples_Li2016': 'Investigating sQTL anal ... deeper look', '20180918_Use_Penalty': 'Adding null component to SuSiE'}
var prototypeArray = ['20180409_Model_DSC', '20171207_Tryout_Parallel', '20171207_MNMASH_Model', '20171203_Vectorized_OLS', '20171129_MNMASH_Model', '20171127_MNMASH_Model', '20171103_MNMASH_Model', '20170628_MR_ASH_Toy_Example', '20170615_MASHR_Benchmark']
var prototypeDict = {'mashr-R-vs.-C++-benchmark-1': '20170615_MASHR_Benchmark', 'mr-ash-example-analysis-1': '20170628_MR_ASH_Toy_Example', 'Prototype-of-VEM-in-M&M-ASH-model-1': '20171103_MNMASH_Model', 'Prototype-of-core-update-in-M&M-ASH-model-1': '20171127_MNMASH_Model', 'M&M-model-VEM-updates-1': '20171207_MNMASH_Model', 'Computing-vectorized-OLS-1': '20171203_Vectorized_OLS', 'Try-out-paralleled-numpy-computations-1': '20171207_Tryout_Parallel', 'Breaking-M&M-prototyping-to-using-DSC-1': '20180409_Model_DSC'}
var prototypeArrayMap = {'20170615_MASHR_Benchmark': 'mashr R vs. C++ benchmark', '20170628_MR_ASH_Toy_Example': 'mr-ash example analysis', '20171103_MNMASH_Model': 'Prototype of VEM in M&M ASH model', '20171127_MNMASH_Model': 'Prototype of core updat ... &M ASH model', '20171129_MNMASH_Model': 'M&M model VEM updates', '20171203_Vectorized_OLS': 'Computing vectorized OLS', '20171207_MNMASH_Model': 'M&M model VEM updates', '20171207_Tryout_Parallel': 'Try out paralleled numpy computations', '20180409_Model_DSC': 'Breaking M&M prototyping to using DSC'}
var dscArray = ['Untitled', 'Untitled1', 'finemapping', 'export.pipeline']
var dscDict = {'dat-readRDS(-susie-comparison-liter-data-liter-data-21.ld-mat.rds-)-1': 'Untitled', 'library(dplyr)-1': 'Untitled1', 'Finemapping-benchmark-1': 'finemapping'}
var dscArrayMap = {'Untitled': "dat = readRDS('susie_co ... ld_mat.rds')", 'Untitled1': 'library(dplyr)', 'finemapping': 'Finemapping benchmark'}
var writeupArray = ['To_Explore', 'Status_MM_ASH', 'Start_Simple', 'Modular_MNMASH', 'Meetings', '20171215_MNMModel_Finemap', '20171203_VEM']
var writeupDict = {'Notes-on-fundamentals-of-EM-and-VEM-1': '20171203_VEM', 'M&M-model-for-fine-mapping-1': '20171215_MNMModel_Finemap', 'Meetings-1': 'Meetings', 'A-modular-approach-to-M&M-ASH-model-1': 'Modular_MNMASH', 'Start-Simple!-1': 'Start_Simple', 'Status-of-m&m-ash-project-1': 'Status_MM_ASH', 'To-explore-1': 'To_Explore'}
var writeupArrayMap = {'20171203_VEM': 'Notes on fundamentals of EM and VEM', '20171215_MNMModel_Finemap': 'M&M model for fine-mapping', 'Meetings': 'Meetings', 'Modular_MNMASH': 'A modular approach to M&M ASH model', 'Start_Simple': 'Start Simple!', 'Status_MM_ASH': 'Status of m&m ash project', 'To_Explore': 'To explore'}