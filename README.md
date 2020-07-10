# Emergent-neutrality

This GutHub repository lists all R scripts used to produce the data analysed in D'Andrea et al. 2020 "Emergent neutrality in consumer-resource dynamics", acceptd for publication in PLoS Computational Biology.

*SimulationScript.R*: numerically simulates the stochastic niche dynamics model described in the paper. Directory paths must be changed. 

*AnalysisGeneralists.R*: reads data generated from *SimulationScript.R* and performs the species abundance distributions fits discussed in the paper.

*AnalysisSpecialists.R*: analogue of *AnalysisGeneralists.R*, except it reads data from the specialists scenarios (from a different data folder, available upon request).

*AnalysisExtinctionTimes.R*: reads data generated from *SimulationScript.R* and performs the analysis of extinction times discussed in the paper.

The *Data* folder includes partial data used in the analysis, namely consumer abundances in the Generalist scenario (please see file *Metadata.docx* for description of the data). 

_Note_: For access to the entire data set generated for this paper, please contact me at rafael.dandrea@stonybrook.edu. All data used in th paper can be generated with SimulationScript.R.

