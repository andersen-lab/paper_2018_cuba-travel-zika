## Data & analysis

**Grubaugh et al., 'International travelers and genomics uncover a ‘hidden’ Zika outbreak'**

---

### Local Zika cases

Monthly Zika cases were digitized from bar graphs published by the Pan American Health Organization (PAHO), an arm of the World Health Organization. See our blog post for a description of how we digitized the data: https://andersen-lab.com/paho-zika-cases/

* Data are located in '/local_Zika_cases'

### Travel data and Zika cases

The Florida Department of Health kindly curated their Weekly Arbovirus Reports from 2016-2018 (last updated in October 2018) to build a database of travel-associated Zika virus cases to monitor for ongoing transmission throughout the Americas. The data is sorted by symptoms onset date (month) and origin of travel. Also included are the air travel volumes from the sources of travel Zika cases into Florida to calculate travel incidence, obtained from the International Air Transport Association.

* Data are located in '/travel_data'

### Phylogenetics

A total of 283 Zika virus genomes collected between 2013 and 2018 from Cuba (n = 10, including 9 generated in this study) and elsewhere from the Pacific and the Americas (n = 273, including 4 generated in this study from Florida, USA) were codon-aligned together using MAFFT and inspected manually. A maximum likelihood phylogeny was reconstructed with RAxML using the general time-reversible (GTR) nucleotide substitution model and gamma-distributed rates amongst sites. A time-scaled phylogenetic was reconstructed using the Bayesian phylogenetic inference framework available in BEAST v1.10.2. Accommodating phylogenetic uncertainty, we used an HKY+Γ4 nucleotide substitution model for each codon position, allowing for relative rates between these positions to be estimated, and an uncorrelated relaxed molecular clock model, with an underlying lognormal distribution, a non-parametric demographic prior and otherwise default priors in BEAUti v1.10.2. The MCMC analysis was run for 800 million iterations, sampling every 100,000th iteration, using the BEAGLE library v2.1.2 to accelerate computation. MCMC performance was inspected for convergence and for sufficient sampling using Tracer v.1.7.1. After discarding the first 200 million iterations as burn-in, virus diffusion over time and space was summarized using a maximum clade credibility (MCC) tree using TreeAnnotator.

* Data and trees are located in '/phylogenetics'

### Estimating Zika cases in Cuba

* Description of data, code, and results are located in '/estimating_cases_in_Cuba'

### Modeling Aedes aegypti-borne virus transmission

* Description of data, code, and results are located in '/mosquito_transmission_model'

---
**Andersen Lab**  
The Scripps Research Institute  
La Jolla, CA, USA  
[data@andersen-lab.com](mailto:data@andersen-lab.com)
