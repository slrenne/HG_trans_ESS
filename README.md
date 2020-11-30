# High grade translocated endometrial stromal sarcoma
The project repository for the HG_trans_ESS study. 

## _Aim_
We aim to understand the role of specific translocations in rearranged high grade endometrial stromal sarcoma.
To do this we will perform a Bayesian analysis of the published series and of our cases (partly presented at the European Congress of Pathology).<sup>[1](#myfootnote1)</sup>
## _Structure_ 
### Retrospective series

_Inclusion criteria_
* Histologically confirmed diagnosis of –BCOR rearranged Endometrial Stromal Sarcoma
* Localized primary, recurrent and/or metastatic disease 

_Exclusion criteria_
* Absence of BCOR rearrangement
* Diagnosis different from Uterine Sarcoma

_Experimental design_

All endometrial stromal sarcomas with -BCOR rearrangement diagnosed at Fondazione IRCCS Istituto Nazionale dei Tumori, Milan, Italy will be included in the study. All available stained slides will be retrieved and reviewed; the presence of the following morphological parameters will be evaluated on the earliest available specimen: tongue like projections, lympho-vascular space involvement, presence of collagen plaques, spindle fascicular appeareance, myxoid matrix, cytologic atypia; moreover, mitotic figures will be counted on 10 high-power field (HPF; an HPF measuring 0.1734 mm²) on the earliest available material. It will also be noted whenever primary untreated tumor will not be the one evaluated. Undescribed morphological characteristics will be also recorded. Eight edition TNM parameters and AJCC Staging Group, will be performed; local and systemic therapies will be gathered from medical records. The immunohistochemical results obtained during original diagnosis will be reviewed. Appropriate control will be used whenever internal control is not available. EnVsion Flex in dako autostainer will be used. Staining will be evaluated as positive, focal (<20%), or negative.  When intensity will be heterogeneous the strongest intensity will be recorded. Fluorescence In Situ Hybridization (FISH) analysis results obtained during original diagnosis will be recorded. Gene rearrangement will be evaluated in formalin-fixed, paraffin-embedded, 4-μm-thick sections. At least 60 non-overlapping nuclei will be analyzed, cases in which 80% of the cells showed abnormal signal will be considered positive.

_Institutional Review Board_

The study was reviewed and approved by the institutional review board with the number INT 19/28.

### Literature Research
We plan to perform a systematic search of the literature of the public databases - PubMed, Word of Science, and Global index medicus - using the following terms:
_((endometrial stromal sarcoma) OR (stromal sarcoma uterus)) AND ((cytogenetics) OR (fusion) OR (rearrangement) OR (bcor) OR (ywhae)) AND ((survival) OR (prognosis) OR (follow-up))_

### Database
Research products will be stored in data/published_data/.
Raw data files extracted by each publication will be stored in data/raw_data/.
The merged dataset with our cases with those selected from the literature will be stored in data/merged_data. 

### Analysis
**Priors**
The prior distributions for the variables will be chosen following the maximum entropy distribution: variables with real values and finite variance will be modeled with a Gaussian distribution, binary events with fixed probability will be modeled with Binomial distribution, multivariate events Dirichlet distribution, non-negative real variables with a mean will be modeled with exponential distribution, real value within an interval with a uniform distribution, low rate events will be modeled with gamma distributions. Further distribution will be motivated if necessary. Prior Regularization will be performed with prior predictive simulation. 
Distributions and plots of prior predictive simulations will be stored in analysis/priors/.

**Models**
Directed acyclic graphs (DAGs) will be designed. Each variable will be modeled in the whole dataset and in the different rearrangement categories. Relevance of common prognosticators will be considered. GLM on survival status and survival analysis exponential models will be performed. 
DAGs and models will be stored in analysis/priors/.

## _Reference_ 
 <a name="myfootnote1">1</a>: Renne SL, Collini P, Dagrada GP, Gronchi A, Radaelli S, Sanfilippo R, Lorusso D, Raspagliesi F, Paolini B. BCOR rearranged endometrial stromal sarcoma (ESS). A clinicopathological study of four cases. Virchows Archiv, 473, 2018 p. S35, ISSN: 0945-6317

