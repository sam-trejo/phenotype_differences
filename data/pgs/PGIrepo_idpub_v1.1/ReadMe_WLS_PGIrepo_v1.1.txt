*****************************************************************************
**                             PGI REPOSITORY                              **
*****************************************************************************

    #-------------------------------------------------------------------#
    #                       VERSION : 1.1                               #
    #                        AUTHOR : Aysu Okbay                        #
    #                        E-MAIL : a.okbay@vu.nl                     #
    #                                                                   #
    #                      FILENAME : WLS_PGIrepo_v1.1.txt              #
    #                  RELEASE DATE : 14 June 2022                      #
    #-------------------------------------------------------------------#


+   DESCRIPTION
=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=

The file contains the following columns: 
    - FID                : Family identifier
    - IID                : Individual identifier
    - PGI_<pheno>-single : Single-trait PGI for phenotype <pheno> (see below)
    - PGI_<pheno>-multi  : Multi-trait PGI for phenotype <pheno> (see below)
    - PC1 - PC20         : Top 20 principal components of genotype data


Information about how the PGIs and PCs were made can be found in Becker et al. (2021). Please also take a look at the User Guide, where we lay out some of the interpretational issues that are likely to arise as researchers begin to use PGIs from the Repository, and outline how we suggest thinking through those issues. The User Guide can be found in Supplementary Information Section 7 of Becker et al. (2021).



+   v1.1 CHANGE LOG
=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=
The single-trait educational attainment PGI was updated. Instead of Lee et al. (2018; EA3), it is now based on the Okbay et al. (2022; EA4) GWAS excluding the following cohorts:
    - National Longitudinal Study of Adolescent to Adult Health (Add Health)
    - Estonian Bioabank (EGCUT)
    - English Longitudinal Study of Ageing (ELSA)
    - Health and Retirement Study (HRS)
    - Minnesota Center for Twin and Family Research (MCTFR)
    - Swedish Twin Registry (STR)
    - Wisconsin Longitudinal Study (WLS)
The sample size for this GWAS is 2,953,156. 
Moreover, following Okbay et al. (2022), this PGI was made using SBayesR (Lloyd-Jones et al., 2019) instead of LDpred (Vilhjalmsson et al., 2015). SBayesR is a Bayesian method that differs from LDpred in that it imposes a flexible finite mixture of normal distributions as the prior on the SNP effects instead of a point-normal mixture distribution. Like LDpred, SBayesR requires an estimate of LD between SNPs. We used the shrunk LD matrix from Lloyd-Jones et al. (2019) that was obtained in a sample ~50k UKB European-genetic-ancestry individuals for a set of 2,865,810 pruned common variants. We excluded the SNPs in the MHC region (Chr6 : 28-34Mb) from the analysis as recommended by Lloyd-Jones et al. (2019), as this was observed to improve model convergence. We ran SBayesR assuming 4 components in the finite mixture model, with initial mixture probabilities pi=(0.95,0.02,0.02,0.01) and fixed gamma=(0.0,0.01,0.1,1), where gamma is a parameter that constrains how the common SNP effect variance scales in each of the four distributions. The MCMC was run for 10,000 iterations with 2,000 taken as burn-in. 



+   PHENOTYPE ABBREVIATIONS
=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=

ACTIVITY        : Physical Activity
ADHD            : Attention Deficit Hyperactivity Disorder
ADVENTURE       : Adventurousness	
AFB             : Age First Birth
ALLERGYCAT      : Allergy - Cat	
ALLERGYDUST     : Allergy - Dust
ALLERGYPOLLEN   : Allergy - Pollen
ASTECZRHI       : Asthma/Eczema/Rhinitis	
ASTHMA          : Asthma
AUDIT           : Alcohol Misuse	
BMI             : Body Mass Index
CANNABIS        : Cannabis Use
COGEMP          : Cognitive Empathy
COPD            : Chronic Obstructive Pulmonary Disease
CPD             : Cigarettes per Day
CP              : Cognitive Performance
DELAYDISC       : Delay Discounting
DEP             : Depressive Symptoms
DPW             : Drinks per Week
EA              : Educational Attainment
EVERSMOKE       : Ever Smoker
EXTRA           : Exraversion
FAMSAT          : Life Satisfaction - Family
FINSAT          : Life Satisfaction - Finance
FRIENDSAT       : Life Satisfaction - Friend
HAYFEVER        : Hayfever (Allergic Rhinitis)
HEIGHT          : Height
HIGHMATH        : Highest Math	
LEFTOUT         : Left Out of Social Activity
LONELY          : Loneliness
MENARCHE        : Age First Menses
MIGRAINE        : Migraine
MORNING         : Morning Person
NARCIS          : Narcissism
NEARSIGHTED     : Nearsightedness
NEBmen          : Number Ever Born (men)
NEBwomen        : Number Ever Born (women)
NEURO           : Neuroticism
OPEN            : Openness
READING         : Childhood Reading
RELIGATT        : Religious Attendance
RISK            : Risk Tolerance
SELFHEALTH      : Self-Rated Health 
SELFMATH        : Self-Rated Math Ability
SWB             : Subjective Well-Being
VOICEDEEP       : Age Voice Deepened
WORKSAT         : Life Satisfaction - Work




+   Citation Instructions
=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=
In any publication that uses one or more Repository PGIs, please cite the published GWAS included in the single-trait or multi-trait input GWAS for the PGI as well as:
Becker, J., Burik, C. A. P. P., Goldman, G., Wang, N., Jayashankar, H., Bennett, M., Belsky, D. W., Karlsson Linner, R., Ahlskog, R., Kleinman, A., Hinds, D. A., Agee, M., Alipanahi, B., Auton, A., Bell, R. K., Bryc, K., Elson, S. L., Fontanillas, P., Furlotte, N. A., ... Okbay, A. (2021). Resource profile and user guide of the Polygenic Index Repository. Nature Human Behaviour. https://doi.org/10.1038/s41562-021-01119-3

You can find the published GWAS included in the single-trait or multi-trait input GWAS in Supplementary Tables 8 and 10 of Becker et al. (2021).



References:
Becker, J., Burik, C. A. P. P., Goldman, G., Wang, N., Jayashankar, H., Bennett, M., Belsky, D. W., Karlsson Linner, R., Ahlskog, R., Kleinman, A., Hinds, D. A., Agee, M., Alipanahi, B., Auton, A., Bell, R. K., Bryc, K., Elson, S. L., Fontanillas, P., Furlotte, N. A., ... Okbay, A. (2021). Resource profile and user guide of the Polygenic Index Repository. Nature Human Behaviour. https://doi.org/10.1038/s41562-021-01119-3

Lee, J. J., Wedow, R., Okbay, A., Kong, E., Maghzian, O., Zacher, M., Nguyen-Viet, T. A., Bowers, P., Sidorenko, J., Karlsson Linner, R., Fontana, M. A., Kundu, T., Lee, C., Li, H., Li, R., Royer, R., Timshel, P. N., Walters, R. K., Willoughby, E. A., ... Cesarini, D. (2018). Gene discovery and polygenic prediction from a genome-wide association study of educational attainment in 1.1 million individuals. Nature Genetics, 50(8), 1112-1121. https://doi.org/10.1038/s41588-018-0147-3

Lloyd-Jones, L. R., Zeng, J., Sidorenko, J., Yengo, L., Moser, G., Kemper, K. E., Wang, H., Zheng, Z., Magi, R., Esko, T., Metspalu, A., Wray, N. R., Goddard, M. E., Yang, J., & Visscher, P. M. (2019). Improved polygenic prediction by Bayesian multiple regression on summary statistics. Nature Communications, 10(1), 5086. https://doi.org/10.1038/s41467-019-12653-0

Okbay, A., Wu, Y., Wang, N., Jayashankar, H., Bennett, M., Nehzati, S. M., Sidorenko, J., Kweon, H., Goldman, G., Gjorgjieva, T., Jiang, Y., Hicks, B., Tian, C., Hinds, D. A., Ahlskog, R., Magnusson, P. K. E., Oskarsson, S., Hayward, C., Campbell, A., ... Young, A. I. (2022). Polygenic prediction of educational attainment within and between families from genome-wide association analyses in 3 million individuals. Nature Genetics 2022 54:4, 54(4), 437-449. https://doi.org/10.1038/s41588-022-01016-z

Vilhjalmsson, B. J., Yang, J., Finucane, H. K., Gusev, A., Lindstrom, S., Ripke, S., Genovese, G., Loh, P. R., Bhatia, G., Do, R., Hayeck, T., Won, H. H., Neale, B. M., Corvin, A., Walters, J. T. R., Farh, K. H., Holmans, P. A., Lee, P., Bulik-Sullivan, B., ... Price, A. L. (2015). Modeling Linkage Disequilibrium Increases Accuracy of Polygenic Risk Scores. American Journal of Human Genetics, 97(4), 576-592. https://doi.org/10.1016/j.ajhg.2015.09.001


