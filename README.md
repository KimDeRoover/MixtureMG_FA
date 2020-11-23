# MixtureMG_FA: R-code for mixture multigroup factor analysis

# MixtureMG_FA_loadings: For finding clusters of groups with metric invariance

Please cite https://www.researchgate.net/publication/343392110_Mixture_multigroup_factor_analysis_for_unraveling_factor_loading_non-invariance_across_many_groups

This version deals with factor loading differences (i.e., finds clusters of groups based on similarity of their factor loadings, given the user-specified number of clusters) and uses EFA.

For model selection, it is advised to use BIC_G (number of groups as sample size) in combination with CHull (see paper, https://www.rdocumentation.org/packages/multichull/versions/1.0.0)

# MixtureMG_FA_intercepts: For finding clusters of groups with scalar invariance

Please cite https://psyarxiv.com/5yr68

This version deals with intercept differences (i.e., finds clusters of groups based on similarity of their intercepts, given the user-specified number of clusters) and can be used with EFA as well as CFA. It builds on invariance of the loadings. Deal with loading non-invariances as described in the Discussion of the paper on the metric invariance variant.

For model selection, it is advised to use BIC_G (number of groups as sample size) in combination with CHull (see paper, https://www.rdocumentation.org/packages/multichull/versions/1.0.0)
