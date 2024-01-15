## Copy and pasted from https://github.com/lrgr/tcsm/blob/master/src/estimate_significance.py
## but making the input files arguments

import pandas as pd, numpy as np
import math
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
import argparse

parser = argparse.ArgumentParser(description="Run estimate_significance.py",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--sigma")
parser.add_argument("-g", "--gamma")
# parser.add_argument("-f", "--feature")
parser.add_argument("-c", "--covariates")
parser.add_argument("-C", "--covariate_of_interest")
parser.add_argument("-o", "--output")
args = parser.parse_args()
config = vars(args)
print(config)

def simulate_logistic_normal(mean, cov, n_draws):
    y = np.random.multivariate_normal(mean, cov, size=n_draws)
    # now figure out how to convert from a multivariate_normal distribution to a logistic_normal distribution
    y = pd.DataFrame(data=y, columns=["Topic{}".format(i) for i in range(1, np.size(mean)+1)])
    y["Topic{}".format(np.size(mean)+1)] = 0
    # first take exp(e) for every element e in the data frame
    y = y.applymap(math.exp)
    # normalize each element by the sum of its row
    return y.divide(y.sum(axis=1), axis='index')



# read input files
sigma_df = pd.read_csv(config['sigma'], sep="\t", index_col=0)
gamma_df = pd.read_csv(config['gamma'], sep="\t", index_col=0)
# feature_df = pd.read_csv(config['feature'], sep="\t", index_col=0)
# set parameters
n_draws = 5000
# set the seed for reproducibility
np.random.seed(int(1234))
# assume that there is only a single binary covariate that we are interested in
covariate_of_interest = config["covariate_of_interest"]
# identify if any other covariates exist (because if so, we need to sample them)
covariate_list = config['covariates'].split("+")
other_covariates = list(set(covariate_list).difference([covariate_of_interest]))
assert len(other_covariates) == 0
cov = sigma_df.values
mean = gamma_df.loc["default"].values
default_exposures = simulate_logistic_normal(mean, cov, n_draws)
mean = gamma_df.loc["default"].values + gamma_df.loc[covariate_of_interest].values
covariate_exposures = simulate_logistic_normal(mean, cov, n_draws)
topics = ["Topic{}".format(i) for i in range(1, np.size(mean, axis=0)+2)]
output = []
for topic in topics:
    ranksum_result = mannwhitneyu(x=covariate_exposures[topic], y=default_exposures[topic], alternative='greater')
    data = [ranksum_result[0], ranksum_result[1], np.mean(covariate_exposures[topic]), np.mean(default_exposures[topic])]
    s = pd.Series(data=data, index=["Mann-Whitney U Statistic", "Mann-Whitney U p-value", "covariate mean", "default mean"])
    output.append(s)
df = pd.concat(output, axis=1).transpose()
# calculate the FDR
df["Wilcoxon FDR"] = fdrcorrection(df["Mann-Whitney U p-value"])[1]
df.to_csv(config['output'], sep="\t")

