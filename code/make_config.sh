#!/bin/bash

name_config="config_PCAWG_do_not_modify.yaml"
name_failed_samples="../text/readme/1_faulty_samples"
echo "Creating config file $name_config ... "

echo "Samples for which there is no VCF file\nThese files appear in the metadata and may have them in the mutccf file, but I don't have their VCF, which is the only file that contains what mutation it is (in mutccf you can have the position and CCF, but not mutation type).\n\n `````` " > $name_failed_samples


## Samples
echo "samples:" > $name_config
flder="../data/restricted/pcawg/pcawg_restricted_snv/"

for i in ${flder}*gz; do
	bs=`(basename $i)`
	echo "    " ${bs::36} ": " $bs >>  $name_config
done


## Sample groups
echo "\nsample_groups:" >>  $name_config
awk '{print "- "$2}'  ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt  | tail -n +2 | sort -u >>  $name_config


## Types of features
echo "\nfeature_types:
-    signatures
-    nucleotidesubstitution1
-    nucleotidesubstitution3" >> $name_config

## To which cancer type samples belong
printf "\n%s" "grouped_samples:" >>  $name_config
for i in `awk '{print $2}'  ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt  | tail -n +2 | sort -u`; do  printf "\n    %s : " "$i" >> $name_config; 
flenames=`grep $i ../data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt | awk '{printf "../data/restricted/pcawg/pcawg_restricted_snv_counts/"$1 " "}'`
for flename in $flenames; do
	nme="../data/restricted/pcawg/pcawg_restricted_snv/"$(basename $flename)".consensus.20160830.somatic.snv_mnv.vcf.gz"
echo $nme
## check if file exists. if it doesn't append to some extra file containing failed samples
if [ -f $nme ]; then
	printf "%s " $flename >> $name_config
else
	echo $flename >> $name_failed_samples
fi
done
done


## Put the simulation results
echo "\n\nsimulation_n:
- 10
- 20
- 50
- 100
- 200
" >> $name_config

echo "simulation_nlambda:
- 30
- 100
- 300
" >> $name_config

echo "simulation_d:
- 3
- 7
" >> $name_config

echo "beta_gamma_shape:
- 0
- 0.02
- 0.1
- 1
- 2.5
" >> $name_config

## Generation C of simulations
echo "\nsimulation_n_generationC:
- 10
- 20
- 50
- 100
" >> $name_config

echo "beta_gamma_shape_generationC:
- 0
- 0.0001
- 0.0002
- 0.0005
- 0.00075
- 0.001
- 0.002
- 0.005
- 0.001
- 0.002
- 0.01
- 0.02
- 0.05
- 0.1
- 0.2
- 0.4
- 0.6
- 0.8
- 2
- 4
" >> $name_config

echo "simulation_d_generationC:
- 3
- 4
- 5
- 6
- 7
" >> $name_config

echo "simulation_nlambda_generationC:
- 100
" >> $name_config

echo "simulation_lambda_generationC:
- 10
- 80
" >> $name_config


## Generation D of simulations
echo "\nsimulation_n_generationD:
- 20
- 25
- 30
- 35
- 40
- 45
" >> $name_config

echo "beta_gamma_shape_generationD:
- 0
- 1
" >> $name_config

echo "simulation_d_generationD:
- 3
- 7
" >> $name_config

echo "simulation_nlambda_generationD:
- 100
" >> $name_config

echo "simulation_lambda_generationD:
- 10
- 80
" >> $name_config


## Generation E of simulations
echo "\nsimulation_n_generationE:
- 10
- 20
- 40
- 80
" >> $name_config

echo "beta_gamma_shape_generationE:
- 0
- 1
- 2
" >> $name_config

echo "simulation_d_generationE:
- 4
- 6
- 10
" >> $name_config

echo "simulation_nlambda_generationE:
- 200
" >> $name_config

echo "simulation_lambda_generationE:
- NA
" >> $name_config

echo "... Config file $name_config created."
echo "\`\`\`" >> $name_failed_samples


## Generation F of simulations
echo "\nsimulation_n_generationF:
- 20
- 50
- 100
" >> $name_config

echo "beta_gamma_shape_generationF:
- 0
- 0.01
- 0.05
- 0.1
- 0.2
- 0.4
- 0.6
- 0.8
- 2
- 4
" >> $name_config

echo "simulation_d_generationF:
- 3
- 4
- 5
- 6
- 7
" >> $name_config

echo "simulation_nlambda_generationF:
- 100
" >> $name_config

echo "simulation_lambda_generationF:
- 10
- 30
- 50
- 80
" >> $name_config
 
##########################################################################################
##########################################################################################

## Generation G of simulations
echo "\nsimulation_n_generationG:
- 20
- 30
- 50
- 100
" >> $name_config

echo "beta_gamma_shape_generationG:
- 0
- 0.01
- 0.1
- 0.6
- 4
" >> $name_config

echo "simulation_d_generationG:
- 3
- 4
- 5
- 6
" >> $name_config

echo "simulation_nlambda_generationG:
- 100
" >> $name_config

echo "simulation_lambda_generationG:
- 80
" >> $name_config

##########################################################################################
##########################################################################################

## Generation Mixture1 of simulations
echo "\nsimulation_n_generationMixture1:
- 30
- 50
- 100
" >> $name_config

echo "beta_gamma_shape_generationMixture1:
- 0
- 0.1
- 0.3
" >> $name_config

echo "simulation_d_generationMixture1:
- 3
- 4
- 5
" >> $name_config

echo "simulation_nlambda_generationMixture1:
- 100
" >> $name_config

echo "simulation_lambda_generationMixture1:
- 0
" >> $name_config

##########################################################################################
##########################################################################################

## Generation GenerationJnormBTwoLambdasOneChangingBeta of simulations
echo "\nsimulation_n_GenerationJnormBTwoLambdasOneChangingBeta:
- 20
- 30
- 50
- 100
" >> $name_config

echo "beta_gamma_shape_GenerationJnormBTwoLambdasOneChangingBeta:
- 0
- 0.01
- 0.1
- 0.2
- 0.3
- 0.6
- 4
" >> $name_config

echo "simulation_d_GenerationJnormBTwoLambdasOneChangingBeta:
- 3
- 4
- 5
- 6
" >> $name_config

echo "simulation_nlambda_GenerationJnormBTwoLambdasOneChangingBeta:
- 100
" >> $name_config

echo "simulation_lambda_GenerationJnormBTwoLambdasOneChangingBeta:
- 80
" >> $name_config

##########################################################################################
#### single change, poisson

echo "\nsimulation_n_GenerationPois:
- 50
- 100
" >> $name_config

echo "beta_gamma_shape_GenerationPois:
- 0
- 0.01
- 0.1
- 0.2
- 0.3
- 0.6
- 4
" >> $name_config

echo "simulation_d_GenerationPois:
- 3
- 4
- 5
- 6
" >> $name_config

echo "simulation_nlambda_GenerationPois:
- 100
" >> $name_config

echo "simulation_lambda_GenerationPois:
- 80
" >> $name_config


##########################################################################################
##########################################################################################

echo "\nsimulation_n_GenerationMixturePCAWG:
- 50
- 100
- 200
" >> $name_config

echo "beta_gamma_shape_GenerationMixturePCAWG:
- -12
- -10
- -8
- -4
- -2
" >> $name_config

echo "simulation_d_GenerationMixturePCAWG:
- 7
" >> $name_config

echo "simulation_nlambda_GenerationMixturePCAWG:
- 200
" >> $name_config

echo "simulation_lambda_GenerationMixturePCAWG:
- 80
" >> $name_config

##########################################################################################
##########################################################################################

echo "\nsimulation_n_GenerationMixturefewersignaturesPCAWG:
- 50
- 100
- 200
" >> $name_config

echo "beta_gamma_shape_GenerationMixturefewersignaturesPCAWG:
- -999
- -12
- -10
- -8
- -4
- -2
" >> $name_config

echo "simulation_d_GenerationMixturefewersignaturesPCAWG:
- 4
" >> $name_config

echo "simulation_nlambda_GenerationMixturefewersignaturesPCAWG:
- 200
" >> $name_config

echo "simulation_lambda_GenerationMixturefewersignaturesPCAWG:
- 80
" >> $name_config

##########################################################################################
##########################################################################################

echo "\nsimulation_n_GenerationMixturefewersignaturespairedstomachPCAWG:
- 50
- 100
- 200
" >> $name_config

echo "beta_gamma_shape_GenerationMixturefewersignaturespairedstomachPCAWG:
- -999
- -12
- -10
- -8
- -4
- -2
- -1
- 0
- 4
- 8
- 12
" >> $name_config

echo "simulation_d_GenerationMixturefewersignaturespairedstomachPCAWG:
- 4
" >> $name_config

echo "simulation_nlambda_GenerationMixturefewersignaturespairedstomachPCAWG:
- 200
" >> $name_config

echo "simulation_lambda_GenerationMixturefewersignaturespairedstomachPCAWG:
- 80
" >> $name_config


##########################################################################################
##########################################################################################

echo "\nsimulation_n_GenerationMixturePCAWG2:
- 50
- 100
- 200
" >> $name_config

echo "beta_gamma_shape_GenerationMixturePCAWG2:
- -999
- -12
- -10
- -8
- -4
- -2
- -1
- 0
- 4
- 8
" >> $name_config

echo "simulation_d_GenerationMixturePCAWG2:
- 4
" >> $name_config

echo "simulation_nlambda_GenerationMixturePCAWG2:
- 200
" >> $name_config

echo "simulation_lambda_GenerationMixturePCAWG2:
- 80
" >> $name_config



