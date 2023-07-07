f = open ("template_GenerationMixturefewersignaturespairedPCAWG.R", "r")
a = f.read()
substitute='Stomach-AdenoCA'

f2 = open ("template_GenerationMixturefewersmallsignaturespairedPCAWG.R", "r")
a2 = f2.read()
substitute2='Kidney-RCC.clearcell'

f3 = open ("template_GenerationMixturefewersignaturespairedObsNmPCAWG.R", "r")
a3 = f3.read()

f4 = open ("template_GenerationMixtureallsignaturespairedObsNmPCAWG.R", "r")
a4 = f4.read()

f5 = open ("template_GenerationMixturefewersignaturespairedObsNmObsDMPCAWG.R", "r")
a5 = f5.read()

f7 = open ("template_GenerationMixturefewersignaturespairedObsNmPoissonSigPCAWG.R", "r")
a7 = f7.read()

f8 = open ("template_GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWG.R", "r")
a8 = f8.read()

f9 = open ("template_GenerationMixturefewersignaturespairedObsNmGaussianPCAWG.R", "r")
a9 = f9.read()

f10 = open ("template_GenerationMixturefewersignaturespairedObsNmPoissonUnifPCAWG.R", "r")
a10 = f10.read()

f11 = open ("template_GenerationMixturefewersignaturespairedObsNmInvPCAWG.R", "r")
a11 = f11.read()


dic = {"BoneOsteosarcPCAWG": "Bone-Osteosarc", "BreastAdenoCAPCAWG": "Breast-AdenoCA", "CNSGBMPCAWG": "CNS-GBM",\
"CNSMedulloPCAWG": "CNS-Medullo", "CNSPiloAstroPCAWG": "CNS-PiloAstro", "ColoRectAdenoCAPCAWG": "ColoRect-AdenoCA",\
"EsoAdenoCAPCAWG": "Eso-AdenoCA", "HeadSCCPCAWG": "Head-SCC", "KidneyChRCCPCAWG": "Kidney-ChRCC",\
"KidneyRCCPCAWG": "Kidney-RCC.clearcell", "KidneyRCCpapillaryPCAWG":  "Kidney-RCC.papillary",\
"LiverHCC": "Liver-HCC", "LungSCCPCAWG": "Lung-SCC", "LymphBNHLPCAWG":  "Lymph-BNHL", "LymphCLLPCAWG":  "Lymph-CLL",\
"OvaryAdenoCAPCAWG": "Ovary-AdenoCA", "PancAdenoCAPCAWG": "Panc-AdenoCA", "PancEndocrinePCAWG": "Panc-Endocrine",\
"ProstAdenoCAPCAWG": "Prost-AdenoCA", "SkinMelanomacutaneousPCAWG": "Skin-Melanoma.cutaneous",\
"StomachPCAWG": "Stomach-AdenoCA", "ThyAdenoCAPCAWG": "Thy-AdenoCA", "UterusAdenoCAPCAWG": "Uterus-AdenoCA" }

for i in dic.keys():
    # o = open ("GenerationMixturefewersignaturespaired"+i+".R", "w")
    # o.write(a.replace(substitute, dic[i]))
    # o.close()

    # o2 = open ("../../../../data/assessing_models_simulation/GenerationMixturefewersignaturespaired"+i, "w")
    # o2.write('Generated automatically using /Users/morril01/Documents/PhD/GlobalDA/code/3_analysis/simulation_model_assessment/generate_datasets_simulations/create_paired_pcawg_datasets.py\n')
    # o2.close()

    # ## signatures of lowest abundance
    # o3 = open ("GenerationMixturefewersmallsignaturespaired"+i+".R", "w")
    # o3.write(a2.replace(substitute2, dic[i]))
    # o3.close()

    ## with number of mutations from observations
    # o4 = open ("GenerationMixturefewersignaturespairedObsNm"+i+".R", "w")
    # o4.write(a3.replace(substitute, dic[i]))
    # o4.close()

    # o4 = open ("../../../../data/assessing_models_simulation/GenerationMixturefewersignaturespairedObsNm"+i, "w")
    # o4.write('Generated automatically using /Users/morril01/Documents/PhD/GlobalDA/code/3_analysis/simulation_model_assessment/generate_datasets_simulations/create_paired_pcawg_datasets.py\n')
    # o4.close()

    # o5 = open ("GenerationMixtureallsignaturespairedObsNm"+i+".R", "w")
    # o5.write(a4.replace(substitute, dic[i]))
    # o5.close()

    # o5 = open ("../../../../data/assessing_models_simulation/GenerationMixtureallsignaturespairedObsNm"+i, "w")
    # o5.write('Generated automatically using /Users/morril01/Documents/PhD/GlobalDA/code/3_analysis/simulation_model_assessment/generate_datasets_simulations/create_paired_pcawg_datasets.py\n')
    # o5.close()

    # o6 = open ("GenerationMixturefewersignaturespairedObsNmObsDM"+i+".R", "w")
    # o6.write(a5.replace(substitute, dic[i]))
    # o6.close()

    # o6 = open ("../../../../data/assessing_models_simulation/GenerationMixturefewersignaturespairedObsNmObsDM"+i, "w")
    # o6.write('Generated automatically using /Users/morril01/Documents/PhD/GlobalDA/code/3_analysis/simulation_model_assessment/generate_datasets_simulations/create_paired_pcawg_datasets.py\n')
    # o6.close()

    for name in ['GenerationMixturefewersignaturespairedObsNmGaussianPCAWG', 'GenerationMixturefewersignaturespairedObsNmGaussianVarPCAWG',\
    'GenerationMixturefewersignaturespairedObsNmPoissonSigPCAWG', 'GenerationMixturefewersignaturespairedObsNmPoissonUnifPCAWG',\
    'GenerationMixturefewersignaturespairedObsNmInvPCAWG']:
        f = open ("template_"+name+".R", "r")
        a = f.read()
        o = open (name+i+".R", "w")
        o.write(a.replace(substitute, dic[i]))
        o.close()

        o = open ("../../../../data/assessing_models_simulation/"+name+i, "w")
        o.write('Generated automatically using /Users/morril01/Documents/PhD/GlobalDA/code/3_analysis/simulation_model_assessment/generate_datasets_simulations/create_paired_pcawg_datasets.py\n')
        o.close()

        if name == 'GenerationMixturefewersignaturespairedObsNmInvPCAWG':
            print("../../../../data/assessing_models_simulation/template_"+name+i)

print(dic.keys())
