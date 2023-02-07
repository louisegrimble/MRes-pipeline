import pandas as pd
import numpy as np
###output files

#   1. Read in egg & koala dfs
### input files
egg_df = snakemake.input[0]
koala_df = snakemake.input[1]
egg_df = pd.read_csv(egg_df, sep='\t', header=0, low_memory=False)
koala_df = pd.read_csv(koala_df, sep='\t', low_memory=False)

##replace 0 to NaN 
egg_df.replace({'-': np.nan}, inplace= True)

############# !!!!! DROP the header rows pls
koala_df.drop(koala_df.index[0:1], inplace= True)
#    2. Print number of rows of proteins in both egg_df & koala_df
# minus the headers!
# remove headers- clean_egg_df = egg_df.drop(egg_df.index[0:1])
total_egg = len(egg_df)
total_egg

#same for koala_df
#remove headers- clean_koala_df = koala_df.drop(koala_df.index[0:1])
total_koala = len(koala_df)
total_koala


#    3. GhostKOALA has duplicates of the same protein.. remove duplicates
koala_df.drop_duplicates(subset=['gene name'], inplace=True)
koala_drop_2 = len(koala_df)
koala_drop_2

#   8. Count empty keggs (egg_df & koala_df)
count_egg = len(egg_df[pd.isnull(egg_df['KEGG_ko'])])

count_kofam = len(koala_df[pd.isnull(koala_df['KO'])])

#   9. Get % of empty KEGGs
egg_no = len(egg_df[pd.isnull(egg_df['KEGG_ko'])])
total_egg = len(egg_df)
egg_percent = (egg_no / total_egg)*100
egg_percent

koala_no = len(koala_df[pd.isnull(koala_df['KO'])])
total_koala = len(koala_df)
koala_percent = (koala_no / total_koala)*100
koala_percent

#####making these into a df
# initialize data of lists.
data = {'Q': ['Total predictions EggNOG', 'Total predictions Kofamscan', 'Total predictions Kofamscan (no duplicates)','Total proteins no KEGG (EggNOG)','Total proteins no KEGG (Kofamscan)','% proteins no KEGG (EggNOG)','% proteins no KEGG (Kofamscan)'],'A': [total_egg, total_koala, koala_drop_2, count_egg, count_kofam, egg_percent, koala_percent]}

stats_df = pd.DataFrame(data)
stats_df
stats_df.to_csv(snakemake.output[0], index=0)
################## these all print their own df ####################
#    4. print/count best tax levels (egg_df)
egg_df["max_annot_lvl"] = egg_df["max_annot_lvl"].apply(lambda x: x.split("|")[1])
egg_tax = egg_df.groupby(['max_annot_lvl'])['max_annot_lvl'].count()
egg_tax

egg_tax_df = pd.DataFrame(egg_tax)
egg_tax_df.to_csv(snakemake.output[1], index=0)


#    5. How many pathway modules? -> keggdecoder *
#this is not finished or correct
#but would be interesting
egg_module = egg_df.groupby(['KEGG_Module'])['KEGG_Module'].count()
egg_module_df = pd.DataFrame(egg_module)

egg_module_df.to_csv(snakemake.output[2], index=0)

#   6. How many COGs?
egg_cog = egg_df.groupby(['COG_category'])['COG_category'].count()
print("EggNOG output: Number of COGs involved:", egg_cog)

egg_cog_df = pd.DataFrame(egg_cog)
egg_cog_df.to_csv(snakemake.output[3], index=0)






#   10. Concat egg_df & koala_df -> no duplicates?
outer_name = pd.merge(egg_df, koala_df, how='outer', left_on='#query', right_on='gene name')
outer_name


## this is to show how to single out certain proteins for the extra annotation
#egg_df['KEGG_ko'].tolist()

#double_kegg = ['ko:K02444,ko:K03655']
#a = outer_name.loc[outer_name['KEGG_ko']== 'ko:K02444,ko:K03655']
#a['#query'].tolist()

###### need to make sure removing duplicates keeps the result with the best E-value
###### need to split up KEGGs that have a list of keggs on each one
###### how can I ask to match of the ones that match? look at line above


#    14. Format egg_df to run through decoder

######### need to regex in spyder instead of on regex101 website

#data = [egg_df['#query_name'], egg_df['KEGG_ko']]
#headers = ["protein_name", "KEGG_ko"]
#egg_decoder = pd.concat(data, axis=1, keys=headers)

########## need to remove 'ko:' pattern from KEGG_ko column in both egg_decoder df & outer_name df
#outer_copy = outer_name.copy()
#txt = str(outer_copy['KEGG_ko'][1485])
#x = txt.split("ko:")
#print(txt)
#print(x)


#    15. Confidence ranking using .apply()

#####
# Function that will set the confidence value
# EDIT THIS FOR MORE CONDITIONS
#####
def test_function(kofam, eggnog, cog):
  if str(kofam) in str(eggnog) and str(kofam) != '0' and str(eggnog) != '0': #egg/kofam=y match + cog
    x = 1
  elif str(eggnog) != '0' and str(kofam) != '0': #egg/kofam=y no match + cog
    x = 2
  elif str(kofam) != '0' and str(eggnog) in str(kofam) and str(cog) != '0': #egg/kofam=y + no cog
    x = 3
  elif str(eggnog) != '0' and str(kofam) in str(eggnog) and str(cog) != '0': #egg/kofam=y + no cog
    x = 3
  elif str(kofam) != '0' and str(eggnog) == '0' and str(cog) == '0': #kofam=y egg=n + no cog
    x = 4
  elif str(kofam) == '0' and str(eggnog) != '0' and str(cog) == '0': #kofam=y egg=n + no cog
    x = 4
  elif str(cog) != '0':
    x = 5
  else:
    x = 6

  return x

subset = outer_name.copy() # creating a test subset - currently the whole outer_name df
subset.rename(columns = {'COG_category':'functional_cog'}, inplace = True)
subset['functional_cog'].fillna(0, inplace=True)
subset['KEGG_ko'].fillna(0, inplace=True)
subset['KO'].fillna(0, inplace=True)
# This is applying the function to each row (quickly! - vectorisation)
# Function returns X which is your confidence score, and sets this to the column "value"
#
subset["confidence_value"] = subset.apply(lambda x: test_function(x.KO, x.KEGG_ko, x.functional_cog), axis = 1)
subset.groupby('confidence_value').count()
subset['confidence_value'].value_counts()
subset

subset.to_csv(snakemake.output[4], index=0)

######################### to make new df for low qual #######################

confidence_count = subset['confidence_value'].value_counts()
confidence_df = pd.DataFrame(confidence_count)
confidence_df.to_csv(snakemake.output[5], index=0)

#subset.loc[subset['confidence_value'] > '3']['#query'] ## or like this??
low_qual = subset.loc[subset['confidence_value'] > 3]['#query']
low_qual_df = pd.DataFrame(low_qual)
low_qual_df.to_csv(snakemake.output[6], index=0)


