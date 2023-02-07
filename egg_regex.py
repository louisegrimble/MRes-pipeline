import re
import pandas as pd

egg_df = snakemake.input[0]
egg_df = pd.read_csv(egg_df, sep='\t', header=0)
egg_df.drop(egg_df.index[0], inplace= True) #remove header??

##regex to remove "ko:" from each row
egg_df.apply


egg_df['KEGG_ko'].fillna('', inplace= True)
egg_index = egg_df['#query'].to_list()
egg_index
r1_egg_167 = [re.sub('[k][o].','', str(x)) for x in egg_df['KEGG_ko']]
##regex to remove any extra keggs
r2_egg_167= [re.sub('[,].*$','', str(x)) for x in r1_egg_167]
r2_egg_167
data2 = {'query':egg_index, 'KEGG':r2_egg_167}
egg_kd_167 = pd.DataFrame(data2)
### will need to rename from r2 to something like egg_decoder_167
egg_kd_167.to_csv(snakemake.output[0], index=0)
