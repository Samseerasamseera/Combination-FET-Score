import pandas as pd


df = pd.read_excel("site example.xlsx", header=None, names=['Data'])

pairs_list = []

for i in range(len(df)):
    current_data = df.at[i, 'mapped_phosphosite']
    remaining_data = df['mapped_phosphosite'].iloc[i+1:].tolist()
    
    pairs = [(current_data, item) for item in remaining_data]
    
    pairs_list.extend(pairs)

result_df = pd.DataFrame(pairs_list, columns=['site', 'site_pair'])


