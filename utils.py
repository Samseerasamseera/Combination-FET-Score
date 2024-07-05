import pandas as pd 
from scipy.stats import fisher_exact
import numpy as np 
import os 
import pymysql
from itertools import combinations
import itertools


import time


# DB_HOST = os.getenv('ciodsdb.cubtgfved0u5.us-west-1.rds.amazonaws.com')
# DB_USER = os.getenv('admin')
# DB_PASSWORD = os.getenv('ciods123')
# DB_NAME = os.getenv('Phosphoproteome_database')


def get_query(k):
    q = f'''SELECT exp_condition , mapped_phosphosite FROM Filter_profile_with_id WHERE mapped_gene = '{k}' '''
    return q

def get_query_d(k,e):

    q = f'''SELECT exp_condition , mapped_phosphosite,expression FROM phospodb_nisar_differential_data WHERE mapped_genesymbol ='{k}'  '''

    return q


def get_d(query):
    connection = pymysql.connect(
        host='ciodsdb.cubtgfved0u5.us-west-1.rds.amazonaws.com',
        port=3306,
        user='admin',
        password='ciods123',
        database='Phosphoproteome_database'

    )
    cursor = connection.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    df_result = pd.DataFrame(result, columns=[desc[0] for desc in cursor.description])
    cursor.close()
    connection.close()
    # df_result.to_excel('jj.xlsx')
    return df_result
    


# def get_combinations(df):

#     pairs_list = []

#     for i in range(len(df)):
#         current_data = df.at[i, 'mapped_phosphosite']
#         remaining_data = df['mapped_phosphosite'].iloc[i+1:].tolist()
        
#         pairs = [(current_data, item) for item in remaining_data]
        
#         pairs_list.extend(pairs)

#     result_df = pd.DataFrame(pairs_list, columns=['site1', 'site2'])

#     mask = result_df['site1'] != result_df['site2']

#     result_df = result_df[mask]

#     result_df
#     print(result_df)
    
#     return result_df







def get_combinations(df):
    pairs_list = []

    
    phosphosites = df['mapped_phosphosite'].tolist()
    print(phosphosites)
    rphosphosites = list(reversed(phosphosites))
    print(rphosphosites)

    
    pairs = list(combinations(phosphosites, 2))
    print(pairs)
    rpairs = list(combinations(rphosphosites, 2))
    print(rpairs)

    pairs = pairs + rpairs

    pairs = list(set(pairs))


    
    result_df = pd.DataFrame(pairs, columns=['site1', 'site2'])
    result_df['sorted'] = result_df.apply(lambda row: sorted([row['site1'], row['site2']]), axis=1)
    result_df = result_df.drop_duplicates(subset='sorted').drop(columns='sorted')


    
    print(result_df)
    

    return result_df


def get_n00(s1,s2,t):
    return list(t - (s1.union(s2))),len(t - (s1.union(s2)))

def get_n01(s1,s2,t):
    return list(s1.intersection(s2)),len(s2 - (s1.intersection(s2)))

def get_n10(s1,s2,t):
    return list(s1 - (s1.intersection(s2))),len(s1 - (s1.intersection(s2)))

def get_n11(s1,s2,t):
    return list(s1.intersection(s2)), len(s1.intersection(s2))

def get_n(s1,s2,df, total):
    s1_con = df.at[s1,'exp_condition']
    s2_con = df.at[s2,'exp_condition']

    s1_con = set(s1_con)
    s2_con = set(s2_con)
    total = set(total)

    n00_exp,n_00 = get_n00(s1_con,s2_con,total) 
    n01_exp,n_01 = get_n01(s1_con,s2_con,total)
    n10_exp,n_10 = get_n10(s1_con,s2_con,total)
    n11_exp,n_11 = get_n11(s1_con,s2_con,total)

    n00_exp = ';'.join(n00_exp)
    n01_exp = ';'.join(n01_exp)
    n10_exp = ';'.join(n10_exp)
    n11_exp = ';'.join(n11_exp)

    data = np.array([[n_00,n_01],[n_10,n_11]])
    _ ,p_value = fisher_exact(data, alternative = 'greater')
    return n00_exp,n01_exp,n10_exp,n11_exp,n_00,n_01,n_10,n_11 , p_value

def generate_matrix(df, CUT_OFF):
    
    total_ex = df["exp_condition"].unique().tolist() 
    df = df.groupby('mapped_phosphosite').agg(pd.Series.tolist).reset_index()

    df['count'] = df["exp_condition"].apply(lambda x:len(x))
    df = df.sort_values(by=['count'], ascending=False)


    df = df.loc[df['count'] >= CUT_OFF]

    all_sites = pd.DataFrame(df['mapped_phosphosite'])
    
    df.set_index('mapped_phosphosite', inplace = True)

    df_comb = get_combinations(all_sites)
    

    df_comb[['n00_exp','n01_exp','n10_exp','n11_exp','n_00','n_01','n_10','n_11','p-Value']] = df_comb.apply(lambda x:get_n(x['site1'],x['site2'],df , total_ex),  axis = 1, result_type='expand')

    # result = df_comb.apply(lambda x:get_n(x['site1'],x['site2'], df , total_ex), axis = 1)

    # df_comb[['n_00','n_01','n_10','n_11','p-Value']] = pd.DataFrame(result.tolist())

    df_comb.drop_duplicates(inplace=True)
    df_comb.sort_values(by=['p-Value'], inplace=True)

    return df_comb




def getCDF(df):
    df['CDF'] = np.arange(1, len(df) + 1) / len(df)
    return df

def get_trypsin(s1,s2, df):
    s1 = int(s1[1:])
    s2 = int(s2[1:])

    peptides = []
    for i,row in df.iterrows():
        if (s1 in row['counts']) and (s2 in row['counts']):
            return [row['Fragments']]
        
        elif s1 in row['counts']:
            peptides.append(row['Fragments'])

        elif s2 in row['counts']:
            peptides.append(row['Fragments'])

        else:
            continue
    
    return peptides
    

# def trypsin_digest(df):
    
#     df_tr = pd.read_excel("trypsin_digestion_PAK1.xlsx")
#     main_count = []
#     couts = []
#     c = 1
#     for i, row in df_tr.iterrows():
#         for _ in range (0, len(row['Fragments'])):
#             couts.append(c)
#             c += 1
#         main_count.append(couts)
#         couts = []
#     df_tr["counts"] = main_count

#     df['peptides'] = df.apply(lambda x:get_trypsin(x['site1'], x['site2'], df_tr), axis =1)

#     df["same_peptide"] = df['peptides'].apply(lambda x:"YES" if len(x) ==1 else "NO")

#     return df






















