import pandas as pd 
from scipy.stats import fisher_exact
import numpy as np 
import os 
from itertools import chain


from utils import get_combinations


def get_n00(s1,s2,t):
    s1 = {i.split('--&&--')[0] for i in s1}
    s2 = {i.split('--&&--')[0] for i in s2}

    return len(t - s1.union(s2))


def get_nud(s1,s2):

    s1 = [i.split('--&&--')[0] for i in s1]
    s2 = [i.split('--&&--')[0] for i in s2]

    s1_f = [i for i in s1 if i not in s2 ]
    s2_f = [i for i in s2 if i not in s1]

    for i in s2_f:
        s1_f.append(i)

    return len(set(s1_f))

def get_ndu(s1,s2):
    list_of_dicts = []

    for s1s in s1:
        c1,e1 = s1s.split('--&&--')
        for s2s in s2:
            c2,e2 = s2s.split('--&&--')
            if (c1 == c2) and (e1 == "Up-regulated") and (e2 == "down-regulated"):
                    list_of_dicts.append(c1)

    for s1s in s1:
        c1,e1 = s1s.split('--&&--')
        for s2s in s2:
            c2,e2 = s2s.split('--&&--')
            if (c1 == c2) and (e1 == "down-regulated") and (e2 == "Up-regulated"):
                list_of_dicts.append(c1)

    return len(set(list_of_dicts))


def get_uudd(s1,s2):
    total_matches = sum(pair in s2 for pair in s1)
    return total_matches

def get_n(s1,s2,df, total):

    s1_con = set(df.at[s1,'condition_exp'])
    s2_con = set(df.at[s2,'condition_exp'])

    n_00 = get_n00(s1_con,s2_con,total)

    n_ud = get_nud(s1_con,s2_con)
    n_du = get_ndu(s1_con,s2_con)
    n_uudd = get_uudd(s1_con,s2_con)
    
    data = np.array([[n_00,n_ud],[n_du,n_uudd]])

    _ ,p_value = fisher_exact(data, alternative = 'greater')
    
    return n_00, n_ud, n_du, n_uudd, p_value

def generate_matrix_dif(df, CUT_OFF):
    
    new_total_ex = set(list(df['exp_condition'].unique()))

    df['condition_exp'] = df[['exp_condition', 'expression']].apply(lambda x: '--&&--'.join(map(str, x)), axis=1)


    df = df.groupby('mapped_phosphosite').agg(pd.Series.tolist).reset_index()

    
    all_sites = df['mapped_phosphosite'].reset_index()

    df.set_index("mapped_phosphosite", inplace = True) 

    df_comb = get_combinations(all_sites)


    df_comb[['n_00','n_uddu_nd','n_uddu','n_uudd','p-Value']] = df_comb.apply(lambda x:get_n(x['site1'],x['site2'],
                                                                                                          df , new_total_ex),  axis = 1, result_type='expand')

    df_comb.drop_duplicates(inplace=True)

    df_comb.sort_values(by=['p-Value'], inplace=True)

    return df_comb
