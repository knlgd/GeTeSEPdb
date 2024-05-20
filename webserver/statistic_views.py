import ast
import json
import os
import pickle
from django.shortcuts import render
from timecourse_app.models import Sampleinfo, clustergene, ModelDF
from django.db.models import Q, Count
import pandas as pd
import numpy as np
from collections import Counter
from django.conf import settings
import re


def get_statistic_page(request):
    total_study_cnt = Sampleinfo.objects.all().count()
    # Number of perturbation datasets in each model organism
    rnaseq_strategy = pd.DataFrame(Sampleinfo.objects.filter(Q(strategy='RNA-Seq')).
                                   values('organism').order_by('organism').
                                   annotate(count=Count('organism'))).rename(columns={'organism':'organism','count':'rnaseq'})

    microarray_strategy = pd.DataFrame(Sampleinfo.objects.filter(Q(strategy='Microarray')).
                                     values('organism').order_by('organism').
                                     annotate(count=Count('organism'))).rename(
        columns={'organism': 'organism', 'count': 'microarray'})

    # tmp = pd.merge(rnaseq_strategy, scrnaseq_strategy, how='outer')
    all = pd.merge(rnaseq_strategy, microarray_strategy, how='left', on='organism')
    all['sum'] = all.iloc[:, 1:].sum(axis=1)
    all = all.fillna(0).sort_values('sum', ascending=False)

    # rnaseq_strategy_cnt_k,rnaseq_strategy_cnt_v = dict2list(rnaseq_strategy_cnt)
    #
    # scrnaseq_strategy = GroupList.objects.filter(Q(library_strategy__library_strategy='scRNA-Seq')).values_list(
    #     'organism__organism', flat=True)
    # scrnaseq_strategy_cnt = Counter(scrnaseq_strategy)
    # scrnaseq_strategy_cnt_k,scrnaseq_strategy_cnt_v = dict2list(scrnaseq_strategy_cnt)
    #
    # microarray_strategy = GroupList.objects.filter(Q(library_strategy__library_strategy='Microarray')).values_list(
    #     'organism__organism', flat=True)
    # microarray_strategy_cnt = Counter(microarray_strategy)
    # microarray_strategy_k, microarray_strategy_v = dict2list(microarray_strategy_cnt)
    # make trace
    # organism_rank = pd.DataFrame.from_dict((rnaseq_strategy_cnt + scrnaseq_strategy_cnt + microarray_strategy_cnt),
    #                                        orient='index',columns=['count']).sort_values('count',ascending=False)
    # all_organism_k,all_organism_v = dict2list(rnaseq_strategy_cnt + scrnaseq_strategy_cnt + microarray_strategy_cnt)
    # rnaseq_strategy_cnt = pd.DataFrame.from_dict(rnaseq_strategy_cnt, orient='index', columns=['count'])
    barplot_organism_strategy = [
        {
            'x': list(all['organism']),
            'y': list(all['rnaseq']),
            'text': list(all['sum']),
            'textposition': 'none',
            'name': 'RNA-Seq',
            'type': 'bar'
        },
        {
            'x': list(all['organism']),
            'y': list(all['microarray']),
            'text': list(all['sum']),
            'textposition': 'outside',
            'name': 'Microarray',
            'type': 'bar'
        }
    ]
    # Proportion of data type for all datasets
    data_type = Sampleinfo.objects.all().values_list('strategy', flat=True)
    data_type_cnt = Counter(data_type)
    # make trace
    pie_data_type = [{
        'values': list(data_type_cnt.values()),
        'labels': list(data_type_cnt.keys()),
        'domain': {'column': 0},
        'name': 'Data Type',
        'hoverinfo': 'label+percent+name',
        'hole': .4,
        'type': 'pie'
    }]
    # Proportion of model type for all datasets
    # progress_model_type = ModelDF.objects.all().values_list('model', flat=True)
    # constant_model_type = clustergene.objects.filter(cluster='cluster0').values_list('model', flat=True)
    #
    # progress_model_type_cnt = Counter(progress_model_type)
    # constant_model_type_cnt = Counter(constant_model_type)
    # # Combine the counts from both counters
    # model_type_cnt = progress_model_type_cnt + constant_model_type_cnt


    model_type = ModelDF.objects.all().values_list('model',flat=True)

    model_type_cnt = Counter(model_type)
    # make trace
    pie_model_type = [{
        'type': "pie",
        'values': list(model_type_cnt.values()),
        'labels': list(model_type_cnt.keys()),
        'textinfo': "label+percent",
        'textposition': "outside",
        'automargin': 'true'
    }]
    # Histogram of number of gene for all perturbation datasets
    # hist_gene = pd.DataFrame(clustergene.objects.all().values_list('gene',flat=True))
    # geneall = pd.DataFrame(clustergene.objects.all().values_list('gene'), columns=['gene'])
    #保存文件为pkl
    # trend_path = os.getcwd() + f'\\static\\data\\studyid_model_gene_counts.txt'
    # df = pd.read_csv(trend_path, sep=' ', header=None, names=['studyid', 'model', 'gene_count'])
    # pkl_file = os.getcwd() + '\\static\\data\\gene_count.pkl'
    # with open(pkl_file, 'wb') as f:
    #     pickle.dump(df, f)
    # pkl_file = os.path.join(os.getcwd(), 'timecourse_app', 'static', 'data', 'gene_count.pkl')
    basedir = settings.STATICFILES_DIRS[0] + 'data/'
    pkl_file = basedir + 'gene_count.pkl'
    # pkl_file = os.getcwd() + '/timecourse_app_static/data/gene_count.pkl'
    df = pd.read_pickle(pkl_file)
    total_counts = df.groupby('studyid')['n.gene'].sum()
    total_counts_list = list(total_counts)
    # bin_size = int(np.ceil((max(total_counts_list) - min(total_counts_list)) / np.sqrt(len(total_counts_list))))
    # bins = np.arange(min(total_counts_list), max(total_counts_list) + bin_size, bin_size)
    # s1 = pd.cut(total_counts_list, bins=bins, include_lowest=True)
    # gene_count_freq = list(s1.value_counts().values)


    bins = np.arange(0, 3001, 50)
    s1 = pd.cut(total_counts_list, bins=bins, include_lowest=True)
    # s1 = pd.cut(total_counts_list, bins=[0, 100, 200, 300, 400, 500, 800, 1000, 1500, 2000, 3000, 4000, 10000])
    gene_count_freq = list(s1.value_counts().values)
    # gene_count_freq_list = [
    #     {
    #         'x': ['0~100', '101~200', '201~300', '301~400', '401~500', '501~800', '801~1000', '1001~1500', '1501~2000', '2001~3000', '3001~4000', '4001~'],
    #         'y': gene_count_freq,
    #         'type': 'histogram'
    #     }
    # ]

    # histogram_data = []
    # for studyid, model_counts in total_counts.items():
    #     total_count = total_counts[studyid]  # 获取该 studyid 的总数
    #     # 添加总数的柱子
    #     bar_data_total = {
    #         'x': [studyid],
    #         'y': [total_count],
    #         'name': 'Total',
    #         'type': 'bar'
    #     }
    #     histogram_data.append(bar_data_total)

    return render(request, 'statistic.html', {
        'total_study_cnt': total_study_cnt,
        'barplot_organism_strategy': barplot_organism_strategy,
        'pie_data_type': pie_data_type,
        'pie_model_type': pie_model_type,
        # 'hist_gene': list(hist_gene['count']),
        # 'histogram_gene': histogram_data,
        'gene_count_freq': gene_count_freq,
    })


def stat_freq():
    strategy_organism = pd.DataFrame(Sampleinfo.objects.all().values_list('strategy', 'organism'),
                                     columns=['strategy', 'organism']).value_counts()
    strategy_organism['strcat'] = strategy_organism['strategy'] + ';' + strategy_organism['organism']
    strategy_organism['strcat'].value_counts()


def dict2list(in_dict):
    k_list=[]
    v_list=[]
    [(k_list.append(k),v_list.append(v)) for k,v in in_dict.items()]
    return k_list,v_list