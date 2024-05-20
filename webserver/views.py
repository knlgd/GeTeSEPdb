import math
import pickle

from django.shortcuts import render
from django.http import HttpResponse
from django.db.models import Q, F, Case, FloatField
from django.core.paginator import Paginator
# from GENE_Database.settings import STATICFILES_DIRS
from djangoProject.settings import STATIC_ROOT
from timecourse_app.models import Sampleinfo, EvaIndicator, ModelDF, GO, KEGG, \
    clustergene  # CTMarkerGene,  KEGGTerm, GOTerm, GSEATerm,  GeneExpVal, DiffGenes
# from gene_edit.models import CellRatio, MarkerGene, SCDiffGenes # ScatterInfo,
# from gene_edit.views.analysis_views import symbol_ensembl_to_id
# , Volcano, AllSortGenes
import numpy as np
import time
import pandas as pd
import json, os, re
from django.conf import settings
import sys
from collections import Counter
# import gseapy as gp
# from gseapy.plot import gseaplot
# from matplotlib_venn import venn2_unweighted, venn2
# from sklearn import preprocessing
# import pickle

# Create your views here.
def hello_world(request):
    return HttpResponse("Hello World!")

def search(request):
    return render(request, 'search.html')

def result(request):
    key = request.GET['query_content']
    search_method = request.GET['search_method']
    result_information = {}
    if key == '':
        return render(request, 'ERROR.html')
    elif search_method == 'gene_search':
        genedf = clustergene.objects.filter(Q(gene__icontains=key) ).values_list(
            'studyid', 'gene', 'cluster').order_by('studyid')
        select_genedf = pd.DataFrame(genedf,columns=['studyid', 'gene', 'cluster']).fillna('Null')
        # select_genedf = pd.DataFrame(clustergene.objects.filter(gene=key).values_list('studyid', 'gene', 'cluster'),
        #                              columns=['studyid', 'gene', 'cluster'])

        if select_genedf['gene'].count() > 1000:
            filtered_studyids = select_genedf['studyid'][0:1000].unique()
        else:
            filtered_studyids = select_genedf['studyid'].unique()
        select_sampleinfo = pd.DataFrame(
            Sampleinfo.objects.filter(studyid__in=filtered_studyids).values_list('studyid', 'tissue_cell_type', 'organism', 'strategy', 'experiment_type'),
            columns=['studyid', 'tissue', 'organism', 'strategy', 'experiment_type'])

        select_modeldf = pd.DataFrame(ModelDF.objects.filter(studyid__in=filtered_studyids).values_list('studyid', 'cluster', 'model'),
                                      columns=['studyid', 'cluster', 'model'])

        cluster = 'cluster0'
        cluster0_data = select_genedf.loc[select_genedf['cluster'] == cluster]
        leftcluster_data = select_genedf.loc[select_genedf['cluster'] != cluster]

        new_modeldf = pd.DataFrame(columns=select_modeldf.columns)
        cluster0_studyids = cluster0_data['studyid']
        cluster0_gene = cluster0_data['gene']

        new_modeldf['studyid'] = cluster0_studyids
        new_modeldf['gene'] = cluster0_gene
        new_modeldf['cluster'] = cluster
        new_modeldf['model'] = 'Constant'

        merged_df = leftcluster_data.merge(select_modeldf, on=['studyid', 'cluster'], how='inner')
        result_df = merged_df[['studyid', 'cluster', 'model', 'gene']]
        merged_df_gene = pd.concat([new_modeldf, result_df], ignore_index=True)

        temp_results = select_sampleinfo.merge(merged_df_gene, on='studyid', how='left')
        temp_results = temp_results.replace({'Logis': 'Logistic', 'Sin': 'Trigonometric'})
        temp_results = temp_results.dropna(axis=0)
        temp_results = temp_results.drop_duplicates(keep=False)

        count_gene = temp_results['gene'].count()
        if count_gene > 1000:
            result_information = temp_results[(temp_results['gene'].str.contains(key))][0:1000]
        else:
            result_information = temp_results[(temp_results['gene'].str.contains(key))]

        result_information = result_information.fillna('')
        result_information = result_information.sort_values('model')

    else:
        result_information = Sampleinfo.objects.filter(
            Q(studyid__exact=key) | Q(accession__icontains=key) | Q(experiment_type__icontains=key) | Q(related_factor__icontains=key) |
            Q(drug__icontains=key) | Q(tissue_cell_type__icontains=key) | Q(organism__icontains=key)).values_list('studyid', 'accession', 'time', 'unit', 'related_factor', 'experiment_type',
                                                                                                                         'drug', 'tissue_cell_type', 'organism').order_by('studyid')
        result_information = pd.DataFrame(result_information,
                                          columns=['studyid', 'accession', 'time', 'unit', 'related_factor', 'experiment_type',
                                                   'drug', 'tissue_cell_type', 'organism']).fillna('Null')

    result_information = result_information.to_dict('records')
        # aa = pd.DataFrame(
        #     Sampleinfo.objects.filter(accession=key).values_list('studyid', 'accession', 'time','unit','related_factor', 'experiment_type',
        #                                            'drug', 'tissue_cell_type','organism'),
        #     columns=['studyid', 'accession', 'time','unit','related_factor', 'experiment_type',
        #                                            'drug', 'tissue_cell_type','organism'])

    return render(request, 'result.html', {
        'result_information': result_information,
        'query_content': key,
        'search_method': search_method,
        # 'study_search_str': study_search_str,
        # 'studyid_list': json.dumps({'studyid_list': studyid_list}),
        # 'show_result_analysis': show_result_analysis,
    })


def detail(request, studyid):
    '''''
    Basic Information
    '''''
    # studyid='TCD0101'

    # start_time = time.time()
    curr_study = Sampleinfo.objects.filter(studyid=studyid).values() \
        .annotate(organism=F('organism')).annotate(data_source=F('data_source')) \
        .annotate(library_strategy=F('strategy'))[0]

    curr_time = curr_study['time']
    curr_time_list = list(set(curr_time.split(';')))
    curr_time_list = sorted(curr_time_list, key=lambda x: (x.isdigit(), x))
    curr_study['time'] = ';'.join(curr_time_list)

    # curr_time = curr_study['time']
    # curr_time_list = list(set(curr_time.split(';')))
    # curr_time_list = [int(x) for x in curr_time_list]
    # curr_time_list.sort()
    # curr_time_list = [str(x) for x in curr_time_list]
    # curr_time_change = ';'.join(curr_time_list)
    # curr_study['time'] = curr_time_change

    # hide sample too long
    if len(curr_study['sample']) > 200:
        curr_study['sample'] = '<details><summary>Click to open</summary><p style="margin: 0;">' + curr_study['sample'] + '</p></details>'

    # pkl_file = os.path.join(os.getcwd(),'timecourse_app','static', 'data', 'result_pkl')
    basedir = settings.STATICFILES_DIRS[0] + 'data/'
    pkl_file = basedir + 'result_pkl/'
    # pkl_file = os.getcwd() + '/timecourse_app_static/data/result_pkl/'
    # pkl_file = pkl_file.replace("\\", "/")
    # Model results
    EvaIndicator_path = pkl_file + studyid + '_EvaIndicators.pkl'
    EvaIndicator_path = EvaIndicator_path.replace("\\", "/")
    cluster = 'cluster0'
    genedata_path = pkl_file + studyid + '_BrowseGene.pkl'
    genedata_path = genedata_path.replace("\\", "/")
    if os.path.exists(genedata_path) and os.path.exists(EvaIndicator_path):
        genedata = pd.read_pickle(genedata_path)
        genedata['mean'] = genedata['mean'].apply(lambda x: round(x, 4))
        genedata['CV'] = genedata['CV'].apply(lambda x: round(x, 4))
        genedata['cluster'] = genedata['cluster'].str.replace('_', '')
        genedf = genedata[(genedata['cluster'] == cluster)]
        studyid=genedf['Studyid'].drop_duplicates()[0]
        c = list(genedf['mean'])
        if len(c) > 5:
            c = c[:5]
            c[-1] = str(c[-1]) + '...'

        e = list(genedf['CV'])
        if len(e) > 5:
            e = e[:5]
            e[-1] = str(e[-1]) + '...'
        EvaIndicator_de = {
            'model': 'Constant',
            'cluster': cluster,
            'rmse':'-',
            'rsq': '-',
            'adjrsq': '-',
            'aic': '-',
            'mean': c,
            'CV': e,
            'Studyid':studyid
        }
        parameters_strc = str(EvaIndicator_de['mean'])  # 转换parameters为字符串形式
        parameters_stre = str(EvaIndicator_de['CV'])
        EvaIndicator_dict = pd.DataFrame([[EvaIndicator_de['model'],EvaIndicator_de['cluster'],EvaIndicator_de['rmse'],EvaIndicator_de['rsq'],EvaIndicator_de['adjrsq'],EvaIndicator_de['aic'],parameters_strc, parameters_stre,EvaIndicator_de['Studyid']]],
            columns=['model','cluster', 'rmse','rsq','adjrsq','aic','mean', 'CV','studyid'])

        EvaIndicator_data = pd.read_pickle(EvaIndicator_path)
        # 删除"model"列中值为"exponential_decay"的行
        EvaIndicator_data = EvaIndicator_data[EvaIndicator_data['model'] != 'exponential_decay']
        EvaIndicator_data['cluster'] = EvaIndicator_data['cluster'].str.replace('_', '')
        EvaIndicator_data = EvaIndicator_data.rename(columns={'Studyid': 'studyid', 'aicc': 'aic'})
        EvaIndicator_data = EvaIndicator_data.replace({'Logis': 'Logistic', 'Sin': 'Trigonometric', 'exponential_growth': 'Exponential'})
        EvaIndicator_data['rmse'] = EvaIndicator_data['rmse'].apply(lambda x: round(x, 4))
        EvaIndicator_data['rsq'] = EvaIndicator_data['rsq'].apply(lambda x: round(x, 4))
        EvaIndicator_data['adjrsq'] = EvaIndicator_data['adjrsq'].apply(lambda x: round(x, 4))
        EvaIndicator_data['aic'] = EvaIndicator_data['aic'].apply(lambda x: round(x, 4))
        EvaIndicator_data['mean'] = "-"
        EvaIndicator_data['CV'] = "-"
        EvaIndicator_data.replace(-float('inf'), -0.0001, inplace=True)
        EvaIndicator_combine=pd.concat([EvaIndicator_dict, EvaIndicator_data])
        EvaIndicator_df = EvaIndicator_combine.to_dict('records')
    elif os.path.exists(genedata_path):
        genedata = pd.read_pickle(genedata_path)
        genedata['mean'] = genedata['mean'].apply(lambda x: round(x, 4))
        genedata['CV'] = genedata['CV'].apply(lambda x: round(x, 4))
        genedata['cluster'] = genedata['cluster'].str.replace('_', '')
        genedf = genedata[(genedata['cluster'] == cluster)]
        studyid = genedf['Studyid'].drop_duplicates()[0]
        c = list(genedf['mean'])
        if len(c) > 5:
            c = c[:5]
            c[-1] = str(c[-1]) + '...'

        e = list(genedf['CV'])
        if len(e) > 5:
            e = e[:5]
            e[-1] = str(e[-1]) + '...'
        EvaIndicator_de = {
            'model': 'Constant',
            'cluster': cluster,
            'rmse': '-',
            'rsq': '-',
            'adjrsq': '-',
            'aic': '-',
            'mean': c,
            'CV': e,
            'Studyid': studyid
        }
        parameters_strc = str(EvaIndicator_de['mean'])  # 转换parameters为字符串形式
        parameters_stre = str(EvaIndicator_de['CV'])
        EvaIndicator_combine = pd.DataFrame([[EvaIndicator_de['model'], EvaIndicator_de['cluster'],
                                              EvaIndicator_de['rmse'], EvaIndicator_de['rsq'], EvaIndicator_de['adjrsq'],EvaIndicator_de['aic'],
                                              parameters_strc, parameters_stre, EvaIndicator_de['Studyid']]],
                                            columns=['model', 'cluster', 'rmse', 'rsq', 'adjrsq', 'aic', 'mean', 'CV', 'studyid'])
        EvaIndicator_df = EvaIndicator_combine.to_dict('records')
    else:
        EvaIndicator_df = {}

    # Other cluster ModelDF_df
    ModelDF_path = pkl_file + studyid + '_ModelDF.pkl'
    ModelDF_path = ModelDF_path.replace("\\", "/")
    if os.path.exists(ModelDF_path):
        ModelDF_df = pd.read_pickle(ModelDF_path)
        ModelDF_df['row_number'] = ModelDF_df.groupby(['cluster', 'modle']).cumcount()
        # 仅保留row_number为0的行
        ModelDF_df = ModelDF_df[ModelDF_df['row_number'] == 0]
        # 删除row_number列
        ModelDF_df = ModelDF_df.drop(columns=['row_number'])
        ModelDF_df['cluster'] = ModelDF_df['cluster'].str.replace('_', '')
        ModelDF_df = ModelDF_df.drop('Studyid', axis=1)
        ModelDF_df['Parameters'] = ModelDF_df['Parameters'].apply(format_parameters)
        ModelDF_df = ModelDF_df.rename(columns={'Parameters': 'parameters','modle': 'model'})
        ModelDF_df = ModelDF_df.replace( {'Logis': 'Logistic', 'Sin': 'Trigonometric'})
    else:
        ModelDF_df = {}


    # Cluster0 ModelDF_df
    cluster = 'cluster0'
    genedata_path = pkl_file + studyid + '_BrowseGene.pkl'
    genedata_path = genedata_path.replace("\\", "/")
    heatmap_path = pkl_file + studyid + '_' + cluster + '_heatmapdata0409.pkl'
    heatmap_path = heatmap_path.replace("\\", "/")
    ModelDF_combine = {}
    if os.path.exists(genedata_path) and os.path.exists(ModelDF_path) and os.path.exists(heatmap_path):
        # genedata = pd.read_pickle(genedata_path)
        # genedata['mean'] = genedata['mean'].apply(lambda x: round(x, 4))
        # genedata['cluster'] = genedata['cluster'].str.replace('_', '')
        # genedf = genedata[(genedata['cluster'] == cluster)]
        # c=list(genedf['mean'])
        # mean_c = sum(c) / len(c)
        # formatted_mean_c = round(mean_c, 4)
        # print(formatted_mean_c)
        # if len(c) > 5:
        #     c = c[:5]
        #     c[-1] = str(c[-1]) + '...'
        cluster0_data = pd.read_pickle(heatmap_path)
        cluster0_data.set_index('Genes', inplace=True)
        gene_means = cluster0_data.mean(axis=1)
        # 计算所有基因均值的均值
        overall_mean = gene_means.mean()
        formatted_mean_c = round(overall_mean, 4)
        ModelDF_dc = {
            'cluster': cluster,
            'model': 'Constant',
            'parameters': formatted_mean_c,
            'exp': 'y=c'
        }
        # parameters_str = 'c= ' + str(ModelDF_dc['parameters'])  # 转换parameters为字符串形式
        parameters_str = formatted_mean_c
        ModelDF_dc_dict = pd.DataFrame(
            [[ModelDF_dc['cluster'], ModelDF_dc['model'], parameters_str, ModelDF_dc['exp']]],
            columns=['cluster', 'model', 'parameters', 'exp'])
        ModelDF_combine = pd.concat([ModelDF_dc_dict, ModelDF_df])
        cluster_list = ModelDF_combine['cluster'].unique().tolist()
        ModelDF_result = ModelDF_combine.to_dict('records')
    elif os.path.exists(heatmap_path):
        # genedata = pd.read_pickle(genedata_path)
        # genedata['mean'] = genedata['mean'].apply(lambda x: round(x, 4))
        # genedata['cluster'] = genedata['cluster'].str.replace('_', '')
        # genedf = genedata[(genedata['cluster'] == cluster)]
        # c = list(genedf['mean'])
        # mean_c = sum(c) / len(c)
        # formatted_mean_c = round(mean_c, 4)
        # print(formatted_mean_c)
        # if len(c) > 5:
        #     c = c[:5]
        #     c[-1] = str(c[-1]) + '...'

        # cluster_list = ModelDF_combine['cluster'].unique().tolist()
        cluster0_data = pd.read_pickle(heatmap_path)
        cluster0_data.set_index('Genes', inplace=True)
        gene_means = cluster0_data.mean(axis=1)
        # 计算所有基因均值的均值
        overall_mean = gene_means.mean()
        formatted_mean_c = round(overall_mean, 4)
        ModelDF_dc = {
            'cluster': cluster,
            'model': 'Constant',
            'parameters': formatted_mean_c,
            'exp': 'y=c'
        }
        # parameters_str = 'c= ' + str(ModelDF_dc['parameters'])  # 转换parameters为字符串形式
        parameters_str = formatted_mean_c
        ModelDF_combine = pd.DataFrame(
            [[ModelDF_dc['cluster'], ModelDF_dc['model'], parameters_str, ModelDF_dc['exp']]],
            columns=['cluster', 'model', 'parameters', 'exp'])
        ModelDF_result = ModelDF_combine.to_dict('records')
    else:
        ModelDF_result={}

    heatmap_data = {}
    trend_data = {}
    tread_data_line = {}
    x_Axis=[]
    if not os.path.exists(ModelDF_path):
        cluster_list = ['cluster0']

    cluster_model_list = [f"{cluster}-{model}" for cluster, model in
                          zip(ModelDF_combine['cluster'], ModelDF_combine['model'])]
    print(cluster_model_list)
    for cluster_model in cluster_model_list:
        cluster= cluster_model.split('-')[0]
        model=cluster_model.split('-')[1]
        basedir = settings.STATICFILES_DIRS[0] + 'data/'
        pkl_file = basedir + 'result_pkl/'
        heatmap_path = pkl_file + studyid + '_' + cluster + '_heatmapdata0409.pkl'
        heatmap_path = heatmap_path.replace("\\", "/")
        if os.path.exists(heatmap_path):
            data = pd.read_pickle(heatmap_path)  # 不使用列索引
            data = data.dropna()
            data = data.set_index(data.columns[0], drop=True)
            heatmap_data[cluster_model] = {
                'x': list(data.columns),
                'y': list(data.index),
                'z': data.values.tolist()
            }
        else:
            heatmap_data[cluster_model] = None

        if os.path.exists(heatmap_path):
            trend_data_all = pd.read_pickle(heatmap_path)
            trend_data_all = trend_data_all.dropna()
            trend_data_all = trend_data_all.set_index(trend_data_all.columns[0], drop=True)
            time = list(trend_data_all.columns)[0:]  # 将获取时间信息的代码提到这里
            trend_data[cluster_model] = {
                'x': time,
                'y': {}
            }
            if not trend_data_all.empty:
                for gene, values in trend_data_all.iterrows():
                    trend_data[cluster_model]['y'][gene] = list(values)[0:]

            Expressions = ModelDF_combine[(ModelDF_combine['cluster'] == cluster) & (ModelDF_combine['model'] == model)]['exp'].iloc[0]
            Parameters = ModelDF_combine[(ModelDF_combine['cluster'] == cluster) & (ModelDF_combine['model'] == model)]['parameters'].iloc[0]

            if 'exp' in Expressions:
                Expressions = Expressions.replace('exp', 'math.exp')
            if 'log' in Expressions:
                Expressions = Expressions.replace('log', 'math.log')
            if 'pow' in Expressions:
                Expressions = Expressions.replace('power', 'math.pow')
            if 'sin' in Expressions:
                Expressions = Expressions.replace('sin', 'math.sin')
                Expressions = Expressions.replace('pi', str(math.pi))

            Expressions = Expressions.replace(' ', '')
            Expressions = Expressions.replace('^', '**')
            # 先将参数赋值
            x_list = np.arange(1, len(time) + 0.001, 0.001)
            y_list = list()

            def get_y(x, Parameters_list):
                for para in Parameters_list:
                    exec(para)
                exec(Expressions)
                local_y = locals()['y']
                return local_y

            if cluster == 'cluster0':
                Parameters_list = [Parameters]
                y_value = float(Parameters)  # 假设 Parameters 包含固定的数值
                y_list = [y_value] * len(x_list)
                tread_data_line[cluster_model] = {
                    'x': x_list.tolist(),
                    'y': y_list
                }
                x_Axis = list(range(1, len(time) + 1))
            else:
                Parameters_list = Parameters.split('; ')
                for k in range(0, len(x_list)):
                    x = x_list[k]
                    local_y = get_y(x, Parameters_list)
                    y_list.append(local_y)

                tread_data_line[cluster_model] = {
                    'x': x_list.tolist(),
                    'y': y_list
                }
                x_Axis = list(range(1, len(time) + 1))
        else:
            trend_data[cluster_model] = {
                'x': [],
                'y': {}
            }
            tread_data_line[cluster_model] = {
                'x': [],
                'y': {}
            }

    # print(time.time() - start_time)
    cluster = cluster_list[0]

    GO_path = pkl_file + studyid + '_GO0409.pkl'
    GO_path = GO_path.replace("\\", "/")
    if os.path.exists(GO_path):
        GO_data = pd.read_pickle(GO_path)
        # GO_data['cluster'] = GO_data['cluster'].str.replace('_', '')
        GO_data = GO_data.rename(columns={'Studyid': 'studyid'})
        # GO_data= pd.DataFrame(GO.objects.filter(studyid=studyid).filter(cluster=cluster).values_list('GO_ID', 'Description', 'ONT',
        #                                                                                   'GeneRatio', 'BgRatio', 'pvalue',
        #                                                                                   'p_adjust', 'qvalue', 'geneID',
        #                                                                                   'Count'),columns=['GO_ID', 'Description', 'ONT',
        #                                                                                   'GeneRatio', 'BgRatio', 'pvalue',
        #                                                                                   'p_adjust', 'qvalue', 'geneID',
        #                                                                                   'Count'])

        GO_data['p.adjust'] = GO_data['p.adjust'].apply(lambda x: 0.00001 if x == 0 else x)
        GO_data = GO_data[GO_data['cluster'] == cluster]
    else:
        GO_data= pd.DataFrame(columns=['GO_ID', 'Description', 'ONT','GeneRatio', 'BgRatio', 'pvalue',
                                       'p_adjust', 'qvalue', 'geneID','Count','studyid','cluster'])

    KEGG_path = pkl_file + studyid + '_KEGG0409.pkl'
    KEGG_path = KEGG_path.replace("\\", "/")
    if os.path.exists(KEGG_path):
        KEGG_data = pd.read_pickle(KEGG_path)
        # GO_data['cluster'] = GO_data['cluster'].str.replace('_', '')
        KEGG_data = KEGG_data.rename(columns={'Studyid': 'studyid'})
    # KEGG_data = pd.DataFrame(KEGG.objects.filter(studyid=studyid).filter(cluster=cluster).values_list('KEGG_ID', 'Description', 'GeneRatio', 'BgRatio',
    #                                         'pvalue', 'p_adjust', 'qvalue', 'geneID', 'Count'),columns=['KEGG_ID', 'Description', 'GeneRatio', 'BgRatio',
    #                                         'pvalue', 'p_adjust', 'qvalue', 'geneID', 'Count'])
        KEGG_data['p.adjust'] = KEGG_data['p.adjust'].apply(lambda x: 0.00001 if x == 0 else x)
        KEGG_data = KEGG_data[KEGG_data['cluster'] == cluster]
    else:
        KEGG_data= pd.DataFrame(columns=['KEGG_ID', 'Description', 'GeneRatio', 'BgRatio',
                                            'pvalue', 'p_adjust', 'qvalue', 'geneID', 'Count','studyid','cluster'])

    '''''
    GO bubble show
    '''''
    if GO_data.empty:
        enrich_obj = None
    else:
        enrich_result = GO_data[GO_data['ONT'] == 'BP']
        if enrich_result.shape[0] != 0:
            enrich_obj = bubble_obj_fun(enrich_result)
        else:
            enrich_obj = None

        # if exist go file
    if GO_data.empty:
        go_exist_or_not = 'Notexist'
    else:
        go_exist_or_not = 'Exist'
    # if exist kegg file
    if KEGG_data.empty:
        kegg_exist_or_not = 'Notexist'
    else:
        kegg_exist_or_not = 'Exist'

    print(go_exist_or_not)
    print(kegg_exist_or_not)


    return render(request, 'detail.html', {
        'curr_study': curr_study,
        'heatmap_data': heatmap_data,
        'cluster_model_list': cluster_model_list,
        'cluster_list': cluster_list,
        'trend_data': trend_data,
        'tread_data_line': tread_data_line,
        'x_Axis': x_Axis,
        'EvaIndicator_df': EvaIndicator_df,
        'ModelDF_df': ModelDF_result,
        'GO_data': GO_data,
        'KEGG_data': KEGG_data,
        'enrich_obj': enrich_obj,
        #'md5': query_md5,
        #'mode': mode,
        'go_exist_or_not': go_exist_or_not,
        'kegg_exist_or_not': kegg_exist_or_not,
    })

def get_bubble(request):
    studyid = request.POST.get('studyid')
    ont = request.POST.get('ont')
    cluster = request.POST.get('cluster')
    # studyid = 'TCD1711'
    # ont = 'kegg'
    # cluster = 'cluster1'
    # print(studyid)
    # print(ont)
    # print(cluster)
    basedir = settings.STATICFILES_DIRS[0] + 'data/'
    pkl_file = basedir + 'result_pkl/'
    # pkl_file = os.getcwd() + '/timecourse_app_static/data/result_pkl/'
    pkl_file = pkl_file.replace("\\", "/")


    GO_path = pkl_file + studyid + '_GO0409.pkl'
    GO_path = GO_path.replace("\\", "/")
    if os.path.exists(GO_path):
        GO_data = pd.read_pickle(GO_path)
        # GO_data['cluster'] = GO_data['cluster'].str.replace('_', '')
        GO_data = GO_data.rename(columns={'Studyid': 'studyid'})
        # GO_data= pd.DataFrame(GO.objects.filter(studyid=studyid).filter(cluster=cluster).values_list('GO_ID', 'Description', 'ONT',
        #                                                                                   'GeneRatio', 'BgRatio', 'pvalue',
        #                                                                                   'p_adjust', 'qvalue', 'geneID',
        #                                                                                   'Count'),columns=['GO_ID', 'Description', 'ONT',
        #                                                                                   'GeneRatio', 'BgRatio', 'pvalue',
        #                                                                                   'p_adjust', 'qvalue', 'geneID',
        #                                                                                   'Count'])

        GO_data['p.adjust'] = GO_data['p.adjust'].apply(lambda x: 0.00001 if x == 0 else x)
        GO_data = GO_data[GO_data['cluster'] == cluster]
    else:
        GO_data = pd.DataFrame(columns=['GO_ID', 'Description', 'ONT', 'GeneRatio', 'BgRatio', 'pvalue',
                                        'p_adjust', 'qvalue', 'geneID', 'Count', 'studyid', 'cluster'])

    KEGG_path = pkl_file + studyid + '_KEGG0409.pkl'
    KEGG_path = KEGG_path.replace("\\", "/")
    if os.path.exists(KEGG_path):
        KEGG_data = pd.read_pickle(KEGG_path)
        # GO_data['cluster'] = GO_data['cluster'].str.replace('_', '')
        KEGG_data = KEGG_data.rename(columns={'Studyid': 'studyid'})
        # KEGG_data = pd.DataFrame(KEGG.objects.filter(studyid=studyid).filter(cluster=cluster).values_list('KEGG_ID', 'Description', 'GeneRatio', 'BgRatio',
        #                                         'pvalue', 'p_adjust', 'qvalue', 'geneID', 'Count'),columns=['KEGG_ID', 'Description', 'GeneRatio', 'BgRatio',
        #                                         'pvalue', 'p_adjust', 'qvalue', 'geneID', 'Count'])
        KEGG_data['p.adjust'] = KEGG_data['p.adjust'].apply(lambda x: 0.00001 if x == 0 else x)
        KEGG_data = KEGG_data[KEGG_data['cluster'] == cluster]
    else:
        KEGG_data = pd.DataFrame(columns=['KEGG_ID', 'Description', 'GeneRatio', 'BgRatio',
                                          'pvalue', 'p_adjust', 'qvalue', 'geneID', 'Count', 'studyid', 'cluster'])
    if ont != 'kegg':
            # enrich_result = pd.DataFrame(
            #     GOTerm.objects.filter(studyid__studyid=studyid).values('description', 'ont', 'pvalue', 'count'),
            #     columns=['description', 'ont', 'pvalue', 'count'])
        if ont == 'undefined':
            ont = 'bp'

        if GO_data.empty:
            enrich_obj = None
        else:
                # enrich_result.to_csv("D:/Project/Database/gene_edit/paper/supplemental_files/S4.csv")
            enrich_result = GO_data[GO_data['ONT'] == ont.upper()]
            enrich_obj = bubble_obj_fun(enrich_result)

    else:
            # enrich_result = pd.DataFrame(
            #     KEGGTerm.objects.filter(studyid__studyid=studyid).values('description', 'pvalue', 'count'),
            #     columns=['description', 'pvalue', 'count'])
        if KEGG_data.empty:
            enrich_obj = None
        else:
            enrich_obj = bubble_obj_fun(KEGG_data)
    return HttpResponse(json.dumps(enrich_obj))

def bubble_obj_fun(enrich_result):
    enrich_result['nlog10pvalue'] = -np.log10(enrich_result['p.adjust'])
    enrich_result_top10 = enrich_result.iloc[0:10, ].sort_values(by='nlog10pvalue')
    x = enrich_result_top10['nlog10pvalue'].tolist()
    y = enrich_result_top10['Description'].tolist()
    count = enrich_result_top10['Count'].tolist()
    if max(enrich_result_top10['Count']) > 30:
        size = (enrich_result_top10['Count'] * 30 / max(enrich_result_top10['Count'])).tolist()
    elif max(enrich_result_top10['Count']) < 5:
        size = [(_ + 10) for _ in count]
    else:
        size = count
    enrich_obj = [{
        'mode': 'markers',
        # 'name': 'industry1',
        'type': 'scatter',
        'x': x,
        'y': y,
        'marker': {
            'line': {
                'width': 1.3
            },
            'color': 'rgba(55, 128, 191, 1.0)',
            'symbol': 'dot',
            'opacity': 0.8,
            'size': size
        },
        'text': ['size: ' + str(c) for c in count],
        'textfont': {
            'color': '#4D5663'
        }
    }]
    return (enrich_obj)

def get_lineplot(request):
    studyid = request.POST.get('studyid')
    cluster = request.POST.get('cluster')

    basedir = settings.STATICFILES_DIRS[0] + 'data/'
    pkl_file = basedir + 'result_pkl/'
    trend_path = pkl_file + studyid + '_' + cluster +'_heatmapdata0409.pkl'
    trend_path = trend_path.replace("\\", "/")
    if os.path.exists(trend_path):
        trend_data_all = pd.read_pickle(trend_path)  # 不使用列索引
        trend_data_all = trend_data_all.set_index(trend_data_all.columns[0], drop=True)
        time = list(trend_data_all.columns)[1:]
        trend_data = {}
        for gene, values in trend_data_all.iterrows():
            trend_data[gene] = list(values)[1:]
    else:
        trend_data = None

    geneData = trend_data
    traceData_list = []
    for gene in geneData:
        trace = {
            'x': time,
            'y': geneData[gene],
            'name': gene,
            'type': 'scatter'
        }
        traceData_list.append(trace)

    return HttpResponse(json.dumps(traceData_list))

def process_numbers(text):
    numbers = re.findall(r'[-+]?\d*\.\d+|\d+', text)  # 提取数字
    rounded_numbers = [round(float(num), 4) for num in numbers]  # 保留四位小数
    return rounded_numbers

# 分割参数字符串并保留4位小数
def format_parameters(parameters_str):
    params = parameters_str.split('; ')
    formatted_params = []
    for param in params:
        key, value = param.split('=')
        formatted_value = f'{float(value):.4f}'
        formatted_params.append(f'{key} = {formatted_value}')
    return '; '.join(formatted_params)

"""
公共方法
"""
"""
说明：页面底部页码栏标号方法
输入：request获取的page参数,要分页的元素列表element_list
输出：当前页元素，上一页页码，下一页页码，页码栏的开始和结束标号
"""
def divide_page(page, element_list):
    if page:
        page = int(page)
    else:
        page = 1
    paginator = Paginator(element_list, 10)
    page_nums = paginator.num_pages
    if page > page_nums:
        page = page_nums
    elif page < 1:
        page = 1
    cur_element_list = paginator.page(page)
    previous_page, next_page, begin_page, end_page = page_label(page, page_nums, cur_element_list)
    curr_page = page
    return cur_element_list, curr_page, previous_page, next_page, begin_page, end_page


"""
说明：页面底部页码栏标号方法
输入：当前页，总页数，当前页内容列表
输出：上一页页码，下一页页码，页码栏的开始和结束标号
"""
def page_label(page, page_nums, curr_page_list):
    if curr_page_list.has_next():
        next_page = page + 1
    else:
        next_page = page
    if curr_page_list.has_previous():
        previous_page = page - 1
    else:
        previous_page = page
    # 开始页和结束页
    if page - 3 <= 0:
        begin_page = 1
    else:
        begin_page = page - 3
    if page_nums <= 10:
        end_page = page_nums
    elif page + 6 >= page_nums:
        begin_page = page_nums - 9
        end_page = page_nums
    else:
        end_page = begin_page + 9
    return previous_page, next_page, begin_page, end_page


# 500报错页面
def server_error(request):
    return render(request, "500.html")

# 404 找不到页面
def page_not_found(request):
    return render(request, "404.html")

# 400错误页面
def bad_request(request):
    return render(request, "400.html")
