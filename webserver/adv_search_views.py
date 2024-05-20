import os
import random
from django.shortcuts import render
from django.core.paginator import Paginator ##Django分页组件
from django.core import serializers
from django.db.models import Q
from django.http import HttpResponse, HttpResponseBadRequest
from django.conf import settings
import pandas as pd
import time
import numpy as np
import math, json
# from sklearn.preprocessing import StandardScaler
from timecourse_app.models import Sampleinfo, clustergene, ModelDF,geneinfo

# def get_adv_search_page(request):
#     allsampleinfo = pd.DataFrame(
#         Sampleinfo.objects.all().values_list('studyid', 'accession', 'time', 'tissue_cell_type', 'organism', 'strategy',
#                                              'experiment_type'),
#         columns=['studyid', 'accession', 'time', 'tissue', 'organism', 'strategy', 'experiment_type']).drop_duplicates()
#
#     modeldf = pd.DataFrame(ModelDF.objects.all().values_list('studyid', 'cluster', 'model', 'parameters', 'exp'),
#                            columns=['studyid', 'cluster', 'model', 'parameters', 'exp'])
#
#     # #     #organism
#     # organism = allsampleinfo[allsampleinfo['organism'].notna()]
#     # #all organism
#     # all_organism = list(set(list(organism['tissue'])))
#     # # Tissue
#     # all_tissues = allsampleinfo[allsampleinfo['tissue'].notna()]['tissue'].unique()
#     # all_tissues.sort()
#     # all_tissue_option = "".join(['<option>' + t + '</option>' for t in all_tissues])
#     # # Model Type
#     # all_model_types = modeldf['model'].unique()
#     # all_model_types = np.append(all_model_types, 'Constant')
#     # all_model_types.sort()
#     # all_modeltype_option = "".join(['<option>' + t + '</option>' for t in all_model_types])
#     # # Strategy
#     # all_strategies = allsampleinfo[allsampleinfo['strategy'].notna()]['strategy'].unique()
#     # all_strategies.sort()
#     # all_strategy_option = "".join(['<option>' + t + '</option>' for t in all_strategies])
#     # # Experiment Type
#     # all_experiment_types = allsampleinfo[allsampleinfo['experiment_type'].notna()]['experiment_type'].unique()
#     # all_experiment_types.sort()
#     # all_experiment_option = "".join(['<option>' + t + '</option>' for t in all_experiment_types])
#     #organism
#     organism = allsampleinfo[allsampleinfo['organism'].notna()]
#     # all organism
#     all_organism = list(set(list(organism['tissue'])))
#     # tissue
#     tissue = allsampleinfo[allsampleinfo['tissue'].notna()][['tissue', 'organism']]
#     # all tissue
#     all_tissue = list(set(list(tissue['tissue'])))
#     # single tissue
#     all_tissue.sort()
#     all_tissue_option = "".join(['<option>' + t + '</option>' for t in all_tissue])
#     # single organism tissue
#     organism_dict_tissue = {}
#     for o in tissue['organism'].unique():
#         single_tissue = list(set(list(tissue[tissue['organism'] == o]['tissue'])))
#         single_tissue.sort()
#         organism_dict_tissue[o] = "".join(['<option>' + t + '</option>' for t in single_tissue])
#
#     # # Model Type
#     all_model_types = modeldf['model'].unique()
#     all_model_types = np.append(all_model_types, 'Constant')
#     all_model_types.sort()
#     all_modeltype_option = "".join(['<option>' + t + '</option>' for t in all_model_types])
#     # # single organism modeltype
#     # organism_dict_mt = {}
#     # for o in modeltype['organism'].unique():
#     #     single_mt = list(set(list(modeltype[modeltype['organism'] == o]['model'])))
#     #     single_mt.sort()
#     #     organism_dict_mt[o] = "".join(['<option>' + t + '</option>' for t in single_mt])
#
#     # library strategy
#     strategy = allsampleinfo[allsampleinfo['strategy'].notna()][['strategy', 'organism']]
#     # all strategy
#     all_strategy = list(set(list(strategy['strategy'])))
#     # single strategy
#     all_strategy.sort()
#     all_strategy_option = "".join(['<option>' + t + '</option>' for t in all_strategy])
#     # single organism strategy
#     organism_dict_ls = {}
#     for o in strategy['organism'].unique():
#         single_ls = list(set(list(strategy[strategy['organism'] == o]['strategy'])))
#         single_ls.sort()
#         organism_dict_ls[o] = "".join(['<option>' + t + '</option>' for t in single_ls])
#
#     # experiment strategy
#     experiment = allsampleinfo[allsampleinfo['experiment_type'].notna()][['experiment_type', 'organism']]
#     # all experiment
#     all_experiment = list(set(list(experiment['experiment_type'])))
#     # single experiment
#     all_experiment.sort()
#     all_experiment_option = "".join(['<option>' + t + '</option>' for t in all_experiment])
#     # single organism strategy
#     organism_dict_ex = {}
#     for o in experiment['organism'].unique():
#         single_ex = list(set(list(experiment[experiment['organism'] == o]['experiment_type'])))
#         single_ex.sort()
#         organism_dict_ex[o] = "".join(['<option>' + t + '</option>' for t in single_ex])
#
#     return render(request, 'adv_search.html', {
#         'all_tissue_option': all_tissue_option,
#         'organism_dict_tissue': organism_dict_tissue,
#         'all_modeltype_option': all_modeltype_option,
#         # 'organism_dict_mt':organism_dict_mt,
#         'all_strategy_option': all_strategy_option,
#         'organism_dict_ls': organism_dict_ls,
#         'all_experiment_option': all_experiment_option,
#         'organism_dict_ex': organism_dict_ex
#     })
# def get_adv_search_page(request):
#     # filter match gene id
#     allsampleinfo = pd.DataFrame(
#         Sampleinfo.objects.all().values_list('studyid', 'accession', 'time', 'tissue_cell_type', 'organism', 'strategy','experiment_type'),
#         columns=['studyid', 'accession', 'time', 'tissue', 'organism', 'strategy','experiment_type']).drop_duplicates()
#
#     genedf = pd.DataFrame(clustergene.objects.all().values_list('studyid', 'gene', 'cluster', 'mean'),
#                           columns=['studyid', 'gene', 'cluster', 'parameters'])
#     modeldf = pd.DataFrame(ModelDF.objects.all().values_list('studyid', 'cluster', 'model', 'parameters', 'exp'),
#                            columns=['studyid', 'cluster', 'model', 'parameters', 'exp'])
#     cluster = 'cluster0'
#     cluster0_data = genedf.loc[genedf['cluster'] == cluster]
#     leftcluster_data=genedf.loc[genedf['cluster'] != cluster]
#     leftcluster_data=leftcluster_data.drop(['parameters'], axis=1)
#     lefttemp_results = modeldf.merge(leftcluster_data, on=["studyid", "cluster"], how='left')
#
#     # 构建新的DataFrame保存cluster0数据
#     new_modeldf = pd.DataFrame(columns=modeldf.columns)
#
#     for index, row in cluster0_data.iterrows():
#         c = row['parameters']
#         studyid = row['studyid']
#         ModelDF_dd = {
#             'studyid': studyid,
#             'cluster': cluster,
#             'model': 'constant',
#             'parameters': c,
#             'exp': 'y=c'
#         }
#         new_modeldf = new_modeldf.append(ModelDF_dd, ignore_index=True)
#
#     cluster0_results = pd.merge(new_modeldf, cluster0_data, on=["studyid", "cluster", "parameters"], how="left")
#
#     merged_df = pd.concat([lefttemp_results, cluster0_results], ignore_index=True)
#     temp_results = allsampleinfo.merge(merged_df, on='studyid', how='left')
#
#
#     #organism
#     # all_organism = temp_results[temp_results['organism'].notna()]
#     # all organism
#     # all_organism = list(set(list(organism['tissue'])))
#     # tissue
#     tissue = temp_results[temp_results['tissue'].notna()][['tissue', 'organism']]
#     # all tissue
#     all_tissue = list(set(list(tissue['tissue'])))
#     # single tissue
#     all_tissue.sort()
#     all_tissue_option = "".join(['<option>' + t + '</option>' for t in all_tissue])
#     # single organism tissue
#     organism_dict_tissue = {}
#     for o in tissue['organism'].unique():
#         single_tissue = list(set(list(tissue[tissue['organism'] == o]['tissue'])))
#         single_tissue.sort()
#         organism_dict_tissue[o] = "".join(['<option>' + t + '</option>' for t in single_tissue])
#
#     # modeltype
#     modeltype = temp_results[temp_results['model'].notna()][['model', 'organism']]
#     # all modeltype
#     all_modeltype = list(set(list(modeltype['model'])))
#     # single modeltype
#     all_modeltype.sort()
#     all_modeltype_option = "".join(['<option>' + t + '</option>' for t in all_modeltype])
#     # single organism modeltype
#     organism_dict_mt = {}
#     for o in modeltype['organism'].unique():
#         single_mt = list(set(list(modeltype[modeltype['organism'] == o]['model'])))
#         single_mt.sort()
#         organism_dict_mt[o] = "".join(['<option>' + t + '</option>' for t in single_mt])
#
#     # library strategy
#     strategy = temp_results[temp_results['strategy'].notna()][['strategy', 'organism']]
#     # all strategy
#     all_strategy = list(set(list(strategy['strategy'])))
#     # single strategy
#     all_strategy.sort()
#     all_strategy_option = "".join(['<option>' + t + '</option>' for t in all_strategy])
#     # single organism strategy
#     organism_dict_ls = {}
#     for o in strategy['organism'].unique():
#         single_ls = list(set(list(strategy[strategy['organism'] == o]['strategy'])))
#         single_ls.sort()
#         organism_dict_ls[o] = "".join(['<option>' + t + '</option>' for t in single_ls])
#
#     # experiment strategy
#     experiment = temp_results[temp_results['experiment_type'].notna()][['experiment_type', 'organism']]
#     # all experiment
#     all_experiment = list(set(list(experiment['experiment_type'])))
#     # single experiment
#     all_experiment.sort()
#     all_experiment_option = "".join(['<option>' + t + '</option>' for t in all_experiment])
#     # single organism strategy
#     organism_dict_ex = {}
#     for o in experiment['organism'].unique():
#         single_ex = list(set(list(experiment[experiment['organism'] == o]['experiment_type'])))
#         single_ex.sort()
#         organism_dict_ex[o] = "".join(['<option>' + t + '</option>' for t in single_ex])
#
#     return render(request, 'adv_search.html', {
#         'all_tissue_option': all_tissue_option,
#         'organism_dict_tissue': organism_dict_tissue,
#         'all_modeltype_option': all_modeltype_option,
#         'organism_dict_mt':organism_dict_mt,
#         'all_strategy_option': all_strategy_option,
#         'organism_dict_ls': organism_dict_ls,
#         'all_experiment_option': all_experiment_option,
#         'organism_dict_ex': organism_dict_ex
#     })

# def get_datatype(request):
#     organism = request.POST.get('organism')
#     gene_name = request.POST.get('gene_name')
#     model_type = request.POST.get('model_type')
#     tissue = request.POST.get('tissue')
#     data_type = request.POST.get('data_type')
#
#     # organism='All Organisms'
#     # gene_name = ''
#     # model_type = 'All Model Types'
#     # tissue = 'All Tissues'
#     # data_type = 'All Data Types'
#     print(organism,gene_name,model_type,tissue,data_type)
#
#     if (organism == 'All Organisms' and gene_name == '' and model_type == 'All Model Types' and tissue == 'All Tissues/Cell' and data_type == 'All Data Types'):
#         organism='Mus_musculus'
#         model_type = 'Logis'
#         tissue = 'brain'
#         data_type = 'RNA-Seq'
#         # gene_name='ACP1'
#         #start_time = time.time()
#         allsampleinfo = pd.DataFrame(Sampleinfo.objects.filter(organism=organism,strategy=data_type,tissue_cell_type=tissue).values_list('studyid',
#                          'tissue_cell_type', 'organism', 'strategy', 'experiment_type'),columns=['studyid','tissue', 'organism', 'strategy', 'experiment_type'])
#
#         modeldf = pd.DataFrame(ModelDF.objects.filter(model=model_type).values_list('studyid', 'cluster', 'model'),
#                            columns=['studyid', 'cluster', 'model'])
#
#         filtered_studyids = allsampleinfo['studyid'].unique()
#
#         genedf = pd.DataFrame(clustergene.objects.filter(studyid__in=filtered_studyids).values_list('studyid', 'gene','cluster'),
#             columns=['studyid', 'gene', 'cluster'])
#
#         cluster = 'cluster0'
#         cluster0_data = genedf.loc[genedf['cluster'] == cluster]
#         leftcluster_data = genedf.loc[genedf['cluster'] != cluster]
#         # leftcluster_data = leftcluster_data.drop('parameters', axis=1)
#
#         # 创建新的 DataFrame 保存 cluster0 数据
#         new_modeldf = pd.DataFrame(columns=modeldf.columns)
#
#         # 提取 cluster0_data 的相关列
#         cluster0_studyids = cluster0_data['studyid']
#         # cluster0_parameters = cluster0_data['parameters']
#         cluster0_gene = cluster0_data['gene']
#         # 构建新的 ModelDF 数据
#         new_modeldf['studyid'] = cluster0_studyids
#         new_modeldf['gene'] = cluster0_gene
#         new_modeldf['cluster'] = cluster
#         new_modeldf['model'] = 'constant'
#         # new_modeldf['parameters'] = cluster0_parameters
#         # new_modeldf['exp'] = 'y=c'
#
#         # 将新的 ModelDF 数据与 不包含cluster0的ModelDF 进行合并
#         merged_df = leftcluster_data.merge(modeldf, on=['studyid', 'cluster'], how='inner')
#         result_df = merged_df[['studyid', 'cluster', 'model', 'gene']]
#         merged_df_gene = pd.concat([new_modeldf, result_df], ignore_index=True)
#         # part_results = select_genedf.merge(merged_df, on=['studyid', 'cluster'], how='left')
#
#         temp_results = allsampleinfo.merge(merged_df_gene, on='studyid', how='left')
#         # print(time.time() - start_time)
#         temp_results = temp_results.dropna(subset=['gene'])
#         result_datatype_json = temp_results.to_dict('records')
#     else:
#         select_sampleinfo = pd.DataFrame(Sampleinfo.objects.filter(strategy=data_type,organism=organism,tissue_cell_type=tissue).values_list('studyid',
#                            'tissue_cell_type', 'organism', 'strategy', 'experiment_type'),columns=['studyid','tissue', 'organism', 'strategy', 'experiment_type'])
#         select_genedf = pd.DataFrame(clustergene.objects.filter(gene=gene_name).values_list('studyid', 'gene', 'cluster'),
#                               columns=['studyid', 'gene', 'cluster'])
#         select_modeldf = pd.DataFrame(ModelDF.objects.filter(model=model_type).values_list('studyid', 'cluster', 'model'),
#                                columns=['studyid', 'cluster', 'model'])
#         cluster = 'cluster0'
#         cluster0_data = select_genedf.loc[select_genedf['cluster'] == cluster]
#         leftcluster_data = select_genedf.loc[select_genedf['cluster'] != cluster]
#         # leftcluster_data = leftcluster_data.drop('parameters', axis=1)
#
#         # 创建新的 DataFrame 保存 cluster0 数据
#         new_modeldf = pd.DataFrame(columns=select_genedf.columns)
#
#         # 提取 cluster0_data 的相关列
#         cluster0_studyids = cluster0_data['studyid']
#         # cluster0_parameters = cluster0_data['parameters']
#         cluster0_gene = cluster0_data['gene']
#         # 构建新的 ModelDF 数据
#         new_modeldf['studyid'] = cluster0_studyids
#         new_modeldf['gene'] = cluster0_gene
#         new_modeldf['cluster'] = cluster
#         new_modeldf['model'] = 'constant'
#         # new_modeldf['parameters'] = cluster0_parameters
#         # new_modeldf['exp'] = 'y=c'
#
#         merged_df = leftcluster_data.merge(select_modeldf, on=['studyid', 'cluster'], how='inner')
#         result_df = merged_df[['studyid', 'cluster', 'model', 'gene']]
#         merged_df_gene = pd.concat([new_modeldf, result_df], ignore_index=True)
#         # part_results = select_genedf.merge(merged_df, on=['studyid', 'cluster'], how='left')
#         temp_results = select_sampleinfo.merge(merged_df_gene, on='studyid', how='left')
#         result_datatype_json = temp_results.to_dict('records')
#
#     return HttpResponse(json.dumps(result_datatype_json))



def plot_trend(request):
    studyid = request.POST.get('studyid')
    cluster = request.POST.get('cluster')
    gene = request.POST.get('gene')
    model = request.POST.get('model')

    basedir = settings.STATICFILES_DIRS[0] + 'data/'
    genedata_path = basedir + 'result_pkl/' + studyid + '_BrowseGene.pkl'
    trend_path = basedir + 'result_pkl/' + studyid + '_' + cluster + '_heatmapdata0409.pkl'
    if os.path.exists(genedata_path):
        genedata = pd.read_pickle(genedata_path)
        genedata['cluster'] = genedata['cluster'].str.replace('_', '')

    time = []
    trend_data = {}
    tread_data_line = {}
    x_Axis = []
    if cluster == 'cluster0':
        if os.path.exists(trend_path):
            trend_data_all = pd.read_pickle(trend_path)  # 不使用列索引
            trend_data_all = trend_data_all.set_index(trend_data_all.columns[0], drop=True)
            original_data = trend_data_all.loc[gene]
            c = original_data.mean()
            selected_data = genedata[(genedata['cluster'] == cluster) & (genedata['Genes'] == gene)]
            c_rounded = round(c, 4)
            cv_value = selected_data['CV'].values[0]
            cv_value_rounded = round(cv_value, 4)
            # cv_value = clustergene.objects.filter(studyid=studyid, cluster=cluster, gene=gene).values('cv').first()['cv']
            ModelDF_dd = {
                'cluster': cluster,
                'model': 'Constant',
                'parameters': f'c={c_rounded}',
                'exp': 'y=c',
                'cv': cv_value_rounded
            }
            ModelDF_df = [ModelDF_dd]

            time = list(trend_data_all.columns)[0:]
            gene_exp = list(trend_data_all.loc[gene])
            trend_data = {
                'x': time,
                'y': gene_exp,
                'z': [gene]
            }
            y_value = float(c_rounded)  # 假设 Parameters 包含固定的数值
            x_list = np.arange(1, len(time) + 0.001, 0.001)
            y_list = [y_value] * len(x_list)
            tread_data_line = {
                'x': x_list.tolist(),
                'y': y_list
            }
            x_Axis = list(range(1, len(time) + 1))
        else:
            trend_data = {
                'x': [],
                'y': [],
                'z': []
            }
            tread_data_line = {
                'x': [],
                'y': []
            }
    else:
        ModelDF_path = basedir + 'result_pkl/' + studyid + '_Generegression.pkl'
        # ModelDF_path = os.getcwd() + '/timecourse_app_static/data/result_pkl/' + studyid + '_ModelDF.pkl'
        # ModelDF_path = ModelDF_path.replace('/', '\\')
        if os.path.exists(ModelDF_path):
            ModelDF_dd = pd.read_pickle(ModelDF_path)
            ModelDF_dd = ModelDF_dd.drop_duplicates()
            ModelDF_dd['cluster'] = ModelDF_dd['cluster'].str.replace('_', '')
            ModelDF_dd.fillna('Null', inplace=True)
            ModelDF_dd = ModelDF_dd[ModelDF_dd['cluster'] == cluster]
            ModelDF_dd['Parameters'] = ModelDF_dd['Parameters'].apply(format_parameters)
            ModelDF_dd = ModelDF_dd.rename(columns={'Parameters': 'parameters', 'best_model': 'model'})
            ModelDF_dd = ModelDF_dd.replace({'Logis': 'Logistic', 'Sin': 'Trigonometric'})
            ModelDF_dd = ModelDF_dd[(ModelDF_dd['Genes'] == gene) & (ModelDF_dd['model'] == model)]
            ModelDF_df = ModelDF_dd.to_dict('records')
        else:
            ModelDF_df = {}

        if os.path.exists(trend_path):
            trend_data_all = pd.read_pickle(trend_path)  # 不使用列索引
            trend_data_all = trend_data_all.set_index(trend_data_all.columns[0], drop=True)
            time = list(trend_data_all.columns)[0:]
            if isinstance(trend_data_all, pd.DataFrame) and not isinstance(trend_data_all.loc[gene], pd.Series):
                gene_exp = list(trend_data_all.loc[gene].iloc[0].values)
            else:
                gene_exp = trend_data_all.loc[gene].tolist()
            trend_data = {
                'x': time,
                'y': gene_exp,
                'z': [gene]
            }
            # 先将参数赋值
            x_list = np.arange(1, len(time) + 0.001, 0.001)
            y_list = list()

            def get_y(x, Parameters_list):
                for para in Parameters_list:
                    exec(para)
                exec(Expressions)
                local_y = locals()['y']
                return local_y

            ###构建函数和虚拟的曲线值
            Expressions = ModelDF_df[0]['exp']
            Parameters = ModelDF_df[0]['parameters']

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
            Parameters_list = Parameters.split('; ')
            for k in range(0, len(x_list)):
                x = x_list[k]
                local_y = get_y(x, Parameters_list)
                y_list.append(local_y)

            tread_data_line = {
                'x': x_list.tolist(),
                'y': y_list,
                'z': [gene]
            }
            x_Axis = list(range(1, len(time) + 1))
        else:
            trend_data = {
                'x': [],
                'y': []
            }
            tread_data_line = {
                'x': [],
                'y': []
            }
    plot_data = {
        'trend_data': trend_data,
        'tread_data_line': tread_data_line,
        'x_Axis': x_Axis
    }
    # print(plot_data)
    return HttpResponse(json.dumps(plot_data))

def gene_detail(request,gene):
    studyid = request.GET.get('studyid')
    cluster = request.GET.get('cluster')
    gene = request.GET.get('gene')
    # gene='TPP'
    # studyid='TCD0071'
    # cluster='cluster2'

    cluster_type=cluster
    allsampleinfo = pd.DataFrame(
        Sampleinfo.objects.filter(studyid=studyid).values_list('studyid', 'time','unit','tissue_cell_type',
                                             'experiment_type','accession','pmid','data_source','organism'),
        columns=['studyid', 'time','unit','tissue_cell_type',
                                             'experiment_type','accession','pmid','data_source','organism']).drop_duplicates()

    # pkl_file = os.path.join(os.getcwd(), 'timecourse_app','static', 'data', 'result_pkl')
    basedir = settings.STATICFILES_DIRS[0] + 'data/'
    genedata_path = basedir + 'result_pkl/' + studyid + '_BrowseGene.pkl'
    # genedata_path = os.getcwd() + '/timecourse_app_static/data/result_pkl/' + studyid + '_BrowseGene.pkl'
    # genedata_path = genedata_path.replace("\\", "/")
    if os.path.exists(genedata_path):
        genedata = pd.read_pickle(genedata_path)
        genedata['cluster'] = genedata['cluster'].str.replace('_', '')
        genedf = genedata[(genedata['cluster'] == cluster)]
        genedf = genedf.rename(columns={'Studyid': 'studyid'})
        curr_study_data = allsampleinfo.merge(genedf, on='studyid', how='left')
        # leftcluster_df = curr_study_data.drop(['studyid'], axis=1)
        leftcluster_df = curr_study_data
        curr_study=leftcluster_df[leftcluster_df['Genes'] == gene]
        specie = leftcluster_df.loc[leftcluster_df['Genes'] == gene, 'organism'].unique()
        curr_study = curr_study.to_dict('records')

        curr_geneinfo=pd.DataFrame(geneinfo.objects.filter(gene=gene,specie__in=specie).values_list('geneid', 'gene',
                                                                                  'synonyms','description','specie'),
                                   columns=['geneid', 'gene', 'synonyms','description','specie'])
        curr_geneinfo = curr_geneinfo.to_dict('records')

    # trend_path = os.getcwd() + f'\\timecourse_app\\static\\data\\result_pkl\\{studyid}_{cluster}_heatmapdata.pkl'
    trend_path = basedir + 'result_pkl/' + studyid + '_' +cluster + '_heatmapdata0409.pkl'
    model_list=[]
    ModelDF_df=[]
    time=[]
    trend_data = {}
    tread_data_line = {}
    x_Axis = []
    # trend_path = os.getcwd() + '/timecourse_app_static/data/result_pkl/' + studyid + cluster + '_heatmapdata.pkl'
    # trend_path=trend_path.replace('/', '\\')
    if cluster=='cluster0':
        if os.path.exists(trend_path):
            trend_data_all = pd.read_pickle(trend_path)  # 不使用列索引
            trend_data_all = trend_data_all.set_index(trend_data_all.columns[0], drop=True)
            original_data = trend_data_all.loc[gene]
            c = original_data.mean()
            selected_data = genedata[(genedata['cluster'] == cluster) & (genedata['Genes'] == gene)]
            c_rounded = round(c, 4)
            cv_value = selected_data['CV'].values[0]
            cv_value_rounded = round(cv_value, 4)
            # cv_value = clustergene.objects.filter(studyid=studyid, cluster=cluster, gene=gene).values('cv').first()['cv']
            ModelDF_dd = {
                'cluster': cluster,
                'model': 'Constant',
                'parameters': f'c={c_rounded}',
                'exp': 'y=c',
                'cv': cv_value_rounded,
                'RMSE': '-',
                'Rsq': '-',
                'Adj_Rsq': '-',
                'AIC': '-',
            }
            ModelDF_df = [ModelDF_dd]

            time = list(trend_data_all.columns)[0:]
            gene_exp = list(trend_data_all.loc[gene])
            trend_data['Constant'] = {
                'x': time,
                'y': gene_exp,
                'z': [gene]
            }

            # c_value = Parameters.split('=')[1]
            y_value = float(c_rounded)  # 假设 Parameters 包含固定的数值
            x_list = np.arange(1, len(time) + 0.001, 0.001)
            y_list = [y_value] * len(x_list)
            tread_data_line['Constant'] = {
                'x': x_list.tolist(),
                'y': y_list
            }
            model_list=['Constant']
        else:
            trend_data['Constant'] = {
                'x': [],
                'y': {},
                'z': []
            }
            tread_data_line['Constant'] = {
                'x': [],
                'y': {}
            }
        x_Axis = list(range(1, len(time) + 1))
    else:
        ModelDF_path = basedir + 'result_pkl/' + studyid + '_Generegression.pkl'
        # ModelDF_path = os.getcwd() + '/timecourse_app_static/data/result_pkl/' + studyid + '_ModelDF.pkl'
        # ModelDF_path = ModelDF_path.replace('/', '\\')
        if os.path.exists(ModelDF_path):
            ModelDF_dd = pd.read_pickle(ModelDF_path)
            ModelDF_dd = ModelDF_dd.drop_duplicates()
            ModelDF_dd['cluster'] = ModelDF_dd['cluster'].str.replace('_', '')
            ModelDF_dd.fillna('Null', inplace=True)
            ModelDF_dd = ModelDF_dd[ModelDF_dd['cluster'] == cluster]
            ModelDF_dd['Parameters'] = ModelDF_dd['Parameters'].apply(format_parameters)
            ModelDF_dd = ModelDF_dd.rename(columns={'Parameters': 'parameters','best_model':'model'})
            ModelDF_dd = ModelDF_dd.replace({'Logis': 'Logistic', 'Sin': 'Trigonometric'})
            ModelDF_dd = ModelDF_dd[ModelDF_dd['Genes'] == gene]
            ModelDF_dd['RMSE'] = ModelDF_dd['RMSE'].apply(lambda x: round(x, 4))
            ModelDF_dd['Rsq'] = ModelDF_dd['Rsq'].apply(lambda x: round(x, 4))
            ModelDF_dd['Adj_Rsq'] = ModelDF_dd['Adj_Rsq'].apply(lambda x: round(x, 4))
            ModelDF_dd['AIC'] = ModelDF_dd['AIC'].apply(lambda x: round(x, 4))
            ModelDF_dd.replace(-float('inf'), -0.0001, inplace=True)
            model_list = ModelDF_dd['model'].tolist()
            ModelDF_df = ModelDF_dd.to_dict('records')
        else:
            ModelDF_df = {}

        if os.path.exists(trend_path):
            trend_data_all = pd.read_pickle(trend_path)  # 不使用列索引
            trend_data_all = trend_data_all.set_index(trend_data_all.columns[0], drop=True)
            # print(trend_data_all)
            for model in model_list:
                # 获取时间和数据值
                # trend_data = {}
                time = list(trend_data_all.columns)[0:]
                if isinstance(trend_data_all, pd.DataFrame) and not isinstance(trend_data_all.loc[gene], pd.Series):
                    gene_exp = list(trend_data_all.loc[gene].iloc[0].values)
                else:
                    gene_exp = trend_data_all.loc[gene].tolist()
                trend_data[model] = {
                    'x': time,
                    'y': gene_exp,
                    'z': [gene]
                }
                # 先将参数赋值
                x_list = np.arange(1, len(time) + 0.001, 0.001)
                y_list = list()
                def get_y(x, Parameters_list):
                    for para in Parameters_list:
                        exec(para)
                    exec(Expressions)
                    local_y = locals()['y']
                    return local_y
                ###构建函数和虚拟的曲线值
                Expressions = ModelDF_dd.loc[ModelDF_dd['model'] == model, 'exp'].iloc[0]
                Parameters = ModelDF_dd.loc[ModelDF_dd['model'] == model, 'parameters'].iloc[0]

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
                Parameters_list = Parameters.split('; ')
                for k in range(0, len(x_list)):
                    x = x_list[k]
                    local_y = get_y(x, Parameters_list)
                    y_list.append(local_y)

                tread_data_line[model] = {
                    'x': x_list.tolist(),
                    'y': y_list,
                    'z': [gene]
                }
            x_Axis = list(range(1, len(time) + 1))
        else:
            trend_data[model] = {
                'x': [],
                'y': {}
            }
            tread_data_line[model] = {
                'x': [],
                'y': {}
            }
    # print(trend_data)
    # print(tread_data_line)
    print(x_Axis)
    return render(request, 'gene_detail.html', {
        'curr_study': curr_study,
        'curr_geneinfo': curr_geneinfo,
        'model_list': model_list,
        'cluster': cluster_type,
        'trend_data': trend_data,
        'ModelDF_df': ModelDF_df,
        'tread_data_line': tread_data_line,
        'x_Axis': x_Axis,
        # 'heatmap_data': heatmap_data,
    })
def format_parameters(parameters_str):
    params = parameters_str.split('; ')
    formatted_params = []
    for param in params:
        key, value = param.split('=')
        formatted_value = f'{float(value):.4f}'
        formatted_params.append(f'{key} = {formatted_value}')
    return '; '.join(formatted_params)