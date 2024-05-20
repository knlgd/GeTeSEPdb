from django.shortcuts import render
from django.core.paginator import Paginator  ##Django分页组件
from django.core import serializers
from django.db.models import Q
from django.http import HttpResponse, HttpResponseBadRequest
from django.conf import settings
import pandas as pd
import time
import math, json
import logging

from timecourse_app.models import Sampleinfo, clustergene, ModelDF,geneinfo
import os, pickle
from django.conf import settings
from django.db import connection
from django.db import connections

# 声明全局变量
showdf = None
# print(showdf)

def load_showdf():
    global showdf
    # pkl_file = os.path.join(os.getcwd(), 'timecourse_app', 'static', 'data', 'search_showdata.pkl')
    basedir = settings.STATICFILES_DIRS[0] + 'data/'
    pkl_file = basedir + 'search_showdata.pkl'
    # pkl_file = pkl_file.replace("\\", "/")
    showdf = pd.read_pickle(pkl_file)
    showdf = showdf.fillna('Null')
    showdf['model'] = showdf['model'].str.capitalize()

def get_search_page(request):
    organism = request.POST.get('organism')

    return render(request, 'adv_search.html')

def get_experimenttype(request):
    organism = request.POST.get('organism')
    logging.info('get_experimenttype')
    querySet = Sampleinfo.objects.filter(organism=organism).values('experiment_type').distinct()
    result_experimenttype_json = list(map(lambda x: x['experiment_type'], querySet))

    # # organism = 'Mus_musculus'
    # # print(organism)
    # global showdf
    # # 检查showdf是否为None（未加载）
    # if showdf is None:
    #     load_showdf()
    # print(showdf[showdf['organism']==organism])
    # print(showdf[showdf['organism']==organism]['experiment_type'])
    #
    # # pkl_file = os.getcwd() + f'\\timecourse_app\\static\\data\\search_showdata.pkl'
    # # showdf = pd.read_pickle(pkl_file)
    # search_study_info = showdf[showdf['organism']==organism]['experiment_type'].drop_duplicates().dropna()
    # result_experimenttype_json = list(set(list(search_study_info)))
    # result_experimenttype_json.sort()
    print(result_experimenttype_json)
    return HttpResponse(json.dumps(result_experimenttype_json))



def get_modeltype(request):
    organism = request.POST.get('organism')
    experiment_type = request.POST.get('experiment_type')

    # # organism = 'Mus_musculus'
    # # experiment_type='Differentiation or development'
    # # print(organism)
    # global showdf
    # # 检查showdf是否为None（未加载）
    # if showdf is None:
    #     load_showdf()
    # # pkl_file = os.getcwd() + f'\\timecourse_app\\static\\data\\search_showdata.pkl'
    # # showdf = pd.read_pickle(pkl_file)
    # search_study_info = showdf[(showdf['organism']==organism)  & (showdf['experiment_type']==experiment_type)]['model'].drop_duplicates().dropna()
    # result_modeltype_json = list(set(list(search_study_info)))
    # result_modeltype_json.sort()
    sql = '''
    select DISTINCT 1 as id, g.model as model from timecourse_app_sampleinfo s left join timecourse_app_clustergene g on s.studyid = g.studyid
    WHERE s.organism = %s
    and s.experiment_type = %s
    and g.model is not null
    '''
    print(123)
    querySet = clustergene.objects.raw(sql, [organism, experiment_type])
    print(999)
    result_modeltype_json = list(map(lambda x: x.model, list(querySet)))
    return HttpResponse(json.dumps(result_modeltype_json))


def get_tissue(request):
    organism = request.POST.get('organism')
    experiment_type = request.POST.get('experiment_type')
    model_type = request.POST.get('model_type')
    # # print('organism',organism)
    # # print('experiment_type',experiment_type)
    # # print('model_type',model_type)
    # # organism = 'Mus_musculus'
    # # model_type = 'Logis'
    # print(organism,experiment_type,model_type)
    # global showdf
    # # 检查showdf是否为None（未加载）
    # if showdf is None:
    #     load_showdf()
    # # pkl_file = os.getcwd() + f'\\timecourse_app\\static\\data\\search_showdata.pkl'
    # # showdf = pd.read_pickle(pkl_file)
    # search_study_info = showdf[(showdf['organism']==organism) & (showdf['experiment_type']==experiment_type) & (showdf['model']==model_type)]['tissue_cell_type'].drop_duplicates()
    # # search_study_info = search_study_info.fillna('')
    # result_tissue_json = list(set(list(search_study_info)))
    # result_tissue_json.sort()
    # # print(result_tissue_json)

    sql = ''' 
    select DISTINCT 1 as studyid, s.tissue_cell_type as tissue_cell_type
    from timecourse_app_sampleinfo s
    left join timecourse_app_clustergene g on s.studyid = g.studyid
    WHERE s.organism = %s and s.experiment_type = %s
    and g.model = %s and s.tissue_cell_type is not null
    '''
    querySet = Sampleinfo.objects.raw(sql, [organism, experiment_type, model_type])
    result_tissue_json = list(map(lambda x: x.tissue_cell_type, list(querySet)))

    return HttpResponse(json.dumps(result_tissue_json))


def get_datatype(request):
    organism = request.POST.get('organism')
    experiment_type = request.POST.get('experiment_type')
    model_type = request.POST.get('model_type')
    tissue = request.POST.get('tissue_cell_type')
    # print(organism, studytype, tissue)
    # organism = 'Mus_musculus'
    # model_type = 'Logis'
    # tissue = 'Null'
    # global showdf
    # # 检查showdf是否为None（未加载）
    # if showdf is None:
    #     load_showdf()
    # # pkl_file = os.getcwd() + f'\\timecourse_app\\static\\data\\search_showdata.pkl'
    # # showdf = pd.read_pickle(pkl_file)
    #
    # search_study_info = showdf[(showdf['organism'] == organism) & (showdf['experiment_type']==experiment_type)
    #                            & (showdf['model'] == model_type) & (showdf['tissue_cell_type'] == tissue)]['strategy'].drop_duplicates().dropna()
    # result_datatype_json = list(set(list(search_study_info)))
    # result_datatype_json.sort()

    sql = ''' 
    select DISTINCT 1 as studyid , s.strategy as strategy
    from timecourse_app_sampleinfo s
    left join timecourse_app_clustergene g on s.studyid = g.studyid
    WHERE s.organism =  %s and s.experiment_type =  %s
    and g.model =  %s and s.tissue_cell_type =  %s and   s.strategy is not null
    '''
    querySet = Sampleinfo.objects.raw(sql, [organism, experiment_type, model_type, tissue])
    result_datatype_json = list(map(lambda x: x.strategy, list(querySet)))

    return HttpResponse(json.dumps(result_datatype_json))


def get_gene(request):
    organism = request.POST.get('organism')
    experiment_type = request.POST.get('experiment_type')
    model_type = request.POST.get('model_type')
    tissue = request.POST.get('tissue_cell_type')
    data_type = request.POST.get('data_type')

    # global showdf
    # # 检查showdf是否为None（未加载）
    # if showdf is None:
    #     load_showdf()
    # # pkl_file = os.getcwd() + f'\\timecourse_app\\static\\data\\search_showdata.pkl'
    # # showdf = pd.read_pickle(pkl_file)
    # # organism = 'Arabidopsis_thaliana'
    # # # gene_name = ''
    # # model_type = 'Logis'
    # # tissue = 'leaf'
    # # data_type = 'RNA-Seq'
    # # organism = [organism] if organism else []
    # # model_type = [model_type] if model_type else []
    # # tissue = [tissue] if tissue else []
    # # data_type = [data_type] if data_type else []
    #
    # search_study_info = showdf[(showdf['organism'] == organism) & (showdf['experiment_type']==experiment_type) & (showdf['model'] == model_type) &
    #                            (showdf['tissue_cell_type'] == tissue) & (showdf['strategy'] == data_type)]['gene'].drop_duplicates().dropna()
    # result_gene_json = list(set(list(search_study_info)))
    # result_gene_json.sort()

    sql = ''' 
    select DISTINCT 1 as id , g.gene as gene
    from timecourse_app_sampleinfo s
    left join timecourse_app_clustergene g on s.studyid = g.studyid
    WHERE s.organism =  %s and s.experiment_type =  %s
    and g.model =  %s and s.tissue_cell_type =  %s and   s.strategy = %s
    and g.gene is not null
    '''
    querySet = clustergene.objects.raw(sql, [organism, experiment_type, model_type, tissue, data_type ])
    result_gene_json = list(map(lambda x: x.gene, list(querySet)))

    return HttpResponse(json.dumps(result_gene_json))

def get_datatable(request):
    organism = request.POST.get('organism')
    experiment_type = request.POST.get('experiment_type')
    gene_name = request.POST.get('gene_name')
    model_type = request.POST.get('model_type')
    tissue = request.POST.get('tissue_cell_type')
    data_type = request.POST.get('data_type')

    # organism="Caenorhabditis_elegans"
    # experiment_type ="Differentiation or development"
    # gene_name = 'math-5'
    # model_type = "Trigonometric"
    # tissue = ''
    # data_type = ''
    # print(organism,gene_name,model_type,tissue,data_type)
    # start_time = time.time()
    # [{"studyid": "TCD0086", "cluster": "cluster1", "gene": "MEE49", "model": "Logistic", "accession": "GSE19271",
    #   "time": "49;49;53;53;57;57;61;61;65;65;69;69;73;73;77;77;81;81;85;85;89;89;93;93", "unit": "hour",
    #   "related_factor": "collected time", "tissue_cell_type": "seedling",
    #   "experiment_type": "Differentiation or development", "organism": "Arabidopsis_thaliana",
    #   "strategy": "Microarray"}]

    # organism= [organism] if organism else []
    # experiment_type = [experiment_type] if experiment_type else []
    # model_type = [model_type] if model_type else []
    # tissue = [tissue] if tissue else []
    # data_type = [data_type] if data_type else []
    # gene_name = [gene_name] if gene_name else []

    # global showdf
    # # 检查showdf是否为None（未加载）
    # if showdf is None:
    #     load_showdf()
    # # pkl_file = os.getcwd() + f'\\timecourse_app\\static\\data\\search_showdata.pkl'
    # # showdf = pd.read_pickle(pkl_file)
    #
    # # search_study_info = showdf[(showdf['organism'] == organism) & (showdf['model'] == model_type) & (showdf['tissue'] == tissue) & (
    # #                 showdf['gene'] == gene_name)]
    # print(organism,experiment_type,model_type,tissue,data_type,gene_name)
    # search_study_info = showdf[(showdf['organism'] == organism[0]) & (showdf['experiment_type'] == experiment_type[0])
    #                            & (showdf['model'] == model_type[0]) & (showdf['tissue_cell_type'] == tissue[0]) & (showdf['gene'] == gene_name[0])]
    # # search_study_info = showdf[(showdf['organism'] == organism) & (showdf['model'] == model_type) & (showdf['tissue'] == tissue) & (showdf['gene'] == gene_name)]
    # result_genetable_json = search_study_info.to_dict('records')
    sql = ''' 
        select DISTINCT 1 as id , 
        g.gene as gene, 
        g.experiment_type as experiment_type,
        g.organism as organism,
        g.strategy as strategy,
        g.tissue_cell_type as tissue_cell_type,
        g.cluster as cluster,
        g.model as model,
        g.studyid as studyid
        from timecourse_app_sampleinfo s
        left join timecourse_app_clustergene g on s.studyid = g.studyid
        WHERE 
        g.gene is not null
    '''
    params = []
    if organism is not None and len(organism) > 0:
        sql = sql + " and g.organism =  %s "
        params.append(organism)
    if experiment_type is not None and len(experiment_type) > 0:
        sql = sql + " and g.experiment_type =  %s "
        params.append(experiment_type)
    if data_type is not None and len(data_type) > 0:
        sql = sql + " and g.strategy =  %s "
        params.append(data_type)
    if tissue is not None and len(tissue) > 0:
        sql = sql + " and g.tissue_cell_type =  %s "
        params.append(tissue)
    if model_type is not None and len(model_type) > 0:
        sql = sql + " and g.model =  %s "
        params.append(model_type)
    if gene_name is not None and len(gene_name) > 0:
        sql = sql + " and g.gene =  %s "
        params.append(gene_name)

    sql = sql + ' limit 10000'

    rows = []
    with connection.cursor() as cursor:
        cursor.execute(sql, params)
        columns = [col[0] for col in cursor.description]
        rows = [
            dict(zip(columns, row))
            for row in cursor.fetchall()
        ]

    # querySet = clustergene.objects.raw(sql, [
    #     select_dict['organism'],
    #     select_dict['experiment_type'],
    #     select_dict['model'],
    #     select_dict['tissue_cell_type'],
    #     select_dict['strategy'],
    #     select_dict['gene']
    # ])
    # result_gene_json = list(map(lambda x: {'gene': x.gene, 'tissue_cell_type': x.tissue_cell_type}, list(querySet)))

    return HttpResponse(json.dumps(rows))

def search_adv_select_items(request):
    type = request.POST.get('type')
    result_json = []
    logging.info('get_experimenttype')
    if type == 'organism':
        # get_datatable(request)
        result_json = do_get_sample_info(request, 'organism')
        # result_json = list(map(lambda x: x.organism, list(querySet)))
    if type == 'experiment_type':
        result_json = do_get_sample_info(request, 'experiment_type')
        # result_json = list(map(lambda x: x.experiment_type, list(querySet)))
    if type == 'tissue':
        result_json = do_get_sample_info(request, 'tissue_cell_type')
        # result_json = list(map(lambda x: x.tissue_cell_type, list(querySet)))
    if type == 'data_type':
        result_json = do_get_sample_info(request, 'strategy')
        # result_json = list(map(lambda x: x.strategy, list(querySet)))
    if type == 'model_type':
        result_json = do_get_sample_info(request, 'model')
        # result_json = list(map(lambda x: x.model, list(querySet)))
    if type == 'gene':
        result_json = do_get_sample_info(request, 'gene')
        # result_json = list(map(lambda x: x.model, list(querySet)))
    return HttpResponse(json.dumps({'items': result_json}), content_type='application/json')


def do_get_sample_info(request, column_name):
    organism = request.POST.get('organism')
    experiment_type = request.POST.get('experiment_type')
    tissue = request.POST.get('tissue')
    data_type = request.POST.get('data_type')
    model_type = request.POST.get('model_type')
    input_gene = request.POST.get('input_')

    # sql = 'select DISTINCT g.studyid, g.' + column_name + ' as ' + column_name \
    #       + ' from cluster_gene g WHERE  g.' + column_name + ' is not null '
    #
    # queryParams = []
    # if organism is not None and len(organism) > 0:
    #     queryParams.append(organism)
    #     sql = sql + ' and g.organism = %s '
    # if experiment_type is not None and len(experiment_type) > 0:
    #     queryParams.append(experiment_type)
    #     sql = sql + ' and g.experiment_type = %s '
    # if tissue is not None and len(tissue) > 0:
    #     queryParams.append(tissue)
    #     sql = sql + ' and g.tissue_cell_type = %s '
    # if data_type is not None and len(data_type) > 0:
    #     queryParams.append(data_type)
    #     sql = sql + ' and g.strategy = %s '
    # if model_type is not None and len(model_type) > 0:
    #     queryParams.append(model_type)
    #     sql = sql + ' and g.model = %s  '
    # if model_type is not None and len(model_type) > 0:
    #     queryParams.append(model_type)
    #     sql = sql + ' and g.model = %s  '

    # sql = '''
    #         select DISTINCT
    #         g.gene as gene,
    #         g.experiment_type as experiment_type,
    #         g.organism as organism,
    #         g.strategy as strategy,
    #         g.tissue_cell_type as tissue_cell_type,
    #         g.cluster as cluster,
    #         g.model as model,
    #         from  timecourse_app_clustergene
    #         WHERE g.gene is not null
    #     '''
    sql = 'select DISTINCT g.' + column_name + ' from cluster_gene g where g.' + column_name + ' is not null'
    params = []
    if organism is not None and len(organism) > 0:
        sql = sql + " and g.organism =  %s "
        params.append(organism)
    if experiment_type is not None and len(experiment_type) > 0:
        sql = sql + " and g.experiment_type =  %s "
        params.append(experiment_type)
    if data_type is not None and len(data_type) > 0:
        sql = sql + " and g.strategy =  %s "
        params.append(data_type)
    if tissue is not None and len(tissue) > 0:
        sql = sql + " and g.tissue_cell_type =  %s "
        params.append(tissue)
    if model_type is not None and len(model_type) > 0:
        sql = sql + " and g.model =  %s "
        params.append(model_type)
    if input_gene is not None and len(input_gene) > 0:
        sql = sql + " and g.gene  like concat('%%', %s, '%%') "
        params.append(input_gene)

    print(sql)
    rows = []
    with connections['postgresql'].cursor() as cursor:
        cursor.execute(sql, params)
        columns = [col[0] for col in cursor.description]
        rows = [
            dict(zip(columns, row))
            for row in cursor.fetchall()
        ]

    data = []



    for d in rows:
        data.append(d[column_name])
    return data

def do_get_gene_info(request, column_name):
    organism = request.POST.get('organism')
    experiment_type = request.POST.get('experiment_type')
    tissue = request.POST.get('tissue')
    data_type = request.POST.get('data_type')
    model_type = request.POST.get('model_type')
    sql = 'select DISTINCT g.id as id, g.' + column_name + ' as ' + column_name \
          + ' from  timecourse_app_clustergene g  WHERE  g.id is not null '

    queryParams = []
    if organism is not None and len(organism) > 0:
        queryParams.append(organism)
        sql = sql + ' and g.organism = %s '
    if experiment_type is not None and len(experiment_type) > 0:
        queryParams.append(experiment_type)
        sql = sql + ' and g.experiment_type = %s '
    if tissue is not None and len(tissue) > 0:
        queryParams.append(tissue)
        sql = sql + ' and g.tissue_cell_type = %s '
    if data_type is not None and len(data_type) > 0:
        queryParams.append(data_type)
        sql = sql + ' and g.strategy = %s '
    if model_type is not None and len(model_type) > 0:
        queryParams.append(model_type)
        sql = sql + ' and g.model = %s  '

    return clustergene.objects.raw(sql, queryParams)

def get_gene_names(request):
    organism = request.POST.get('organism')
    experiment_type = request.POST.get('experiment_type')
    model_type = request.POST.get('model_type')
    tissue = request.POST.get('tissue')
    data_type = request.POST.get('data_type')
    input_= request.POST.get('input_')
    data = do_get_sample_info(request, 'gene')
    return HttpResponse(json.dumps({'data_list': data}))
    # organism = 'Arabidopsis_thaliana'
    # experiment_type = 'Differentiation or development'
    # model_type = 'Logistic'
    # tissue = 'seedling'
    # data_type = 'Microarray'
    # input_ = 'MGP1'


    # global showdf
    # # 检查showdf是否为None（未加载）
    # if showdf is None:
    #     load_showdf()
    # search_study_info = showdf[(showdf['organism'] == organism) & (showdf['experiment_type'] == experiment_type) & (showdf['model'] == model_type) &
    #                            (showdf['tissue_cell_type'] == tissue) & (showdf['strategy'] == data_type)]['gene'].drop_duplicates().dropna()
    # data_list = list(set(list(search_study_info)))
    #
    # print(data_list)
    # top_data_list = []
    # for i in data_list:
    #     if i.upper().find(input_.upper()) != -1:
    #         top_data_list.append(i)
    #     if len(top_data_list) == 30:
    #         break
    # # print(1234156)
    # print(top_data_list)

    # sql = '''
    # select DISTINCT 1 as id , g.gene as gene
    # from timecourse_app_sampleinfo s
    # left join timecourse_app_clustergene g on s.studyid = g.studyid
    # WHERE
    # g.gene is not null
    # '''
    # query_params = []
    # if organism is not None and len(organism) > 0:
    #     query_params.append(organism)
    #     sql = sql + 'and s.organism =  %s '
    # if experiment_type is not None and len(experiment_type) > 0:
    #     query_params.append(experiment_type)
    #     sql = sql + ' and s.experiment_type =  %s '
    # if model_type is not None and len(model_type) > 0:
    #     query_params.append(model_type)
    #     sql = sql + ' and g.model =  %s '
    # if tissue is not None and len(tissue) > 0:
    #     query_params.append(tissue)
    #     sql = sql + ' and s.tissue_cell_type =  %s  '
    # if data_type is not None and len(data_type) > 0:
    #     query_params.append(data_type)
    #     sql = sql + ' and s.strategy = %s  '
    # if input_ is not None and len(input_) > 0:
    #     query_params.append(input_)
    #     sql = sql + " and g.gene like concat('%%', %s, '%%')  "
    #
    # sql = sql + " limit 20"
    # querySet = clustergene.objects.raw(sql, query_params)
    # top_data_list = list(map(lambda x: x.gene, list(querySet)))
    # return HttpResponse(json.dumps({'data_list': top_data_list}))
