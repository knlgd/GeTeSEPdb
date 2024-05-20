import math
import os
import shutil
import pickle
import pandas as pd
import json
from django.http import FileResponse, HttpResponse, Http404, StreamingHttpResponse
from django.shortcuts import render
from django.core.files.storage import FileSystemStorage
from djangoProject import settings
from concurrent.futures.thread import ThreadPoolExecutor
import uuid
import numpy as np
import time
from django.core.mail import send_mail
import re
import chardet
from zipfile import ZipFile

# 定义全局线程池
global_thread_pool = ThreadPoolExecutor(100)

def get_analysis_page(request):
    return render(request, 'analysis.html')

def single_gene_analysis_page(request):
    return render(request, "analysis_singlegene.html")

def analysis_retrieve_page(request):
    return render(request, "analysis_retrieve.html")


def single_gene_analysis_submit(request):
    exp_file = request.FILES.get('file1')
    info_file = request.FILES.get('file2')
    timeSeriesTextOpen = request.POST.get("timeSeriesTextOpen")
    timeInfoTextOpen = request.POST.get("timeInfoTextOpen")
    analysis_id = request.POST.get("analysisResultId")
    timeSeriesText = request.POST.get("timeSeriesText")
    timeInfoText = request.POST.get("timeInfoText")

    temp_cache_directory = settings.STATICFILES_DIRS[0] + 'cache/' + analysis_id + '/'
    if analysis_id == 'example':
        analysis_id = 'example'
        temp_cache_directory = settings.STATICFILES_DIRS[0] + 'exampledata/' + analysis_id + '/'
    elif len(analysis_id) <= 0:
        analysis_id = str(uuid.uuid1())
        temp_cache_directory = settings.STATICFILES_DIRS[0] + 'cache/' + analysis_id + '/'
    data = {
        # status 1 代表成功
        'status': 0,
        'msg': "失败",
        'data': {
            'analysis_id': analysis_id
        }
    }
    try:
        # 文件保存的临时目录
        temp_cache_directory = temp_cache_directory.replace("\\", "/")
        gene_filepath = temp_cache_directory + 'gene_data.txt'
        info_filepath = temp_cache_directory + 'sampleinfo_data.txt'
        # 将文件保存到临时目录
        if not os.path.exists(temp_cache_directory):
            os.makedirs(temp_cache_directory)
            if not os.path.exists(temp_cache_directory + 'type_single_gene'):
                f = open(temp_cache_directory + 'type_single_gene', "w")
                f.close()
            with open(temp_cache_directory + 'gene_data.txt', 'wb') as file:
                if timeSeriesTextOpen == 'true':
                    file.write(bytes(timeSeriesText, encoding='UTF-8'))
                else:
                    byteData = exp_file.read()
                    encoding = chardet.detect(byteData)['encoding']
                    file.write(bytes(str(byteData, encoding=encoding), encoding='UTF-8'))
            with open(temp_cache_directory + 'sampleinfo_data.txt', 'wb') as file:
                if timeInfoTextOpen == 'true':
                    file.write(bytes(timeInfoText, encoding='UTF-8'))
                else:
                    byteData = info_file.read()
                    encoding = chardet.detect(byteData)['encoding']
                    file.write(bytes(str(byteData, encoding=encoding), encoding='UTF-8'))

            if not os.path.exists(temp_cache_directory + 'result.zip'):
                r_script_path = settings.STATICFILES_DIRS[0] + 'script/fit_function.R'
                command = f'Rscript "{r_script_path}" "{gene_filepath}" "{info_filepath}" "{temp_cache_directory}"'
                os.system(command)
                with ZipFile(os.path.join(temp_cache_directory,'result.zip'), 'w') as zip:
                    for root,dirs,files in os.walk(temp_cache_directory):
                        for file in files:
                            filepath=os.path.join(root,file)
                            arcname=os.path.relpath(filepath,temp_cache_directory)
                            zip.write(filepath,arcname)
                print('successfully!')
        data = get_single_gene_result_json(analysis_id, temp_cache_directory)
    except Exception as ee:
        print(ee)


    return HttpResponse(json.dumps(data), content_type='application/json')


# 获取single gene result json
def get_single_gene_result_json(analysis_id, temp_cache_directory):
    gene_filepath = temp_cache_directory + 'gene_data.txt'
    info_filepath = temp_cache_directory + 'sampleinfo_data.txt'
    EvaIndicator_path = temp_cache_directory + 'EvaIndicators.csv'
    EvaIndicator_path = EvaIndicator_path.replace("\\", "/")
    if os.path.exists(EvaIndicator_path):
        EvaIndicator_data = pd.read_csv(EvaIndicator_path)
        EvaIndicator_data.columns = ['model', 'rmse', 'rsq', 'adjrsq','aicc']
        EvaIndicator_data = EvaIndicator_data.replace({'exponential': 'Exponential'})
        EvaIndicator_data['rmse'] = EvaIndicator_data['rmse'].apply(lambda x: round(x, 4))
        EvaIndicator_data['rsq'] = EvaIndicator_data['rsq'].apply(lambda x: round(x, 4))
        EvaIndicator_data['adjrsq'] = EvaIndicator_data['adjrsq'].apply(lambda x: round(x, 4))
        EvaIndicator_data['aicc'] = EvaIndicator_data['aicc'].apply(lambda x: round(x, 4))
        EvaIndicator_data.replace(-float('inf'), -0.0001, inplace=True)
        EvaIndicator_de = {
            'model': 'Constant',
            'rmse': '-',
            'rsq': '-',
            'adjrsq': '-',
            'aicc': '-',
        }
        EvaIndicator_dict = pd.DataFrame([[EvaIndicator_de['model'],EvaIndicator_de['rmse'], EvaIndicator_de['rsq'], EvaIndicator_de['adjrsq'],
                                           EvaIndicator_de['aicc']]],
                                         columns=['model', 'rmse', 'rsq', 'adjrsq', 'aicc'])
        EvaIndicator_data1=pd.concat([EvaIndicator_dict, EvaIndicator_data])
        EvaIndicator_df = EvaIndicator_data1.to_dict('records')
    else:
        EvaIndicator_df = {}

    ConstantDF_path = temp_cache_directory + 'ConstantDF.csv'
    ConstantDF_path = ConstantDF_path.replace("\\", "/")
    ModelDF_path = temp_cache_directory + 'ModelDF.csv'
    ModelDF_path = ModelDF_path.replace("\\", "/")
    trend_data={}
    tread_data_line={}
    if os.path.exists(ConstantDF_path):
        ConstantDF_data = pd.read_csv(ConstantDF_path)
        ConstantDF_data['mean'] = ConstantDF_data['mean'].apply(lambda x: round(x, 4))
        ConstantDF_result = ConstantDF_data.to_dict('records')
        ModelDF_result = {}
        time = pd.read_csv(info_filepath, delimiter='\t')
        gene_exp = pd.read_csv(gene_filepath, delimiter='\t')

        time_list = [int(col) for col in time.columns[0:].tolist()]
        gene_exp_list = [float(col) for col in gene_exp.columns[1:].tolist()]

        # Prepare trend_data
        trend_data['Constant'] = {
            'x': time_list,
            'y': gene_exp_list,
            'z': gene_exp.columns[0]
        }
        Parameters = ConstantDF_data['mean']
        x_list = np.arange(1, len(time_list) + 0.01, 0.01)
        y_list = np.tile(Parameters, len(x_list))

        tread_data_line['Constant'] = {
            'x': x_list.tolist(),
            'y': y_list,
            'z': gene_exp.columns[0]
        }
        x_Axis = list(range(1, len(time_list) + 1))

    elif os.path.exists(ModelDF_path):
        ConstantDF_result = {}
        # Load ModelDF_data, time, and gene_exp
        ModelDF_data = pd.read_csv(ModelDF_path)
        ModelDF_data.columns = ['model', 'parameters', 'exp']
        ModelDF_data['parameters'] = ModelDF_data['parameters'].apply(format_parameters)
        ModelDF_result = ModelDF_data.to_dict('records')
        model_list = ModelDF_data['model'].tolist()
        time = pd.read_csv(info_filepath, delimiter='\t')
        gene_exp = pd.read_csv(gene_filepath, delimiter='\t')

        time_list = [int(col) for col in time.columns[0:].tolist()]
        gene_exp_list = [float(col) for col in gene_exp.columns[1:].tolist()]

        for model in model_list:
            # Prepare trend_data
            trend_data[model] = {
                'x': time_list,
                'y': gene_exp_list,
                'z': gene_exp.columns[0]
            }

            Expressions = ModelDF_data[ModelDF_data['model'] == model]['exp'].tolist()
            Parameters = ModelDF_data[ModelDF_data['model'] == model]['parameters'].tolist()

            # Modify Expressions and Parameters
            for i in range(len(Expressions)):
                Expressions[i] = Expressions[i].replace('exp', 'math.exp')
                Expressions[i] = Expressions[i].replace('log', 'math.log')
                Expressions[i] = Expressions[i].replace('power', 'math.pow')
                Expressions[i] = Expressions[i].replace('sin', 'math.sin')
                Expressions[i] = Expressions[i].replace('pi', str(math.pi))
                Expressions[i] = Expressions[i].replace(' ', '')
                Expressions[i] = Expressions[i].replace('^', '**')

            Parameters_list = [para.split('; ') for para in Parameters]

            x_list = np.arange(1, len(time_list) + 0.01, 0.01)
            y_list = []

            def get_y(x, para_list):
                for para in para_list:
                    exec(para)
                exec(Expressions[i])
                local_y = locals()['y']
                return local_y

            for k in range(0, len(x_list)):
                x = x_list[k]
                local_y_list = []  # Create an empty list for each x value
                for i in range(len(Parameters_list)):
                    local_y = get_y(x, Parameters_list[i])
                    y_list.append(local_y)

            tread_data_line[model] = {
                'x': x_list.tolist(),
                'y': y_list,
                'z': gene_exp.columns[0]
            }

            x_Axis = list(range(1, len(time_list) + 1))
    else:
        ConstantDF_result = {}
        ModelDF_result = {}
        trend_data = {}
        tread_data_line = {}

    plot_data = {
        'trend_data': trend_data,
        'tread_data_line': tread_data_line,
        'x_Axis': x_Axis
    }
    data = {
        # status 1 代表成功
        'status': 1,
        'msg': "成功",
        'data': {
            'analysisType': 'singleGene',
            'analysis_id': analysis_id,
            'model_list':model_list,
            'ConstantDF_result': ConstantDF_result,
            'ModelDF_df': ModelDF_result,
            'EvaIndicator_df': EvaIndicator_df,
            'plot_data': plot_data,
        }
    }
    return data


def retrieve_result(request):
    analysis_id = request.POST.get("analysisResultId")

    temp_cache_directory = settings.STATICFILES_DIRS[0] + 'cache/' + analysis_id + '/'

    if analysis_id == 'example' or analysis_id == 'example2':
        temp_cache_directory = settings.STATICFILES_DIRS[0] + 'exampledata/' + analysis_id + '/'


    data = {
        # status 1 代表成功
        'status': 0,
        'msg': "结果不存在",
        'data': {
            'analysis_id': analysis_id,
            'ConstantDF_result': [],
            'ModelDF_df': [],
            'EvaIndicator_df': [],
            'plot_data': [],
        }
    }

    if os.path.exists(temp_cache_directory + 'result.zip'):
        if analysis_id == 'example' or os.path.exists(temp_cache_directory + '/type_single_gene'):
            data = get_single_gene_result_json(analysis_id, temp_cache_directory)
        elif analysis_id == 'example2' or (not os.path.exists(temp_cache_directory + '/type_single_gene')):
            data = get_omic_result_json(analysis_id)


    return HttpResponse(json.dumps(data), content_type='application/json')

def get_analysis_result_page(request):
    heatmap_data = {}
    cluster_list = {}
    trend_data = {}
    EvaIndicator_df = {}
    ModelDF_result = {}
    analysisId = request.GET.get("analysisId")
    return render(request, 'analysis_result.html', {
        'analysisId': analysisId,
        'heatmap_data': heatmap_data,
        'cluster_list': cluster_list,
        'trend_data': trend_data,
        'EvaIndicator_df': EvaIndicator_df,
        'ModelDF_df': ModelDF_result,
        # 'result_zip': result_zip,
    })

def get_analysis_result_data(request):

    # 分析ID
    analysisId = request.GET.get("analysisId")

    data = get_omic_result_json(analysisId)
    return HttpResponse(json.dumps(data), content_type='application/json')


def get_omic_result_json(analysisId):
    # 文件保存的临时目录
    temp_cache_directory = settings.STATICFILES_DIRS[0] + 'cache/' + analysisId + '/'
    if analysisId == 'example2':
        temp_cache_directory = settings.STATICFILES_DIRS[0] + 'exampledata/' + analysisId + '/'
    temp_cache_directory = temp_cache_directory.replace("\\", "/")
    # 通过分析ID从结果目录中提取数据
    heatmap_data = {}
    cluster_list = {}
    ## 数据提取
    EvaIndicator_path = temp_cache_directory + 'EvaIndicators.csv'
    EvaIndicator_path = EvaIndicator_path.replace("\\", "/")
    cluster = 'cluster_0'
    genedata_path = temp_cache_directory + 'BrowseGene.csv'
    genedata_path = genedata_path.replace("\\", "/")
    if os.path.exists(genedata_path) and os.path.exists(EvaIndicator_path):
        genedata = pd.read_csv(genedata_path)
        genedata['mean'] = genedata['mean'].apply(lambda x: round(x, 4))
        genedata['CV'] = genedata['CV'].apply(lambda x: round(x, 4))
        genedata['cluster'] = genedata['cluster'].str.replace('_', '')
        genedf = genedata[(genedata['cluster'] == cluster)]
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
            'aicc': '-',
            'mean': c,
            'CV': e
        }
        parameters_strc = str(EvaIndicator_de['mean'])  # 转换parameters为字符串形式
        parameters_stre = str(EvaIndicator_de['CV'])
        EvaIndicator_dict = pd.DataFrame([[EvaIndicator_de['model'], EvaIndicator_de['cluster'],
                                           EvaIndicator_de['rmse'], EvaIndicator_de['rsq'],EvaIndicator_de['adjrsq'], EvaIndicator_de['aicc'],
                                           parameters_strc, parameters_stre]],
                                         columns=['model', 'cluster', 'rmse', 'rsq','adjrsq', 'aicc', 'mean', 'CV'])

        EvaIndicator_data = pd.read_csv(EvaIndicator_path)
        EvaIndicator_data['rmse'] = EvaIndicator_data['rmse'].apply(lambda x: round(x, 4))
        EvaIndicator_data['rsq'] = EvaIndicator_data['rsq'].apply(lambda x: round(x, 4))
        EvaIndicator_data['adjrsq'] = EvaIndicator_data['adjrsq'].apply(lambda x: round(x, 4))
        EvaIndicator_data['aicc'] = EvaIndicator_data['aicc'].apply(lambda x: round(x, 4))
        EvaIndicator_data['mean'] = "-"
        EvaIndicator_data['CV'] = "-"
        EvaIndicator_data.replace(-float('inf'), -0.0001, inplace=True)
        EvaIndicator_combine = pd.concat([EvaIndicator_dict, EvaIndicator_data])
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
                                              EvaIndicator_de['rmse'], EvaIndicator_de['rsq'],EvaIndicator_de['adjrsq'], EvaIndicator_de['aic'],
                                              parameters_strc, parameters_stre, EvaIndicator_de['Studyid']]],
                                            columns=['model', 'cluster', 'rmse', 'rsq','adjrsq', 'aic', 'mean', 'CV', 'studyid'])
        EvaIndicator_df = EvaIndicator_combine.to_dict('records')
    else:
        EvaIndicator_df = {}
    # if os.path.exists(EvaIndicator_path):
    #     EvaIndicator_data = pd.read_csv(EvaIndicator_path)
    #     # EvaIndicator_data = EvaIndicator_data[EvaIndicator_data['model'] != 'exponential_decay']
    #     # EvaIndicator_data.columns = ['model', 'rmse', 'rsq', 'aic', 'cluster']
    #     # EvaIndicator_data = EvaIndicator_data.replace(
    #     #     {'Logis': 'Logistic', 'Sin': 'trigonometric', 'exponential_growth': 'Exponential'})
    #     EvaIndicator_data['rmse'] = EvaIndicator_data['rmse'].apply(lambda x: round(x, 4))
    #     EvaIndicator_data['rsq'] = EvaIndicator_data['rsq'].apply(lambda x: round(x, 4))
    #     EvaIndicator_data['aicc'] = EvaIndicator_data['aicc'].apply(lambda x: round(x, 4))
    #     EvaIndicator_df = EvaIndicator_data.to_dict('records')
    # else:
    #     EvaIndicator_df = {}
    cluster = 'cluster_0'
    ModelDF_path = temp_cache_directory + 'ModelDF.csv'
    ModelDF_path = ModelDF_path.replace("\\", "/")
    heatmap_path = temp_cache_directory + cluster + '_heatmapdata.tsv'
    heatmap_path = heatmap_path.replace("\\", "/")
    ModelDF_combine={}
    if os.path.exists(ModelDF_path) and os.path.exists(heatmap_path):
        ModelDF_data = pd.read_csv(ModelDF_path)
        ModelDF_data.columns = ['cluster', 'model', 'parameters', 'exp']
        ModelDF_data['parameters'] = ModelDF_data['parameters'].apply(format_parameters)
        ####cluster0 ModelDF_df
        # temp_cache_directory=r'F:/project/newproject/timecourse_app/static/cache/'
        cluster0_data = pd.read_csv(heatmap_path, delimiter='\t', index_col=0)

        # cluster0_data.set_index('Genes', inplace=True)
        gene_means = cluster0_data.mean(axis=1)
        # 计算所有基因均值的均值
        overall_mean = gene_means.mean()
        formatted_mean_c = round(overall_mean, 4)
        # genedf = pd.read_csv(temp_cache_directory + 'BrowseGene.csv')
        # selected_df = genedf[genedf['cluster'] == cluster]
        # selected_df['mean'] = selected_df['mean'].apply(lambda x: round(x, 4))
        # selected_df['CV'] = selected_df['CV'].apply(lambda x: round(x, 4))
        # c = list(selected_df['CV'])
        # mean_c = sum(c) / len(c)
        # formatted_mean_c = round(mean_c, 4)
        # if len(c) > 5:
        #     c = c[:5]
        #     c[-1] = str(c[-1]) + '...'
        #c = round(c, 4)
        ModelDF_dc = {
            'cluster': cluster,
            'model': 'constant',
            'parameters': formatted_mean_c,
            'exp': 'y=c'
        }
        parameters_str = formatted_mean_c  # 转换parameters为字符串形式
        ModelDF_dc_dict = pd.DataFrame(
            [[ModelDF_dc['cluster'], ModelDF_dc['model'], parameters_str, ModelDF_dc['exp']]],
            columns=['cluster', 'model', 'parameters', 'exp'])
        ModelDF_combine = pd.concat([ModelDF_dc_dict, ModelDF_data])

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
        cluster0_data = pd.read_csv(heatmap_path, delimiter='\t', index_col=0)
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
        ModelDF_result = {}
    # result_zip = None
    heatmap_data = {}
    trend_data = {}
    tread_data_line = {}
    x_Axis = []
    cluster_model_list = [f"{cluster}-{model}" for cluster, model in
                          zip(ModelDF_combine['cluster'], ModelDF_combine['model'])]
    for cluster_model in cluster_model_list:
        cluster = cluster_model.split('-')[0]
        model = cluster_model.split('-')[1]
        heatmap_path = temp_cache_directory + cluster + '_heatmapdata.tsv'
        heatmap_path = heatmap_path.replace("\\", "/")
        if os.path.exists(heatmap_path):
            data = pd.read_csv(heatmap_path, delimiter='\t', index_col=0)
            heatmap_data[cluster_model] = {
                'x': list(data.columns),
                'y': list(data.index),
                'z': data.values.tolist()
            }
        else:
            heatmap_data[cluster_model] = None

        if os.path.exists(heatmap_path):
            trend_data_all = pd.read_csv(heatmap_path, delimiter='\t', index_col=0)
            time = list(trend_data_all.columns)[0:]
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
            x_list = np.arange(1, len(time) + 0.001, 0.001)
            y_list = list()

            def get_y(x, Parameters_list):
                for para in Parameters_list:
                    exec(para)
                exec(Expressions)
                local_y = locals()['y']
                return local_y

            # 先将参数赋值
            if cluster == 'cluster_0':
                Parameters_list = [Parameters]
                y_value = float(Parameters)  # 假设 Parameters 包含固定的数值
                y_list = [y_value] * len(x_list)
                tread_data_line[cluster_model] = {
                    'x': x_list.tolist(),
                    'y': y_list
                }
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
    # 返回json数据给页面
    status = 0
    msg = '运行中'
    if os.path.exists(temp_cache_directory +'rc'):
        status = 2
    if os.path.exists(temp_cache_directory + 'result.zip'):
        status = 1
        msg = '成功'
    data = {
        # status 1 代表成功
        'status': status,
        'msg': msg,
        'data': {
            'analysisType': 'omicData',
            'analysis_id': analysisId,
            'heatmap_data': heatmap_data,
            'tread_data_line': tread_data_line,
            'cluster_model_list': cluster_model_list,
            'trend_data': trend_data,
            'x_Axis': x_Axis,
            'plot_data': {
                'heatmap_data': heatmap_data,
                'trend_data': trend_data,
                'tread_data_line': tread_data_line,
                'cluster_list': cluster_list
            },
            'EvaIndicator_df': EvaIndicator_df,
            'ModelDF_df': ModelDF_result,
        }
    }
    return data


def submit_analysis(request):
    # 生成分析ID
    analysis_id = str(uuid.uuid1())
    endpoint = '{scheme}://{host}'.format(scheme=request.scheme, host=request.get_host())


    exp_file = request.FILES.get('file1')
    info_file = request.FILES.get('file2')
    timeSeriesTextOpen = request.POST.get("timeSeriesTextOpen")
    timeInfoTextOpen = request.POST.get("timeInfoTextOpen")
    timeSeriesText = request.POST.get("timeSeriesText")
    timeInfoText = request.POST.get("timeInfoText")

    # 邮件
    to_mail = request.POST.get('receiveMail')

    # 文件保存的临时目录
    temp_cache_directory = settings.STATICFILES_DIRS[0] + 'cache/' + analysis_id + '/'
    temp_cache_directory = temp_cache_directory.replace("\\", "/")

    # 将文件保存到临时目录
    if not os.path.exists(temp_cache_directory):
        os.makedirs(temp_cache_directory)
    with open(temp_cache_directory + 'gene_data.txt', 'wb') as file:
        if timeSeriesTextOpen == 'true':
            file.write(bytes(timeSeriesText, encoding='UTF-8'))
        else:
            byteData = exp_file.read()
            encoding = chardet.detect(byteData)['encoding']
            file.write(bytes(str(byteData, encoding=encoding), encoding='UTF-8'))
    with open(temp_cache_directory + 'sampleinfo_data.txt', 'wb') as file:
        if timeInfoTextOpen == 'true':
            file.write(bytes(timeInfoText, encoding='UTF-8'))
        else:
            byteData = info_file.read()
            encoding = chardet.detect(byteData)['encoding']
            file.write(bytes(str(byteData, encoding=encoding), encoding='UTF-8'))

    print('保存的临时目录为：' + temp_cache_directory)

    data = {
        # status 1 代表成功
        'status': 1,
        'msg': '成功了',
        'data': {
            # 分析ID
            'analysis_id': analysis_id
        }
    }

    # 提交任务到线程池
    global_thread_pool.submit(asyncHandle, analysis_id, endpoint, to_mail)
    return HttpResponse(json.dumps(data), content_type='application/json')

# 异步处理的方法
def asyncHandle(analysis_id, endpoint, to_mail):
    print('任务ID:' + analysis_id + '分析完成,tomail=' + to_mail)

    # if re.match('^.*?@.*', to_mail):
    #     subject = '任务运行完成'
    #     content = '<a style="font-size: 28px; font-weight: 700; text-align: center;" href="' + endpoint + '/timecourse/get_analysis_result_page?analysisId=' + analysis_id + '">--》点击查看运行结果《--</a>'
    #     from_mail = '任务运行通知<' + settings.EMAIL_HOST_USER + '>'
    #     try:
    #         send_mail(subject,message=None, from_email=from_mail, recipient_list=[to_mail], fail_silently=False, html_message=content)
    #         print('发件邮件成功！')
    #     except Exception as ee:
    #         print(ee)


    # 文件保存的临时目录
    temp_cache_directory = settings.STATICFILES_DIRS[0] + 'cache/' + analysis_id + '/'
    temp_cache_directory = temp_cache_directory.replace("\\", "/")
    #todo从临时目录中读取上传的文件，并做数据生成，生成完成后将数据写入到文件中
    if not os.path.exists(temp_cache_directory + 'result.zip'):
        try:
            gene_filepath = temp_cache_directory + 'gene_data.txt'
            info_filepath = temp_cache_directory + 'sampleinfo_data.txt'
            # print(gene_filepath)
            # print(info_filepath)
            r_script_path = settings.STATICFILES_DIRS[0] + 'script/workflow.R'
            command = f'Rscript "{r_script_path}" "{gene_filepath}" "{info_filepath}" "{temp_cache_directory}"'
            os.system(command)
            if os.path.exists(temp_cache_directory + 'result.zip'):
                print('Task run complete!' + temp_cache_directory)
                if re.match('^.*?@.*', to_mail):
                    subject = 'Task run complete!'
                    content = '<a style="font-size: 28px; font-weight: 700; text-align: center;" href="' + endpoint + '/timecourse/analysis/omics_data_analysis_page?analysisId=' + analysis_id + '">--》Click to see the results《--</a>'
                    from_mail = 'Task run notification<' + settings.EMAIL_HOST_USER + '>'
                    print('发送邮件' + content )
                    try:
                        send_mail(subject, message=None, from_email=from_mail, recipient_list=[to_mail], fail_silently=False,
                                  html_message=content)
                        print('Sending an email succeeded！')
                    except Exception as ee:
                        print('发送邮件异常')
                        print(ee)
            # 如果跑脚本没生成就认为是失败
            if not os.path.exists(temp_cache_directory + 'result.zip'):
                if not os.path.exists(temp_cache_directory + 'rc'):
                    with open(temp_cache_directory + 'rc', 'wb') as file:
                        file.write(os.system(command %os.getcwd()))
        except Exception as ee:
            print(ee)
            # 如果上面报错，不管什么错就认为失败
            if not os.path.exists(temp_cache_directory + 'rc'):
                with open(temp_cache_directory + 'rc', 'wb') as file:
                    file.write('2')


def download_result(request):
    analysisId = request.GET.get('analysisId')
    if analysisId == 'example' or analysisId == 'example2':
        basedir = settings.STATICFILES_DIRS[0] + 'exampledata/'
        file_path = basedir + analysisId + '/result.zip'
        file_path = file_path.replace("\\", "/")
    else:
        base_dir = settings.STATICFILES_DIRS[0] + 'cache/'
        file_path = base_dir + analysisId + '/result.zip'
        file_path = file_path.replace("\\", "/")
    try:
        response = StreamingHttpResponse(open(file_path, 'rb'))
        response['content_type'] = "application/octet-stream"
        response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(file_path)
        return response
    except Exception:
        print('路径不存在' + file_path)
        raise Http404

def md5_convert(string):
    import hashlib
    m = hashlib.md5()
    m.update(string.encode())
    return m.hexdigest()

def format_parameters(parameters_str):
    params = parameters_str.split('; ')
    formatted_params = []
    for param in params:
        key, value = param.split('=')
        formatted_value = f'{float(value):.4f}'
        formatted_params.append(f'{key} = {formatted_value}')
    return '; '.join(formatted_params)
