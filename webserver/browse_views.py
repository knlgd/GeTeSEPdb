import json
import os
import time
import pickle
from django.http import HttpResponse
from django.shortcuts import render
from timecourse_app.models import Sampleinfo, clustergene,ModelDF
import pandas as pd
import re
from django.conf import settings
def get_browse_page(request):
    select = request.GET.get('select')
    browserType = request.GET.get('browserType')
    print(select)
    # pkl_file = os.getcwd() + '/timecourse_app_static/data/browse_showdata.pkl'
    basedir = settings.STATICFILES_DIRS[0] + 'data/'
    pkl_file = basedir + 'browse_showdata.pkl'
    allsampleinfo = pd.DataFrame(Sampleinfo.objects.all().values('studyid', 'accession', 'time', 'unit', 'related_factor', 'tissue_cell_type',
                                    'experiment_type', 'organism', 'strategy')).fillna('null').drop_duplicates()
    allsampleinfo_data=allsampleinfo.to_dict('records')

    showdf = pd.read_pickle(pkl_file)
    showdf = showdf.fillna('')
    showdata = showdf.to_dict('records')
    if select is None:
        if browserType != 'browser_parttern':
            dataset_tissue, dataset_tissue_list = attr_maker(allsampleinfo, 'tissue_cell_type', prefix='browser_dataset')
            dataset_organism, dataset_organism_list = attr_maker(allsampleinfo, 'organism', prefix='browser_dataset')
            dataset_library_strategy, dataset_library_strategy_list = attr_maker(allsampleinfo, 'strategy', prefix='browser_dataset')
            dataset_experiment_type, dataset_experiment_type_list = attr_maker(allsampleinfo, 'experiment_type', prefix='browser_dataset')
            show_tissue = {}
            show_tissue_list = {}
            show_organism = {}
            show_organism_list = {}
            show_library_strategy = {}
            show_library_strategy_list = {}
            show_model_type = {}
            show_model_type_list = {}
            show_experiment_type = {}
            show_experiment_type_list = {}
        else:
            show_tissue, show_tissue_list = attr_maker(showdf, 'tissue_cell_type', prefix='browser_parttern')
            show_organism, show_organism_list = attr_maker(showdf, 'organism', prefix='browser_parttern')
            show_library_strategy, show_library_strategy_list = attr_maker(showdf, 'strategy', prefix='browser_parttern')
            show_model_type, show_model_type_list = attr_maker(showdf, 'model', prefix='browser_parttern')
            show_experiment_type, show_experiment_type_list = attr_maker(showdf, 'experiment_type', prefix='browser_parttern')
            dataset_tissue = {}
            dataset_tissue_list = {}
            dataset_organism = {}
            dataset_organism_list = {}
            dataset_library_strategy = {}
            dataset_library_strategy_list = {}
            dataset_experiment_type = {}
            dataset_experiment_type_list = {}

        return render(request, 'BROWSE.html', {
            'showdata': showdata,
            'show_tissue': show_tissue,
            'show_organism': show_organism,
            'show_library_strategy': show_library_strategy,
            'show_model_type': show_model_type,
            'show_experiment_type': show_experiment_type,
            'show_tissue_list': show_tissue_list,
            'show_organism_list': show_organism_list,
            'show_library_strategy_list': show_library_strategy_list,
            'show_model_type_list': show_model_type_list,
            'show_experiment_type_list': show_experiment_type_list,
            'allsampleinfo_data': allsampleinfo_data,
            'dataset_tissue': dataset_tissue,
            'dataset_organism': dataset_organism,
            'dataset_library_strategy': dataset_library_strategy,
            'dataset_experiment_type': dataset_experiment_type,
            'dataset_tissue_list': dataset_tissue_list,
            'dataset_organism_list': dataset_organism_list,
            'dataset_library_strategy_list': dataset_library_strategy_list,
            'dataset_experiment_type_list': dataset_experiment_type_list,
            'data_source': pd.DataFrame(['GEO', 'ArrayExpress'], columns=['data_source']).to_dict('records'),
        })

    # select = pd.DataFrame(select.split(','), columns=['select'])
    # print(select)
    # select['data_type'], select['data_name'] = select['select'].str.split(':', 1).str
    # select = select.drop_duplicates()
    # select['data_name'] = select['data_name'].apply(lambda x: x.replace('\'', ''))
    # showstr = select[['data_type', 'data_name']].to_dict('records')
    #
    # select_dict = {}
    # for t in list(select['data_type'].drop_duplicates()):
    #     select_dict[t] = select[select['data_type'] == t]['data_name'].tolist()[0]
    if browserType == 'browser_dataset':
        select = pd.DataFrame(select.split(','), columns=['select'])
        select['data_type'], select['data_name'] = select['select'].str.split(':', 1).str
        select['data_name'] = select['data_name'].apply(lambda x: x.replace('\'', ''))
        showstr = select[['data_type', 'data_name']].to_dict('records')
        select_dict = {}
        for _, row in select.iterrows():
            data_type = row['data_type']
            data_name = row['data_name']
            if data_type in select_dict:
                select_dict[data_type].append(data_name)
            else:
                select_dict[data_type] = [data_name]

        filtered_data = allsampleinfo.copy()
        for field, options in select_dict.items():
            filtered_data = filtered_data[filtered_data[field].isin(options)]
        allsampleinfo_data=filtered_data.fillna('null').to_dict('records')

        # allsampleinfo_data = allsampleinfo[allsampleinfo[list(select_dict)].eq(pd.Series(select_dict)).all(axis=1)].fillna('null').to_dict('records')
        dataset_tissue, dataset_tissue_list = attr_maker(allsampleinfo, 'tissue_cell_type', prefix='browser_dataset')
        dataset_organism, dataset_organism_list = attr_maker(allsampleinfo, 'organism', prefix='browser_dataset')
        dataset_library_strategy, dataset_library_strategy_list = attr_maker(allsampleinfo, 'strategy',
                                                                             prefix='browser_dataset')
        dataset_experiment_type, dataset_experiment_type_list = attr_maker(allsampleinfo, 'experiment_type',
                                                                           prefix='browser_dataset')
        show_tissue = {}
        show_tissue_list = {}
        show_organism = {}
        show_organism_list = {}
        show_library_strategy = {}
        show_library_strategy_list = {}
        show_model_type = {}
        show_model_type_list = {}
        show_experiment_type = {}
        show_experiment_type_list = {}
    else:
        # print(select)
        select_dict = {}
        select_items = select.split(',')
        select = pd.DataFrame(select.split(','), columns=['select'])
        select['data_type'], select['data_name'] = select['select'].str.split(':', 1).str
        select['data_name'] = select['data_name'].apply(lambda x: x.replace('\'', ''))
        showstr = select[['data_type', 'data_name']].to_dict('records')

        for item in select_items:
            key, value = item.split(':')
            key = key.strip()
            value = value.strip()
            if key in select_dict:
                select_dict[key].append(value)
            else:
                select_dict[key] = [value]

        cleaned_select_dict = {key: [value.strip("'") for value in values] for key, values in select_dict.items()}
        filtered_data = showdf.copy()
        for field, values in cleaned_select_dict.items():
            condition = filtered_data[field].isin(values)
            filtered_data = filtered_data[condition]
        showdata = filtered_data.fillna('null').to_dict('records')

        # showdata = showdf[showdf[list(select_dict)].eq(pd.Series(select_dict)).all(axis=1)].fillna('null').to_dict(
        #     'records')
        show_tissue, show_tissue_list = attr_maker(showdf, 'tissue_cell_type', prefix='browser_parttern')
        show_organism, show_organism_list = attr_maker(showdf, 'organism', prefix='browser_parttern')
        show_library_strategy, show_library_strategy_list = attr_maker(showdf, 'strategy', prefix='browser_parttern')
        show_model_type, show_model_type_list = attr_maker(showdf, 'model', prefix='browser_parttern')
        show_experiment_type, show_experiment_type_list = attr_maker(showdf, 'experiment_type', prefix='browser_parttern')
        dataset_tissue = {}
        dataset_tissue_list = {}
        dataset_organism = {}
        dataset_organism_list = {}
        dataset_library_strategy = {}
        dataset_library_strategy_list = {}
        dataset_experiment_type = {}
        dataset_experiment_type_list = {}

    return render(request, 'BROWSE.html', {
        'showdata': showdata,
        'showstr': showstr,
        'show_tissue': show_tissue,
        'show_organism': show_organism,
        'show_library_strategy': show_library_strategy,
        'show_model_type': show_model_type,
        'show_experiment_type': show_experiment_type,
        'show_tissue_list': show_tissue_list,
        'show_organism_list': show_organism_list,
        'show_library_strategy_list': show_library_strategy_list,
        'show_model_type_list': show_model_type_list,
        'show_experiment_type_list': show_experiment_type_list,
        'allsampleinfo_data': allsampleinfo_data,
        'dataset_tissue': dataset_tissue,
        'dataset_organism': dataset_organism,
        'dataset_library_strategy': dataset_library_strategy,
        'dataset_experiment_type': dataset_experiment_type,
        'dataset_tissue_list': dataset_tissue_list,
        'dataset_organism_list': dataset_organism_list,
        'dataset_library_strategy_list': dataset_library_strategy_list,
        'dataset_experiment_type_list': dataset_experiment_type_list,
        'data_source': pd.DataFrame(['GEO', 'ArrayExpress'], columns=['data_source']).to_dict('records'),
    })


def gene_table(request):
    studyid = request.POST.get('studyid')
    cluster = request.POST.get('cluster')
    # studyid='TCD0001'
    # cluster='cluster2'

    # start_time = time.time()
    #time=2.84344220161438
    genedf = pd.DataFrame(clustergene.objects.filter(studyid=studyid, cluster=cluster).values_list('studyid', 'gene', 'cluster','model', 'cv', 'mean'),columns=['studyid', 'gene', 'cluster', 'model', 'cv', 'mean'])
    # basedir = settings.STATICFILES_DIRS[0] + 'data/'
    # pkl_file = basedir + 'search_showdata.pkl'
    # genedf = pd.read_pickle(pkl_file)
    # genedf = genedf[(genedf['studyid'] == studyid) & (genedf['cluster'] == cluster)]
    # genedf = genedf.fillna('')
    # 替换0.0为NaN
    genedf = genedf[genedf['model'] != 'exponential_decay']
    genedf = genedf.replace({'Logis': 'Logistic', 'Sin': 'Trigonometric', 'exponential_growth': 'Exponential'})
    genedf['cv'] = genedf['cv'].apply(lambda x: round(x, 4) if x is not None else None)
    genedf['mean'] = genedf['mean'].apply(lambda x: round(x, 4) if x is not None else None)
    genedf['cv'] = genedf['cv'].replace(0.0, '')
    genedf['mean'] = genedf['mean'].replace(0.0, '')
    genedf.fillna('Null', inplace=True)
    gene_data=genedf.to_dict('records')
    # print(time.time() - start_time)
    # print(gene_data)
    return HttpResponse(json.dumps(gene_data))



def attr_maker(obj, attr, prefix=''):
    this_attr = pd.DataFrame(obj[attr], columns=[attr])[attr].drop_duplicates().dropna().sort_values()
    this_attr = this_attr[(this_attr.str.contains('\\|') == False) & (this_attr.str.contains('\\_') == False)].append(this_attr).drop_duplicates()

    show_selectmenu = pd.DataFrame()
    style_content = "style='padding: 10px 500px 10px 50px;margin: -10px -500px -10px -50px;'"
    show_selectmenu['name'] = [
        ("<span onclick='select_data(this)' onkeydown='keyboard_select_data(this)' title='%s' id='" + prefix + "_" + attr + "_%s' data-name='%s' data-browser-type='" + prefix + "' data-type='" + attr + "'" + style_content + ">%s</span>") % (s, s.replace(",", "_"), s, s)
        for s in this_attr
    ]
    show_selectmenu['id'] = range(show_selectmenu.shape[0])
    show_selectmenu = show_selectmenu.to_dict('records')

    show_pagination = []
    for s in this_attr:
        show_pagination.append([("<span title='%s' onclick='select_data(this)' data-browser-type='" + prefix + "'  class='browser-db-column list-group-item' id='" + prefix + "_" + attr + "_%s' data-name='%s' data-type='" + attr + "'>%s</span>") % (s, s.replace(",", "_"), s, s)])

    return show_selectmenu, show_pagination


# ###############make pkl file
# basedir = settings.STATICFILES_DIRS[0] + 'data/'
# pkl_file = basedir + 'browse_showdata.pkl'
# # # if os.path.isfile(pkl_file):
# # #     # 如果pkl文件已存在，则读取文件
#     showdf = pd.read_pickle(pkl_file)
# # # else:
# allsampleinfo = pd.DataFrame(
#         Sampleinfo.objects.all().values('studyid','tissue_cell_type','experiment_type','organism','strategy')).drop_duplicates()
# genedf = pd.DataFrame(clustergene.objects.all().values_list('studyid', 'cluster','model'),
#                           columns=['studyid', 'cluster','model']).drop_duplicates()
# genedf = genedf[genedf['model'] != 'exponential_growth']
# #
# # #genedf = genedf.replace( {'trigonometric': 'Trigonometric'})
# genedf = genedf.replace( {'Logis': 'Logistic', 'Sin': 'Trigonometric','constant': 'Constant','exponential_decay': 'Exponential'})
# showdf = genedf.merge(allsampleinfo, on='studyid', how='left')
# showdf = showdf[showdf['organism'].notna()]
# showdf = showdf[showdf['strategy'].notna()]
#
# showdf.sort_values(by=['studyid', 'cluster'], inplace=True)
#
# # 保存DataFrame为pkl文件
# with open(pkl_file, 'wb') as f:
#     pickle.dump(showdf, f)
# # #
# # # #################
# allsampleinfo = pd.DataFrame(
#         Sampleinfo.objects.all().values('studyid', 'accession', 'time', 'unit','related_factor','tissue_cell_type','experiment_type','organism','strategy')).drop_duplicates()
# genedf = pd.DataFrame(clustergene.objects.all().values_list('studyid', 'cluster', 'model','gene','cv','mean'),
#                           columns=['studyid', 'cluster', 'model','gene','cv','mean']).drop_duplicates()
# genedf = genedf[genedf['model'] != 'exponential_growth']
# genedf = genedf.replace( {'Logis': 'Logistic', 'Sin': 'Trigonometric','constant': 'Constant','exponential_decay': 'Exponential'})
# showdf = genedf.merge(allsampleinfo, on='studyid', how='left')
# showdf = showdf[showdf['organism'].notna()]
# showdf = showdf[showdf['strategy'].notna()]
# showdf.fillna('Null', inplace=True)
# # # rows_with_empty_values = showdf[showdf['strategy'].isnull()]
# # # rows_with_empty_values = showdf[showdf['organism'].isnull()]
# # # print(rows_with_empty_values)
# basedir = settings.STATICFILES_DIRS[0] + 'data/'
# pkl_file = basedir + 'search_showdata.pkl'
# with open(pkl_file, 'wb') as f:
#     pickle.dump(showdf, f)


#统计每个model的cluster、gene和studyid数目
# model_counts = showdf.groupby('organism').agg({'cluster': 'size', 'gene': 'size', 'studyid': 'nunique'}).reset_index()
# totals = model_counts.sum(numeric_only=True)
# totals['organism'] = 'Total'
# model_counts = model_counts.append(totals, ignore_index=True)
# # 重命名统计结果的列名
# model_counts.columns = ['organism', 'cluster总数', 'gene总数', 'studyid总数']
#
# #nunique
# model_counts = showdf.groupby('model').agg({'cluster': 'size', 'gene': 'size', 'studyid': 'nunique'}).reset_index()
# totals = model_counts.sum(numeric_only=True)
# totals['model'] = 'Total'
# model_counts = model_counts.append(totals, ignore_index=True)
# # 重命名统计结果的列名
# model_counts.columns = ['model', 'cluster总数', 'gene总数', 'studyid总数']

############################统计图
# 统计每个studyid的gene数目
# basedir = settings.STATICFILES_DIRS[0] + 'data/'
# pkl_file = basedir + 'gene_count.pkl'
# studyid_counts = showdf.groupby('studyid')['gene'].nunique().reset_index()
# # 重命名统计结果的列名
# studyid_counts.columns = ['studyid', 'n.gene']
# with open(pkl_file, 'wb') as f:
#     pickle.dump(studyid_counts, f)

#####################################
# basedir = settings.STATICFILES_DIRS[0] + 'data/'
# pkl_file = basedir + 'gene_count.pkl'
# # if os.path.isfile(pkl_file):
# #     # 如果pkl文件已存在，则读取文件
# showdf = pd.read_pickle(pkl_file)
# showdf = showdf[showdf['n.gene'].notna()].drop_duplicates()
# sum_n_gene = showdf['n.gene'].sum()


################统计基因出现次数
# 统计organism为Homo_sapiens，模型为Constant的基因出现的次数
# constant_homo_genes_count = showdf[(showdf['organism'] == 'Homo_sapiens') & (showdf['model'].str.contains('Constant'))]['gene'].value_counts()
# constant_homo_genes_count.to_csv('E:/研究生阶段工作内容/时序分析/housekeeping分析数据/constant_homo_genes_count.csv')
# # 找出出现次数大于等于40次以上的基因
# genes_above_40 = constant_homo_genes_count[constant_homo_genes_count >= 40]
# genes_above_40.to_csv('E:/研究生阶段工作内容/时序分析/housekeeping分析数据/genes_above_40.csv')