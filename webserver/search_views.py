import re

from django.shortcuts import render
from django.http import HttpResponse
from django.db.models import Q, F
from django.core.paginator import Paginator
from djangoProject.settings import STATIC_ROOT
from timecourse_app.models import Sampleinfo, EvaIndicator, ModelDF, GO,KEGG
import numpy as np
import pandas as pd
import json, os, re , time
import sys
from collections import Counter

# Create your views here.

def hello_world(request):
    return HttpResponse("Hello World!")


def search(request):
    return render(request, 'search_bak.html')


def result(request):
    key = request.GET['query_content']
    page = request.GET.get('page')
    if not key :
        message = u'请输入搜索内容'
        return render(request, 'result.html',{'message':message})
    else :
        result_information = Sampleinfo.objects.filter(Q(studyid__icontains=key)|Q(accession=key)|Q(sample__icontains=key)|Q(related_factor__icontains=key)|Q(organism=key)| \
                                                      Q(tissue_cell_type__icontains=key).order_by('studyid'))
        study_search_str = len(result_information)
        page_search_list, curr_page, prev_page, next_page, begin_page, end_page = divide_page(page, result_information)
        return render(request, 'result.html', {
                'result_informarion':page_search_list,
                'page_num': range(begin_page,end_page + 1),
                'curr_page': curr_page,
                'next_page': next_page,
                'prev_page': prev_page,
                'query_content': key,
                'study_search_str': study_search_str,
            })


def detail(request, studyid):
    all_study = Sampleinfo.objects.all()
    curr_study = None
    for study in all_study:
        if study.studyid == studyid:
            curr_study = Sampleinfo.objects.filter(studyid=studyid).values() \
            .annotate(organism=F('organism__organism')).annotate(data_source=F('data_source__datasource')) \
            .annotate(library_strategy=F('library_strategy__library_strategy'))[0]

            EvaIndicator_data = pd.DataFrame(EvaIndicator.objects.filter(studyid__studyid=studyid).
                              values_list('model','rmse','rsq','aic','cluster'),
                              columns=['model','rmse','rsq','aic','cluster']).\
                              round({'rmse': 3, 'rsq': 3, 'aic': 3})


            ModelDF_data = pd.DataFrame(ModelDF.objects.filter(studyid__studyid=studyid).
                              values_list('cluster','model','parameters','exp'),
                              columns=['cluster','model','parameters','exp'])
            break
    return render(request, 'detail_bak.html', {
        'curr_study': curr_study,
        })

# def get_bubble(request):
#     studyid = request.POST.get('studyid')
#     ont = request.POST.get('ont')
#     library_strategy = Sampleinfo.objects.filter(studyid=studyid).values_list('library_strategy__library_strategy', flat=True)[0]
#     GO_data = pd.DataFrame(GO.objects.filter(Q(studyid__studyid=studyid)).
#                 values('GO_ID','Description','ONT','GeneRatio','BgRatio','pvalue','p_adjust','qvalue','geneID','Count','cluster'),
#                 columns=['GO_ID','Description','ONT','GeneRatio','BgRatio','pvalue','p_adjust','qvalue','geneID','Count','cluster']).\
#          round({'qvalue': 3,'pvalue': 3,'p_adjust': 3}).to_dict('records')
#
#     KEGG_data = pd.DataFrame(KEGG.objects.filter(Q(studyid__studyid=studyid)).
#                 values('KEGG_ID','Description','GeneRatio','BgRatio','pvalue','p_adjust','qvalue','geneID','Count','cluster'),
#                 columns=['KEGG_ID','Description','GeneRatio','BgRatio','pvalue','p_adjust','qvalue','geneID','Count','cluster']).\
#          round({'qvalue': 3,'pvalue': 3,'p_adjust': 3}).to_dict('records')
#
#     if ont != 'kegg':
#         if ont == 'undefined':
#             ont = 'bp'
#             enrich_result = GO_data[GO_data['ont'] == ont.upper()]
#             enrich_obj = bubble_obj_fun(enrich_result)
#         else:
#             enrich_obj = None
#
#     else:
#         if os.path.exists(KEGG_path):
#             enrich_obj = bubble_obj_fun(KEGG_data)
#         else:
#             enrich_obj = None
#     return HttpResponse(json.dumps(enrich_obj))
#
# def bubble_obj_fun(enrich_result):
#     enrich_result['nlog10pvalue'] = -np.log10(enrich_result['pvalue'])
#     enrich_result_top10 = enrich_result.iloc[0:10, ].sort_values(by='nlog10pvalue')
#     x = enrich_result_top10['nlog10pvalue'].tolist()
#     y = enrich_result_top10['description'].tolist()
#     count = enrich_result_top10['count'].tolist()
#     if max(enrich_result_top10['count']) > 30:
#         size = (enrich_result_top10['count'] * 30/max(enrich_result_top10['count'])).tolist()
#     elif max(enrich_result_top10['count']) < 5:
#         size = [(_ + 10) for _ in count]
#     else:
#         size = count
#     enrich_obj = [{
#         'mode': 'markers',
#         'name': 'industry1',
#         'type': 'scatter',
#         'x': x,
#         'y': y,
#         'marker': {
#             'line': {
#                 'width': 1.3
#             },
#             'color': 'rgba(55, 128, 191, 1.0)',
#             'symbol': 'dot',
#             'opacity': 0.8,
#             'size': size
#         },
#         'text': ['size: ' + str(c) for c in count],
#         'textfont': {
#             'color': '#4D5663'
#         }
#     }]
#     return(enrich_obj)

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