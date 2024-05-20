from django.shortcuts import render
from timecourse_app.models import Sampleinfo
from django.http import HttpResponse, FileResponse,Http404, StreamingHttpResponse
# from django.db.models import Q
from django.conf import settings

import os
import pandas as pd

def get_download_page(request):
    basedir = settings.STATICFILES_DIRS[0] + 'download/'
    file_path = basedir + 'SampleInfo_mysql20230624.txt'
    # file_path = os.getcwd() + '\\timecourse_app_static\\download\\SampleInfo_mysql20230624.txt'
    # file_path=file_path.replace("\\","/")
    # print(file_path)
    if os.path.exists(file_path):
        # print(os.path.getsize(file_path))
        filesize=float('%.1f' % (os.path.getsize(file_path)/1024/1024))
        # print(filesize)
    else:
        filesize = None
    return render(request, 'download.html', {
        'filesize': filesize,
    })

def timeseries_datasets(request):
    basedir = settings.STATICFILES_DIRS[0] + 'download/'
    file_path = basedir + 'SampleInfo_mysql20230624.txt'
    # file_path = (os.getcwd() + '\\timecourse_app_static\\download\\SampleInfo_mysql20230624.txt').replace('\\', '/')
    if os.path.exists(file_path):
        try:
            response = StreamingHttpResponse(open(file_path, 'rb'))
            response['content_type'] = "application/octet-stream"
            response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(file_path)
            return response
        except Exception:
            raise Http404
    else:
        pd.DataFrame(
            Sampleinfo.objects.all().values_list('studyid', 'time', 'unit','tissue_cell_type','organism','data_source'),
            columns=['studyid', 'time', 'unit','tissue_cell_type','organism','data_source']).fillna(
            'null').rename(columns={'studyid': 'TCDid'}).to_csv(file_path, index=False)
        try:
            response = StreamingHttpResponse(open(file_path, 'rb'))
            response['content_type'] = "application/octet-stream"
            response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(file_path)
            return response
        except Exception:
            raise Http404



def cluster_gene_download(request):
    basedir = settings.STATICFILES_DIRS[0] + 'download/'
    file_path = basedir + 'merged_BrowseGene.csv'
    # file_path = (os.getcwd() + '\\timecourse_app_static\\download\\merged_BrowseGene.csv').replace('\\', '/')
    try:
        response = StreamingHttpResponse(open(file_path, 'rb'))
        response['content_type'] = "application/octet-stream"
        response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(file_path)
        return response
    except Exception:
        raise Http404

def model_info_download(request):
    basedir = settings.STATICFILES_DIRS[0] + 'download/'
    file_path = basedir + 'merged_ModelDF.csv'
    # file_path = (os.getcwd() + '\\timecourse_app_static\\download\\merged_ModelDF.csv').replace('\\', '/')
    try:
        response = StreamingHttpResponse(open(file_path, 'rb'))
        response['content_type'] = "application/octet-stream"
        response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(file_path)
        return response
    except Exception:
        raise Http404

def Expression_download(request):
    basedir = settings.STATICFILES_DIRS[0] + 'download/'
    file_path = basedir + 'all_expression_data.zip'
    # file_path = (os.getcwd() + '\\timecourse_app_static\\download\\merged_ModelDF.csv').replace('\\', '/')
    try:
        response = StreamingHttpResponse(open(file_path, 'rb'))
        response['content_type'] = "application/octet-stream"
        response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(file_path)
        return response
    except Exception:
        raise Http404

def readme(request):
    basedir = settings.STATICFILES_DIRS[0] + 'download/'
    file_path = basedir + 'README.txt'
    # file_path = (
    #             os.getcwd() + '\\timecourse_app_static\\download\\README.txt').replace(
    #     '\\', '/')
    try:
        response = StreamingHttpResponse(open(file_path, 'rb'))
        response['content_type'] = "application/octet-stream"
        response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(file_path)
        return response
    except Exception:
        raise Http404