from django.urls import path
from timecourse_app import views, browse_views, adv_search_views, analysis, statistic_views, download_views, help_views,adv_search_function
urlpatterns=[
    path('hello_world', views.hello_world),
    path('home', views.search),
    path('', views.search),
    # path('get_boxplot', views.get_lineplot),
    # path('adv_search', adv_search_views.get_adv_search_page),
    path('get_search_page', adv_search_function.get_search_page),
    path('get_experimenttype', adv_search_function.get_experimenttype),
    path('get_modeltype', adv_search_function.get_modeltype),
    path('get_tissue', adv_search_function.get_tissue),
    path('get_datatype', adv_search_function.get_datatype),
    path('get_gene', adv_search_function.get_gene),
    path('get_datatable',adv_search_function.get_datatable),
    path('get_gene_names', adv_search_function.get_gene_names),
    path('search_adv_select_items', adv_search_function.search_adv_select_items),

    path('plot_trend',adv_search_views.plot_trend),
    path('gene_detail/<str:gene>', adv_search_views.gene_detail),
    # path('get_modeltype', adv_search_views_bak2.get_modeltype),
    # path('get_tissue', adv_search_views_bak2.get_tissue),
    # path('get_datatype', adv_search_views_bak2.get_datatype),

    path('result', views.result),
    path('detail/<str:studyid>',views.detail),
    path('get_bubble',views.get_bubble),
    path('browse', browse_views.get_browse_page),
    path('gene_table',browse_views.gene_table),
    # 回溯页面
    path('analysis/retrieve_page',analysis.analysis_retrieve_page),
    path('analysis/retrieve_result',analysis.retrieve_result),

    # single gene 页面
    path('analysis/single_gene_analysis_page',analysis.single_gene_analysis_page),
    # single gene提交任务接口
    path('analysis/single_gene_analysis_submit',analysis.single_gene_analysis_submit),
    # path('get_analysis',analysis.get_analysis),

    # 分析结果视图页（废弃）
    path('get_analysis_result_page', analysis.get_analysis_result_page),

    # omics data 分析页面
    path('analysis/omics_data_analysis_page',analysis.get_analysis_page),
    # omics data 获取分析结果数据的接口
    path('analysis/omics_data_analysis_result_data', analysis.get_analysis_result_data),
    # omics data 提交任务的接口
    path('analysis/omics_data_analysis_submit', analysis.submit_analysis),

    # 分析结果下载接口
    path('download_result', analysis.download_result),
    # statistic
    path('statistic',statistic_views.get_statistic_page),
    # # download
    path('download',download_views.get_download_page),
    path('download/timeseries_datasets',download_views.timeseries_datasets),
    path('download/cluster_gene_download',download_views.cluster_gene_download),
    path('download/model_info_download',download_views.model_info_download),
    path('download/Expression_download',download_views.Expression_download),
    path('download/readme',download_views.readme),
    # # help
    path('help',help_views.get_help_page),
]

