from django.shortcuts import render
from django.db.models import Q

def get_help_page(request):
    return render(request, 'help.html')