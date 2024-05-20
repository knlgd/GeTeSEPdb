# Create your models here.
from django.db import models

#Create your models here.
class Sampleinfo(models.Model):
    # StudyID
    studyid = models.CharField(primary_key = True,max_length=254)
    # Accession
    accession = models.CharField(max_length=254, null=True)
    sample = models.TextField(max_length=254, null=True)
    time = models.TextField(max_length=254, null=True)
    phototime = models.TextField(max_length=254, null=True)
    unit = models.CharField(max_length=254, null=True)
    related_factor = models.TextField(null=True)
    drug = models.TextField(null=True)
    drug_concentration = models.TextField(null=True)
    perturbation = models.TextField(null=True)
    tissue_cell_type = models.TextField(null=True)
    organism = models.TextField(null=True)
    data_source = models.TextField(null=True)
    strategy = models.CharField(max_length=254,null=True)
    pmid = models.TextField(null=True)
    disease = models.TextField(null=True)
    other = models.TextField(null=True)
    experiment_type = models.TextField(null=True)
    Data_prepro_info = models.TextField(null=True)

class EvaIndicator(models.Model):
    id = models.AutoField(primary_key=True, db_index=True)
    studyid = models.CharField(max_length=254)
    model = models.CharField(max_length=254)
    rmse = models.FloatField(null=True)
    rsq = models.FloatField(null=True)
    aic = models.FloatField(null=True)
    cluster = models.CharField(max_length=254,null=True)

class ModelDF(models.Model):
    id = models.AutoField(primary_key=True, db_index=True)
    cluster = models.CharField(max_length=254,null=True)
    model = models.CharField(max_length=254,null=True)
    parameters = models.TextField(max_length=254,null=True)
    exp = models.TextField(max_length=254,null=True)
    studyid = models.CharField(max_length=254)

class clustergene(models.Model):
    id = models.AutoField(primary_key=True, db_index=True)
    studyid = models.CharField(max_length=254)
    cluster = models.CharField(max_length=254,null=True)
    model = models.CharField(max_length=254,null=True)
    gene = models.TextField(max_length=254,null=True)
    cv= models.FloatField(max_length=254,null=True)
    mean = models.FloatField(max_length=254,null=True)
    organism = models.TextField(null=True)
    experiment_type = models.TextField(null=True)
    tissue_cell_type = models.TextField(null=True)
    strategy = models.CharField(max_length=254, null=True)



class GO(models.Model):
    id = models.AutoField(primary_key=True, db_index=True)
    GO_ID = models.TextField(max_length=254)
    Description = models.CharField(max_length=254,null=True)
    ONT = models.TextField(max_length=254,null=True)
    GeneRatio = models.TextField(max_length=254,null=True)
    BgRatio = models.TextField(max_length=254, null=True)
    pvalue = models.FloatField(max_length=254,null=True)
    p_adjust = models.FloatField(null=True)
    qvalue = models.FloatField(null=True)
    geneID = models.TextField(max_length=254)
    Count = models.IntegerField(null=True, db_index=True)
    cluster = models.CharField(max_length=254,null=True)
    studyid = models.CharField(max_length=254)

class KEGG(models.Model):
    id = models.AutoField(primary_key=True, db_index=True)
    KEGG_ID = models.TextField(max_length=254)
    Description = models.CharField(max_length=254, null=True)
    GeneRatio = models.TextField(max_length=254, null=True)
    BgRatio = models.TextField(max_length=254, null=True)
    pvalue = models.FloatField(null=True)
    p_adjust = models.FloatField(null=True)
    qvalue = models.FloatField(null=True)
    geneID = models.TextField(max_length=254)
    Count = models.IntegerField(null=True, db_index=True)
    cluster = models.CharField(max_length=254,null=True)
    studyid = models.CharField(max_length=254)


class geneinfo(models.Model):
    id = models.AutoField(primary_key=True, db_index=True)
    geneid = models.IntegerField(null=True, db_index=True)
    gene = models.TextField(max_length=254, null=True)
    synonyms = models.TextField(max_length=254, null=True)
    description = models.TextField(max_length=254, null=True)
    ensembl = models.TextField(max_length=254, null=True)
    specie = models.TextField(max_length=254, null=True)
