#!/usr/bin/env python
# coding: utf-8

# # Annotation with CellAssign

# ## Assigning single-cell RNA-seq data to known cell types

# CellAssign is a probabilistic model that uses prior knowledge of cell-type marker genes to annotate scRNA data into predefined cell types. Unlike other methods for assigning cell types, CellAssign does not require labeled single cell data and only needs to know whether or not each given gene is a marker of each cell type. The original paper and R code are linked below.
# 
# Paper: [Probabilistic cell-type assignment of single-cell RNA-seq for tumor microenvironment profiling, *Nature Methods 2019*](https://www.nature.com/articles/s41592-019-0529-1 )
# 
# Code: https://github.com/Irrationone/cellassign

# This notebook will demonstrate how to use CellAssign on follicular lymphoma and HGSC scRNA data.

# To demonstrate CellAssign, we use the data from the original publication, which we converted into h5ad format. The data are originally available from here:
# 
# https://zenodo.org/record/3372746

# ### import scvi comment out korsi

# In[40]:


# import scvi
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

sc.set_figure_params(figsize=(4, 4))




# In[41]:


from scvi.external import CellAssign
CellAssign.hello_model()


# ## Follicular Lymphoma Data
# 
# Load follicular lymphoma data and marker gene matrix (see Supplementary Table 2 from the original paper).

# In[42]:


adata = sc.read("data/sce_follicular_annotated_final.h5ad")
adata.var_names_make_unique()
adata.obs_names_make_unique()


# In[43]:


# marker_gene_mat = pd.read_csv('data/FL_celltype.csv', index_col=0)
# marker_gene_mat= pd.read_csv("data/follia1/var.csv", usecols =[0])
marker_gene_mat = pd.read_csv('data/follia1/var5.csv', index_col=0)


# ### Create and fit CellAssign model

# The anndata object and cell type marker matrix should contain the same genes, so we index into `adata` to include only the genes from `marker_gene_mat`.

# In[44]:


bdata = adata[:, marker_gene_mat.index].copy()
print(bdata.n_obs)


# Then we setup anndata and initialize a `CellAssign` model. Here we set the `size_factor_key` to "size_factor", which is a column in `bdata.obs`. A size factor may be defined manually as scaled library size (total UMI count) and should not be placed on the log scale, as the model will do this manually. The library size should be computed before any gene subsetting (in this case, technically, a few notebook cells up).
# 
# For example,
# 
# ```python
# lib_size = adata.X.sum(1)
# adata.obs["size_factor"] = lib_size / np.mean(lib_size)
# ```

# In[45]:




from scvi.external import CellAssign


# In[46]:
print(adata.n_obs)
# adata = sc.pp.filter_cells(adata, min_genes=30)

# add some true zeros
# adata.X[adata.X <= 0.3] = 0
# simply compute the number of genes per cell
sc.pp.filter_cells(adata, min_genes=0)
print("adata.obs['n_genes'].min():")
print(adata.obs['n_genes'].min())

print("adata.obs['n_genes'].max():")
print(adata.obs['n_genes'].max())



sc.pp.filter_cells(adata, min_genes=2300)
print("adata.obs['n_genes'].min():")
print(adata.obs['n_genes'].min())


# filter manually
# adata_copy = adata[adata.obs['n_genes'] >= 3]
# adata_copy.obs['n_genes'].min()
# adata.n_obs

# adata.obs['n_genes'].min()

# actually do some filtering
# sc.pp.filter_cells(adata, min_genes=3)
# adata.n_obs
# adata.obs['n_genes'].min()

# sc.pp.filter_cells(adata, min_counts=300)
print(adata.n_obs)
# print(adata)

CellAssign.setup_anndata(adata, size_factor_key="size_factor")


# In[47]:


# from scvi.external import CellAssign
# model = CellAssign(bdata, marker_gene_mat)

model = CellAssign(adata, marker_gene_mat)


# In[48]:


model.train()


# Inspecting the convergence:
