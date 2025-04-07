#Linux 
#设置环境
conda activate tangram-env

#安装相关包
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install scipy
pip install matplotlib
pip install seaborn
pip install jupyterlab
pip install numpy
pip install pandas==1.2.0
pip install scanpy
pip install scikit-learn
pip install tqdm
pip install tangram-sc

#Python
#加载包
import scanpy as sc
import tangram as tg
import gc
import psutil

#加载数据
adata_st = sc.read_h5ad("st_data_3_10.h5ad")
adata_sc = sc.read_h5ad("sc_data_V4.h5ad")

# 检查聚类标签是否存在
if 'integrated_cluster' not in adata_sc.obs:
    raise ValueError("adata_sc 中缺少 'integrated_cluster' 字段，请先进行聚类分析")

# 寻找共同基因
tg.pp_adatas(adata_sc, adata_st, genes=None)

# Cluster Level 映射
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_st,
    mode='clusters',
    cluster_label='integrated_cluster',
    device='cuda' if torch.cuda.is_available() else 'cpu'  # 自动检测GPU
)

# 释放内存
gc.collect()

# 整合并投影基因表达
ad_ge = tg.project_genes(
    ad_map,
    adata_sc,  
    cluster_label='integrated_cluster'
)

# 保存结果
ad_map.write_h5ad("cell_to_space_mapping.h5ad")  # 保存 ad_map 以便复现
ad_ge.write_h5ad("spatial_gene_expression.h5ad")
print("整合完成，结果已保存至 cell_to_space_mapping.h5ad 和 spatial_gene_expression.h5ad"
