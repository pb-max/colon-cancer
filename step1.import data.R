rm(list = ls())#ħ�ò�����һ�����~
options(stringsAsFactors = F)
library(data.table)
library(R.utils)
pca1 <- fread("./GSM4773521_PCa1_gene_counts_matrix.txt.gz",data.table = F)
#�����warning message ��Ҫ�£�����Ϊ��ȡ��ʱ�������������һ�У����Զ����һ�С�V1���������ɾ������һ����������������
pca1[1:4,1:4]
pca2 <- fread('./GSM4773522_PCa2_gene_counts_matrix.txt.gz',data.table = F)
pca2[1:4,1:4]
d1 <- pca1[ ,-1]
rownames(d1) <- pca1[ ,1]

#����SeuratObject
#�ֱ𴴽�SeuratObject
scRNA1 <- CreateSeuratObject(counts = d1,#������һ���������ʿ�
                             min.cells = 3,#�ڲ���������ϸ���б���
                             min.features = 200,#������������200
                             project = "pca1")
#24576 features & 2896 samples

d2 = pca2[ ,-1]
rownames(d2) <- pca2