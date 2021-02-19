
library(clusterProfiler)

tmp <- read.csv('./激酶列表.csv')
eg = bitr(tmp$UPA, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
eg2 = bitr(tmp$磷酸化, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

kk <- enrichKEGG(gene         = eg$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
head(kk)

kk2 <- enrichKEGG(gene         = eg2$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
library(enrichplot)
library(ggplot2)
p1 <- dotplot(kk, showCategory=10,color = "pvalue") + ggtitle("dotplot for IPA") + scale_color_continuous(low='#FF332C', high='#0074B6')

p1
p2 <- dotplot(kk2, showCategory=10,color = "pvalue") + ggtitle("dotplot for protein")+ scale_color_continuous(low='#FF332C', high='#0074B6')

p2
pdf('kkk.pdf',width = 18,height =8)
plot_grid(p1, p2, ncol=2)
dev.off()

upsetplot(kk)


upa <- read.csv('./UPA-link-2.0.csv',header = F)

upa <- upa[,c(1,2)]

upa$V1 <- paste0(substring(upa$V1, 1, 1),tolower(substring(upa$V1, 2)))
upa$V2 <- paste0(substring(upa$V2, 1, 1),tolower(substring(upa$V2, 2)))

table(duplicated(upa))

write.csv(upa,file = 'newupa3.csv',quote = F)


