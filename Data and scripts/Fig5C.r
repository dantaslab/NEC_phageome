library(readxl)

Data = read_excel("Metabolism_early_sum.xlsx")
# rownames(Data) = Data$Metabolism
# Data$Metabolism = NULL
DataT = t(Data[,c(2,3)])
colnames(DataT) = Data$Metabolism
DataT[,-12]

pdf(file="250423_bar_AMG_countu.pdf",
    width=5, height=5)

barplot(DataT[,-12], 
        col= c("#f28a5f", "#3291a4"),
        border=c("#f28a5f", "#3291a4"), 
        legend.text = TRUE,
        horiz = TRUE,
        font.axis=2, 
        space=0.5, 
        xlab="AMG count")
dev.off()
