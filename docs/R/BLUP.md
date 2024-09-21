# BLUP值

最佳线性无偏预测（Best Linear Unbiased Prediction，简称BLUP）可以对多环境数据进行整合，去除环境效应，得到个体稳定遗传的表型。BLUP是表型处理的常用做法。R包lme4中lmer函数是BLUP分析常用的方法，在很多NG文章都引用了该方法。

## 无重复数据

对于多年或者多环境的无重复数据，或是使用了重复数据平均值来代表品种表型值的数据的情况。以下是多年表型计算BLUP值的示例

**数据准备**：总之根据自己的数据情况，整理成类似下面的格式

```R
library(lme4)
phes <- readxl::read_xlsx("~/Desktop/表型汇总.xlsx")
pheName <- "Heading_date"

# 将多年或者多环境数据从宽格式转换为长格式，以符合lmer的需要
phe <- phes %>% 
	select(Number, starts_with(pheName)) %>% 
    pivot_longer(-Number, names_to = "year", values_to = pheName)
head(phe)

##   Number Year            Heading_date
##   <chr>  <chr>                  <dbl>
## 1 C001   Heading_date_11           80
## 2 C001   Heading_date_12           76
## 3 C002   Heading_date_11           96
## 4 C002   Heading_date_12           96
## 5 C003   Heading_date_11           82
## 6 C003   Heading_date_12           78
```

计算多年表型BLUP值

```R
#########BLUP#########

# 构建模型，拟合数据
f <- as.formula(paste(pheName, "~ (1|Number) + (1|year)"))
model <- lmer(f, data = phe)

# 提取随机效应的 BLUP 值
random_effects_blups <- ranef(model)
blups <- ranef(model)
lines <- blups$Number + model@beta # model@beta 为整体均值

# 将得到的BLUP值储存到数据框
res <- tibble(ID = rownames(lines), Heading_date_blup = lines$`(Intercept)`)
head(res)

## # A tibble: 6 × 2
##   ID    Heading_date_blup
##   <chr>             <dbl>
## 1 C001               79.3
## 2 C002               95.9
## 3 C003               81.1
## 4 C004               83.0
## 5 C005              101. 
## 6 C006               89.9
```

