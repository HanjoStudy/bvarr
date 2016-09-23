bvarr
=====


[![Travis-CI Build Status](https://travis-ci.org/bdemeshev/bvarr.svg?branch=master)](https://travis-ci.org/bdemeshev/bvarr)

Пакет `bvarr` может пригодиться для оценки BVAR моделей с сопряжённым нормальным - обратным Уишарта априорным распределением.

Пакет можно установить командами:
```R
install.packages("devtools")
devtools::install_github("bdemeshev/bvarr")
```

Простой пример оценки BVAR с сопряжённым нормальным - обратным Уишарта априорным распределением
```R
library("bvarr")
data(Yraw)
priors <- Carriero_priors(Yraw, p = 4)
model <- bvar_conjugate0(priors = priors)
summary_conjugate(model) 
forecast_conjugate(model, h = 2, output = "wide")
forecast_conjugate(model, out_of_sample = FALSE, include = "mean", level = NULL, type = "credible")
```

[Теория моделей BVAR](https://github.com/bdemeshev/bvar_om/raw/master/text/bvar_mapping/bvar_mapping.pdf)

Презентация [BVAR: Great Grimpen Mire](http://bdemeshev.github.io/bvar_om/bvar_mire_pres.html/) в Нижнем Новгороде 2016-09-22.


Цели пакета:

1. Хорошая документация

2. Гибкость

3. Разумные значения параметров по умолчанию

4. Робастность к мерзким матрицам


Модели BVAR также можно оценивать с помощью пакетов:

- [BMR](http://bayes.squarespace.com/bmr/)

- [MSBVAR](http://cran.r-project.org/web/packages/MSBVAR/) 

- [bvarsv](https://cran.r-project.org/web/packages/bvarsv/index.html) 

## English translation


The package `bvarr` may be useful for estimation BVARs with conjugate Normal-Inverse Wishart prior.

You may install the package usinge the commands:
```R
install.packages("devtools")
devtools::install_github("bdemeshev/bvarr")
```

Basic example of BVAR estimation with forecasting
```R
library("bvarr")
data(Yraw)
priors <- Carriero_priors(Yraw, p = 4)
model <- bvar_conjugate0(priors = priors)
summary_conjugate(model) 
forecast_conjugate(model, h = 2, output = "wide")
forecast_conjugate(model, out_of_sample = FALSE, include = "mean", level = NULL, type = "credible")
```

[Theory behind package](https://github.com/bdemeshev/bvar_om/raw/master/text/bvar_mapping/bvar_mapping.pdf)

Presentation [BVAR: Great Grimpen Mire](http://bdemeshev.github.io/bvar_om/bvar_mire_pres.html/) in Nizhniy Novgorod 2016-09-22.

Goals of the package:

1. Good documentation

2. Versatile 

3. Reasonable default values

4. Robustness for bad matrices


You may also wish look at 

- [BMR](http://bayes.squarespace.com/bmr/)

- [MSBVAR](http://cran.r-project.org/web/packages/MSBVAR/) 

- [bvarsv](https://cran.r-project.org/web/packages/bvarsv/index.html) 




