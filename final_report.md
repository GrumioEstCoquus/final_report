Final Report
================
Joseph Burks, Sandali Wijeratne, Thin Han
2024-05-07

# Data Set

The data set is from 538 a news website owned by ABC News with a focus
on opinion poll analysis, economics and politics. The data set utilized
is a collection of various demographic and economic information per
state, such as percentage of no white citizens, percentage of non
citizens and median incomes, as wells as, hate crime statistics. The
dataset can be obtained here:
<https://github.com/fivethirtyeight/data/tree/master/hate-crimes>

``` r
library(tidyverse)
library(factoextra)
library(factoextra)
library(cluster)
hate_crimes <- readr::read_csv("C:/Users/josep/Downloads/hate_crimes.csv")
head(hate_crimes)
```

    ## # A tibble: 6 × 12
    ##   state     median_household_inc…¹ share_unemployed_sea…² share_population_in_…³
    ##   <chr>                      <dbl>                  <dbl>                  <dbl>
    ## 1 Alabama                    42278                  0.06                    0.64
    ## 2 Alaska                     67629                  0.064                   0.63
    ## 3 Arizona                    49254                  0.063                   0.9 
    ## 4 Arkansas                   44922                  0.052                   0.69
    ## 5 Californ…                  60487                  0.059                   0.97
    ## 6 Colorado                   60940                  0.04                    0.8 
    ## # ℹ abbreviated names: ¹​median_household_income, ²​share_unemployed_seasonal,
    ## #   ³​share_population_in_metro_areas
    ## # ℹ 8 more variables: share_population_with_high_school_degree <dbl>,
    ## #   share_non_citizen <dbl>, share_white_poverty <dbl>, gini_index <dbl>,
    ## #   share_non_white <dbl>, share_voters_voted_trump <dbl>,
    ## #   hate_crimes_per_100k_splc <dbl>, avg_hatecrimes_per_100k_fbi <dbl>

# Exlporatory Analysis

## Summary Statitics

``` r
summary_stats <- summarytools::descr(hate_crimes, round.digits = 2, transpose = TRUE)
summarytools::view(summary_stats, method = "render")
```

    ## Non-numerical variable(s) ignored: state

<div class="container st-container">
<h3>Descriptive Statistics</h3>
<h4>hate_crimes</h4>
<strong>N</strong>: 51
<br/>
<table class="table table-bordered table-striped st-table st-table-bordered st-table-striped st-descr-table ">
  <thead>
    <tr>
      <th></th>
      <th align="center" class="st-protect-top-border">Mean</th>
      <th align="center" class="st-protect-top-border">Std.Dev</th>
      <th align="center" class="st-protect-top-border">Min</th>
      <th align="center" class="st-protect-top-border">Q1</th>
      <th align="center" class="st-protect-top-border">Median</th>
      <th align="center" class="st-protect-top-border">Q3</th>
      <th align="center" class="st-protect-top-border">Max</th>
      <th align="center" class="st-protect-top-border">MAD</th>
      <th align="center" class="st-protect-top-border">IQR</th>
      <th align="center" class="st-protect-top-border">CV</th>
      <th align="center" class="st-protect-top-border">Skewness</th>
      <th align="center" class="st-protect-top-border">SE.Skewness</th>
      <th align="center" class="st-protect-top-border">Kurtosis</th>
      <th align="center" class="st-protect-top-border">N.Valid</th>
      <th align="center" class="st-protect-top-border">Pct.Valid</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>
        <strong>avg_hatecrimes_per_100k_fbi</strong></td>
      <td><span>2.37</span></td>
      <td><span>1.71</span></td>
      <td><span>0.27</span></td>
      <td><span>1.28</span></td>
      <td><span>1.99</span></td>
      <td><span>3.20</span></td>
      <td><span>10.95</span></td>
      <td><span>1.37</span></td>
      <td><span>1.89</span></td>
      <td><span>0.72</span></td>
      <td><span>2.49</span></td>
      <td><span>0.34</span></td>
      <td><span>10.08</span></td>
      <td><span>50</span></td>
      <td><span>98.04</span></td>
    </tr>
    <tr>
      <td>
        <strong>gini_index</strong></td>
      <td><span>0.45</span></td>
      <td><span>0.02</span></td>
      <td><span>0.42</span></td>
      <td><span>0.44</span></td>
      <td><span>0.45</span></td>
      <td><span>0.47</span></td>
      <td><span>0.53</span></td>
      <td><span>0.02</span></td>
      <td><span>0.03</span></td>
      <td><span>0.05</span></td>
      <td><span>0.95</span></td>
      <td><span>0.33</span></td>
      <td><span>2.13</span></td>
      <td><span>51</span></td>
      <td><span>100.00</span></td>
    </tr>
    <tr>
      <td>
        <strong>hate_crimes_per_100k_splc</strong></td>
      <td><span>0.30</span></td>
      <td><span>0.25</span></td>
      <td><span>0.07</span></td>
      <td><span>0.14</span></td>
      <td><span>0.23</span></td>
      <td><span>0.36</span></td>
      <td><span>1.52</span></td>
      <td><span>0.15</span></td>
      <td><span>0.21</span></td>
      <td><span>0.83</span></td>
      <td><span>2.64</span></td>
      <td><span>0.35</span></td>
      <td><span>9.29</span></td>
      <td><span>47</span></td>
      <td><span>92.16</span></td>
    </tr>
    <tr>
      <td>
        <strong>median_household_income</strong></td>
      <td><span>55223.61</span></td>
      <td><span>9208.48</span></td>
      <td><span>35521.00</span></td>
      <td><span>48060.00</span></td>
      <td><span>54916.00</span></td>
      <td><span>60730.00</span></td>
      <td><span>76165.00</span></td>
      <td><span>8619.84</span></td>
      <td><span>12062.00</span></td>
      <td><span>0.17</span></td>
      <td><span>0.19</span></td>
      <td><span>0.33</span></td>
      <td><span>-0.58</span></td>
      <td><span>51</span></td>
      <td><span>100.00</span></td>
    </tr>
    <tr>
      <td>
        <strong>share_non_citizen</strong></td>
      <td><span>0.05</span></td>
      <td><span>0.03</span></td>
      <td><span>0.01</span></td>
      <td><span>0.03</span></td>
      <td><span>0.04</span></td>
      <td><span>0.08</span></td>
      <td><span>0.13</span></td>
      <td><span>0.02</span></td>
      <td><span>0.05</span></td>
      <td><span>0.57</span></td>
      <td><span>0.56</span></td>
      <td><span>0.34</span></td>
      <td><span>-0.77</span></td>
      <td><span>48</span></td>
      <td><span>94.12</span></td>
    </tr>
    <tr>
      <td>
        <strong>share_non_white</strong></td>
      <td><span>0.32</span></td>
      <td><span>0.16</span></td>
      <td><span>0.06</span></td>
      <td><span>0.19</span></td>
      <td><span>0.28</span></td>
      <td><span>0.42</span></td>
      <td><span>0.81</span></td>
      <td><span>0.16</span></td>
      <td><span>0.22</span></td>
      <td><span>0.52</span></td>
      <td><span>0.66</span></td>
      <td><span>0.33</span></td>
      <td><span>0.05</span></td>
      <td><span>51</span></td>
      <td><span>100.00</span></td>
    </tr>
    <tr>
      <td>
        <strong>share_population_in_metro_areas</strong></td>
      <td><span>0.75</span></td>
      <td><span>0.18</span></td>
      <td><span>0.31</span></td>
      <td><span>0.63</span></td>
      <td><span>0.79</span></td>
      <td><span>0.90</span></td>
      <td><span>1.00</span></td>
      <td><span>0.19</span></td>
      <td><span>0.26</span></td>
      <td><span>0.24</span></td>
      <td><span>-0.62</span></td>
      <td><span>0.33</span></td>
      <td><span>-0.45</span></td>
      <td><span>51</span></td>
      <td><span>100.00</span></td>
    </tr>
    <tr>
      <td>
        <strong>share_population_with_high_school_degree</strong></td>
      <td><span>0.87</span></td>
      <td><span>0.03</span></td>
      <td><span>0.80</span></td>
      <td><span>0.84</span></td>
      <td><span>0.87</span></td>
      <td><span>0.90</span></td>
      <td><span>0.92</span></td>
      <td><span>0.04</span></td>
      <td><span>0.06</span></td>
      <td><span>0.04</span></td>
      <td><span>-0.37</span></td>
      <td><span>0.33</span></td>
      <td><span>-1.13</span></td>
      <td><span>51</span></td>
      <td><span>100.00</span></td>
    </tr>
    <tr>
      <td>
        <strong>share_unemployed_seasonal</strong></td>
      <td><span>0.05</span></td>
      <td><span>0.01</span></td>
      <td><span>0.03</span></td>
      <td><span>0.04</span></td>
      <td><span>0.05</span></td>
      <td><span>0.06</span></td>
      <td><span>0.07</span></td>
      <td><span>0.01</span></td>
      <td><span>0.02</span></td>
      <td><span>0.22</span></td>
      <td><span>0.03</span></td>
      <td><span>0.33</span></td>
      <td><span>-0.79</span></td>
      <td><span>51</span></td>
      <td><span>100.00</span></td>
    </tr>
    <tr>
      <td>
        <strong>share_voters_voted_trump</strong></td>
      <td><span>0.49</span></td>
      <td><span>0.12</span></td>
      <td><span>0.04</span></td>
      <td><span>0.41</span></td>
      <td><span>0.49</span></td>
      <td><span>0.58</span></td>
      <td><span>0.70</span></td>
      <td><span>0.12</span></td>
      <td><span>0.16</span></td>
      <td><span>0.24</span></td>
      <td><span>-0.92</span></td>
      <td><span>0.33</span></td>
      <td><span>2.13</span></td>
      <td><span>51</span></td>
      <td><span>100.00</span></td>
    </tr>
    <tr>
      <td>
        <strong>share_white_poverty</strong></td>
      <td><span>0.09</span></td>
      <td><span>0.02</span></td>
      <td><span>0.04</span></td>
      <td><span>0.07</span></td>
      <td><span>0.09</span></td>
      <td><span>0.10</span></td>
      <td><span>0.17</span></td>
      <td><span>0.01</span></td>
      <td><span>0.02</span></td>
      <td><span>0.27</span></td>
      <td><span>0.61</span></td>
      <td><span>0.33</span></td>
      <td><span>0.68</span></td>
      <td><span>51</span></td>
      <td><span>100.00</span></td>
    </tr>
  </tbody>
</table>
<p>Generated by <a href='https://github.com/dcomtois/summarytools'>summarytools</a> 1.0.1 (<a href='https://www.r-project.org/'>R</a> version 4.3.1)<br/>2024-05-07</p>
</div>

## Distribution of Variables

``` r
par(mfrow =c(3,3))
colnames <- dimnames(hate_crimes)[[2]]
for (i in 2:12) {
    hist(as.numeric(unlist(hate_crimes[,i])),  main=colnames[i], col="cyan", xlab = colnames[i])
}
```

![](final_report_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

## Variables by State

``` r
for (i in 2:12){
  print(ggplot(hate_crimes) +
    aes(x = reorder(state, as.numeric(unlist(hate_crimes[,i]))), y = as.numeric(unlist(hate_crimes[,i]))) +
    geom_bar(position="dodge",stat="identity",fill="blue",cex=0.75) + 
    coord_flip() +
    labs(title = paste(colnames[i], "by State"),
         x = "State", y = colnames[i]))
}
```

![](final_report_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-6.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-7.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-8.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-9.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-10.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-4-11.png)<!-- -->

## Missing Values

``` r
missing.values <- hate_crimes |>
  gather(key = "key", value = "val") |>
  mutate(is.missing = is.na(val)) |>
  group_by(key, is.missing) |>
  summarise(num.missing = n()) |>
  filter(is.missing == T) |>
  select(-is.missing) |>
  arrange(desc(num.missing))
```

    ## `summarise()` has grouped output by 'key'. You can override using the `.groups`
    ## argument.

``` r
missing.values |>
  ggplot() +
  geom_bar(aes(x=key, y = num.missing), stat = "identity",fill = "orange") +
  labs(x = "variable", y = "number of missing values", title="Number of missing values") +
  theme(axis.text.x = element_text(angle = 45, hjust =1))
```

![](final_report_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> Luckily
the data set does not contain many missing values, and in fact many
seems to be the result of human error. 538’s given source for share of
non citizens, Kaiser Family Foundation does contain the percentage for
the non citizens for the missing states, which can be added to the data
set. Similarly all the missing data for hate crimes per 100k from the
Southern Poverty Law Center is reported, just that in the 4 states with
missing data the Southern Poverty Law Center actually reported 0
instincances of hate crime. Although this is likely due to how there
data collection relied on people reporting directly to them. The last
missing data point is a result of Hawaii not sharing hate crime
information with the FBI, meaning that the data is not missing
completely at random. Since Hawaii is the only state with missing data,
we decided to drop the state. We while also drop DC, since it is not a
state and including it drastically changes the results.

``` r
no_na_crimes <- hate_crimes
no_na_crimes[20,6] <- 0.01
no_na_crimes[25,6] <- 0.01
no_na_crimes[42,6] <- 0.03

no_na_crimes[is.na(no_na_crimes)] <- 0
no_na_crimes[12,12] <- NA
no_na_crimes <- no_na_crimes[-c(9),]

missing.values2 <- no_na_crimes |>
  gather(key = "key", value = "val") |>
  mutate(is.missing = is.na(val)) |>
  group_by(key, is.missing) |>
  summarise(num.missing = n()) |>
  filter(is.missing == T) |>
  select(-is.missing) |>
  arrange(desc(num.missing))
```

    ## `summarise()` has grouped output by 'key'. You can override using the `.groups`
    ## argument.

``` r
missing.values2
```

    ## # A tibble: 1 × 2
    ## # Groups:   key [1]
    ##   key                         num.missing
    ##   <chr>                             <int>
    ## 1 avg_hatecrimes_per_100k_fbi           1

# Multiple Linear Regression

``` r
lm_crimes <- lm(avg_hatecrimes_per_100k_fbi ~ . -state , data = no_na_crimes)
summary(lm_crimes)
```

    ## 
    ## Call:
    ## lm(formula = avg_hatecrimes_per_100k_fbi ~ . - state, data = no_na_crimes)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9975 -0.5946  0.1179  0.6704  2.4025 
    ## 
    ## Coefficients:
    ##                                            Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                              -2.224e+01  1.545e+01  -1.440   0.1581
    ## median_household_income                   3.864e-05  4.113e-05   0.939   0.3535
    ## share_unemployed_seasonal                 9.279e+00  2.249e+01   0.413   0.6822
    ## share_population_in_metro_areas          -6.511e-01  1.558e+00  -0.418   0.6784
    ## share_population_with_high_school_degree  9.709e+00  1.145e+01   0.848   0.4019
    ## share_non_citizen                         2.284e+01  1.010e+01   2.261   0.0295
    ## share_white_poverty                       1.065e+01  1.362e+01   0.782   0.4392
    ## gini_index                                2.414e+01  1.445e+01   1.670   0.1031
    ## share_non_white                          -3.846e+00  2.330e+00  -1.651   0.1071
    ## share_voters_voted_trump                  3.173e+00  2.921e+00   1.086   0.2841
    ## hate_crimes_per_100k_splc                 1.518e+00  1.153e+00   1.317   0.1956
    ##                                           
    ## (Intercept)                               
    ## median_household_income                   
    ## share_unemployed_seasonal                 
    ## share_population_in_metro_areas           
    ## share_population_with_high_school_degree  
    ## share_non_citizen                        *
    ## share_white_poverty                       
    ## gini_index                                
    ## share_non_white                           
    ## share_voters_voted_trump                  
    ## hate_crimes_per_100k_splc                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.083 on 38 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.3524, Adjusted R-squared:  0.1819 
    ## F-statistic: 2.068 on 10 and 38 DF,  p-value: 0.05261

The only significant predictor is share_non_citizen with a p-value
0.0295 and slope 22.84. The R^2 is around 0.3523 and was obtained from
fitting hate crimes from fbi vs all other variables. This R^2 is very
small, and performing stepwise variable selection, will only result in a
regression model with a worse R^2.

``` r
plot(lm_crimes)
```

![](final_report_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->![](final_report_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

The qq plot looks mostly alright then is some deviation from the line in
the right tail. The Residual plot does have a u shape, indicating that
either the assumption of linearity or homoscedasticity was violated.

``` r
plot(avg_hatecrimes_per_100k_fbi ~ share_non_citizen, data = no_na_crimes)
```

![](final_report_files/figure-gfm/unnamed-chunk-9-1.png)<!-- --> There
is no obvious pattern in the avg_hatecrimes_per_100k_fbi vs
share_non_citizen so polynomial regression is would likely not improve
upon linear regression. None of the variables seem like they would be
good predictors of hate crimes, indicating that regression is likely not
a good technique, at least with respect to hate crimes.

# Principle Component Analysis

``` r
no_na_crimes <- drop_na(no_na_crimes)
scaled_crimes <- apply(no_na_crimes[,2:12], 2, scale)
crimes_cov <- cov(scaled_crimes)
crimes_eigen <- eigen(crimes_cov)
str(crimes_eigen)
```

    ## List of 2
    ##  $ values : num [1:11] 3.757 3.342 1.388 0.799 0.481 ...
    ##  $ vectors: num [1:11, 1:11] -0.0471 0.3417 0.4118 -0.3466 0.4089 ...
    ##  - attr(*, "class")= chr "eigen"

According to Kaiser’s rule the first three principle components are
enough since they are the only ones with eigenvalues greater than 1

``` r
phi <-crimes_eigen$vectors[,1:3]

colnames(phi) <- c("PC1","PC2","PC3")

pc1 <- scaled_crimes %*% phi[,1]
pc2 <- scaled_crimes %*% phi[,2]
pc3 <- scaled_crimes %*% phi[,3]

PC <- data.frame(State = no_na_crimes$state,pc1,pc2,pc3)# REMOVED 3RD PC
head(PC)
```

    ##        State          pc1        pc2        pc3
    ## 1    Alabama  0.447603881 -3.1780209  0.1609234
    ## 2     Alaska -0.797816329  1.2190309 -1.7994691
    ## 3    Arizona  2.376928173 -0.1527681  0.1960283
    ## 4   Arkansas  0.012088819 -2.6144792 -0.6038205
    ## 5 California  4.296824306  0.8709089 -0.3089248
    ## 6   Colorado -0.004830473  1.8753252  0.2237174

``` r
results <- princomp(scaled_crimes, fix_sign = FALSE)
fviz_pca_biplot(results)
```

![](final_report_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

 The biggest contributions for the first principle component are
positive contributions from the variables median_income, both hate
crimes statistics and negative contributions from share_white_poverty.
This indicates that the first principle component is a kind of measure
hate crime rate and the average family economic status. The biggest
contributions for the second component share_nonwhite,
share_non_citizen, share_population_in_metro_areas, Together the first
two principle components explain the nearly 65 percent of the
variability in the data.

``` r
ggplot(PC, aes(pc1, pc2)) + 
  modelr::geom_ref_line(h = 0) +
  modelr::geom_ref_line(v = 0) +
  geom_text(aes(label = State), size = 3) +
  xlab("First Principal Component") + 
  ylab("Second Principal Component") + 
  ggtitle("First Two Principal Components of USArrests Data")
```

![](final_report_files/figure-gfm/unnamed-chunk-13-1.png)<!-- --> West
Virginia appears to have below household economic status and rate of
hate crimes. It also has a below average share of non white citizens and
non citizens. Maryland has above average economy and above average
hate_crimes, percentage of non citizens and urbanization.

``` r
PVE <- crimes_eigen$values/sum(crimes_eigen$values)
round(PVE,3)
```

    ##  [1] 0.342 0.304 0.126 0.073 0.044 0.032 0.027 0.021 0.012 0.011 0.009

``` r
crimes_eigen$values
```

    ##  [1] 3.75668189 3.34231967 1.38757129 0.79936774 0.48120612 0.34985911
    ##  [7] 0.30176244 0.23161315 0.13017449 0.12371693 0.09572716

The first 3 Principle components explain roughly 77 percent of the
variability of the data and are the only pcs with eigenvalues greater
than 1, so 3 principle components are enough to represent the data set
while significantly reducing the dimensions.

# Cluster Analysis

# K-Means Clustering

``` r
hate_crimes <- read.csv("C:/Users/josep/Downloads/hate_crimes.csv",header=TRUE,row.names="state")
hate_crimes[20,5] <- 0.01
hate_crimes[25,5] <- 0.01
hate_crimes[42,5] <- 0.03

#no_na_crimes <- drop_na(no_na_crimes)
#no_na_crimes
#scaled_crimes <- scale(no_na_crimes)
head(scaled_crimes)
```

    ##      median_household_income share_unemployed_seasonal
    ## [1,]              -1.3890165                 1.0081757
    ## [2,]               1.4616133                 1.3933656
    ## [3,]              -0.6045901                 1.2970681
    ## [4,]              -1.0917081                 0.2377958
    ## [5,]               0.6585208                 0.9118782
    ## [6,]               0.7094590                -0.9177739
    ##      share_population_in_metro_areas share_population_with_high_school_degree
    ## [1,]                      -0.5772940                               -1.3769565
    ## [2,]                      -0.6323279                                1.3265294
    ## [3,]                       0.8535865                               -0.7664919
    ## [4,]                      -0.3021247                               -1.2897472
    ## [5,]                       1.2388236                               -1.8130026
    ## [6,]                       0.3032478                                0.7160648
    ##      share_non_citizen share_white_poverty gini_index share_non_white
    ## [1,]        -0.9944839           1.1214842  1.0905787      0.35147457
    ## [2,]        -0.3447544          -1.3954345 -1.7190090      0.83563431
    ## [3,]         1.6044341          -0.1369752  0.1353189      1.31979406
    ## [4,]        -0.3447544           1.1214842  0.3038942     -0.27101653
    ## [5,]         2.5790284          -0.1369752  1.0343870      2.14978219
    ## [6,]         0.3049751          -0.9759481  0.2477024      0.07481186
    ##      share_voters_voted_trump hate_crimes_per_100k_splc
    ## [1,]               1.30002752               -0.72949195
    ## [2,]               0.27589009               -0.63259858
    ## [3,]              -0.03135115               -0.19103327
    ## [4,]               0.99278629               -1.03681376
    ## [5,]              -1.77238479               -0.02602561
    ## [6,]              -0.64583361                0.70315910
    ##      avg_hatecrimes_per_100k_fbi
    ## [1,]                  -0.3224643
    ## [2,]                  -0.4475386
    ## [3,]                   1.0205224
    ## [4,]                  -1.1054414
    ## [5,]                   0.1717623
    ## [6,]                   0.5115387

``` r
scaled_crimes <- scaled_crimes[-9,]
```

We will use Total Within Sum of Square vs. Number of clusters graph to
find out the number of clusters that are suitable for our dataset.

``` r
fviz_nbclust(scaled_crimes, kmeans, nstart = 25, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)
```

![](final_report_files/figure-gfm/unnamed-chunk-16-1.png)<!-- --> From
the graph, we conclude that the best number of clusters is four.
Therefore, we apply k-means algorithm on the data set with four clusters
as a parameter. The following print-outs are cluster assignments found.

``` r
set.seed(314)
km.res <- kmeans(scaled_crimes,4,nstart=25)
print(km.res$cluster)
```

    ##  [1] 3 1 2 3 2 4 4 4 2 1 2 1 1 1 3 3 1 4 4 4 4 3 1 1 1 2 1 2 2 2 3 1 1 3 4 1 2 3
    ## [39] 1 3 2 1 1 4 4 3 1 1

We can visualize the clusters as following in two dimensional space.

``` r
fviz_cluster(km.res, 
             data = scaled_crimes,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal()
)
```

![](final_report_files/figure-gfm/unnamed-chunk-18-1.png)<!-- --> \#
K-Medoids Clustering

We will apply one more method, K-Medoids clustering to find the clusters
in our data set.

``` r
fviz_nbclust(scaled_crimes, pam, method = "silhouette") +
  theme_classic()
```

![](final_report_files/figure-gfm/unnamed-chunk-19-1.png)<!-- --> Using
Silhouette Width method, we get the same result that four clusters is
the best for the data set.

``` r
pam.res <- pam(scaled_crimes, 4)
kmedoid.df <- cbind(scaled_crimes, cluster = pam.res$cluster)
fviz_cluster(pam.res, 
             data = kmedoid.df,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal()
)
```

![](final_report_files/figure-gfm/unnamed-chunk-20-1.png)<!-- --> From
KMedoids clustering, South Carolina, Wisconsin, Virginia and Washington
comes out as medoids.

``` r
med_sil <- eclust(scaled_crimes, "pam", k = 4, hc_metric = "euclidean",
                 hc_method = "ward.D2", graph = FALSE)
fviz_silhouette(med_sil, palette = "jco",
                ggtheme = theme_classic())
```

    ##   cluster size ave.sil.width
    ## 1       1   10          0.32
    ## 2       2   18          0.21
    ## 3       3    9          0.23
    ## 4       4   11          0.16

![](final_report_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

The average silhouette width of KMedoids is 0.23, which is very small,
and indicates that clustering may not be the best option for explaining
the variance in the data.

## Hierarchical Clustering on First Two Principle Components

``` r
row.names(PC) <- (PC$State)
PC <- subset(PC, select = -c(State))
pc.dist <- dist(PC[,1:2], method = "euclidean")

as.matrix(pc.dist)[1:6,1:6]
```

    ##              Alabama   Alaska  Arizona  Arkansas California Colorado
    ## Alabama    0.0000000 4.570026 3.588098 0.7122167   5.586621 5.073559
    ## Alaska     4.5700258 0.000000 3.458444 3.9181304   5.106521 1.029344
    ## Arizona    3.5880980 3.458444 0.000000 3.4135739   2.175756 3.128248
    ## Arkansas   0.7122167 3.918130 3.413574 0.0000000   5.523304 4.489836
    ## California 5.5866207 5.106521 2.175756 5.5233041   0.000000 4.417362
    ## Colorado   5.0735593 1.029344 3.128248 4.4898363   4.417362 0.000000

``` r
pc.hc <- hclust(d = pc.dist, method = "ward.D2")

fviz_dend(pc.hc, cex = 0.5)
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.
    ## ℹ The deprecated feature was likely used in the factoextra package.
    ##   Please report the issue at <https://github.com/kassambara/factoextra/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](final_report_files/figure-gfm/unnamed-chunk-22-1.png)<!-- --> Using
Ward’s method, there seems to be four main clusters of states.

``` r
grp <- cutree(pc.hc, k = 4)
fviz_cluster(list(data = PC[,1:2], cluster = grp),
             pallete = c("blue","orange","red","pink"),
             ellipse.type = "convex",
             repel = TRUE,
             show.clust.cent = FALSE, ggtheme = theme_minimal())
```

![](final_report_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
res.coph <- cophenetic(pc.hc)
cor(pc.dist,res.coph)
```

    ## [1] 0.653878

The correlation between the cophenetic distance and the original
distance is around 0.62 which is not the largest amount of correlation,
so the clustering solution may not accurately reflect the data.

``` r
res.hc2 <- hclust(pc.dist,method = "average")
cor(pc.dist, cophenetic(res.hc2))
```

    ## [1] 0.684254

``` r
fviz_dend(res.hc2, cex = 0.5)
```

![](final_report_files/figure-gfm/unnamed-chunk-25-1.png)<!-- --> The
correlation between cophenetic distance from using average linkage
method and the original distance is around 0.64 is larger than the one
from Ward’s method, although it is still not that large. Also the
dendogram obtained with average linkage does not have as clear a place
to “cut” the tree.

## Kmeans with First 2 Principle Components

``` r
fviz_nbclust(PC[,1:2], kmeans, nstart = 25, method = "wss")
```

![](final_report_files/figure-gfm/unnamed-chunk-26-1.png)<!-- --> The
graph of the function of WSS vs K, indicates that 4 is the optimal
number of clusters for kmeans

``` r
set.seed(42)
km.res <- kmeans(PC, 4, nstart = 25)
print(km.res)
```

    ## K-means clustering with 4 clusters of sizes 10, 11, 10, 18
    ## 
    ## Cluster means:
    ##          pc1        pc2        pc3
    ## 1  2.5411565 -0.1175661 -0.5308790
    ## 2  0.7258250  2.0951723  0.6852901
    ## 3  0.2082508 -2.7048290  0.3607785
    ## 4 -1.9710082  0.2876142 -0.3242881
    ## 
    ## Clustering vector:
    ##        Alabama         Alaska        Arizona       Arkansas     California 
    ##              3              4              1              3              1 
    ##       Colorado    Connecticut       Delaware        Florida        Georgia 
    ##              2              2              2              1              1 
    ##          Idaho       Illinois        Indiana           Iowa         Kansas 
    ##              4              1              4              4              4 
    ##       Kentucky      Louisiana          Maine       Maryland  Massachusetts 
    ##              3              3              4              2              2 
    ##       Michigan      Minnesota    Mississippi       Missouri        Montana 
    ##              2              2              3              4              4 
    ##       Nebraska         Nevada  New Hampshire     New Jersey     New Mexico 
    ##              4              1              4              2              1 
    ##       New York North Carolina   North Dakota           Ohio       Oklahoma 
    ##              1              3              4              4              3 
    ##         Oregon   Pennsylvania   Rhode Island South Carolina   South Dakota 
    ##              2              4              1              3              4 
    ##      Tennessee          Texas           Utah        Vermont       Virginia 
    ##              3              1              4              4              2 
    ##     Washington  West Virginia      Wisconsin        Wyoming 
    ##              2              3              4              4 
    ## 
    ## Within cluster sum of squares by cluster:
    ## [1] 16.77693 34.56740 27.81483 53.21098
    ##  (between_SS / total_SS =  67.5 %)
    ## 
    ## Available components:
    ## 
    ## [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    ## [6] "betweenss"    "size"         "iter"         "ifault"

The clusters obtained explain roughly 77.4 percent of the variance in
the data. This is a decent size, however, it is important to remember
that clustering was done with respect to the first 2 principle
components which explain less than 70 percent of the variance

``` r
fviz_cluster(km.res, data = PC[,1:2], xlab = "PC1 (34.2%)", ylab = "PC2 (30.4%)")
```

![](final_report_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

The clusters on the plane spanned by the first 2 principle component, do
good job separating observations. However, the the first two principle
components do not explain enough of the variance in the data to
clustering with only them in mind.

``` r
mean_sil <- eclust(PC[,1:2], "kmeans", k = 4, hc_metric = "euclidean",
                 hc_method = "ward.D2", graph = FALSE)
fviz_silhouette(mean_sil, palette = "jco",
                ggtheme = theme_classic())
```

    ##   cluster size ave.sil.width
    ## 1       1   11          0.40
    ## 2       2   10          0.46
    ## 3       3   18          0.41
    ## 4       4   10          0.52

![](final_report_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

The average silhouette width is 0.44 and no clusters have not assigned
points “incorrectly” according to this metric. The width is not the
largest suggesting that the clusters might not be the most defined.

## Kmedoids on Principal components

``` r
# Performing k-medoids clustering on the principal component scores
set.seed(123) 
k <- 4
kmedoids_clusters <- pam(PC[,1:2], k = k)

#Clusters
fviz_cluster(kmedoids_clusters, data = PC_scores, geom = "point",
             stand = FALSE, ellipse.type = "convex", ellipse = TRUE,
             repel = TRUE)
```

![](final_report_files/figure-gfm/unnamed-chunk-30-1.png)<!-- --> The
clusters provide a good amount of separation, and it is similar to the
clusters obtained through kmeans, although a few points are assigned to
different clusters.

``` r
PC_with_clusters <- cbind(PC[,1:2], Cluster = as.factor(kmedoids_clusters$clustering))
head(PC_with_clusters)
```

    ##                     pc1        pc2 Cluster
    ## Alabama     0.447603881 -3.1780209       1
    ## Alaska     -0.797816329  1.2190309       2
    ## Arizona     2.376928173 -0.1527681       3
    ## Arkansas    0.012088819 -2.6144792       1
    ## California  4.296824306  0.8709089       3
    ## Colorado   -0.004830473  1.8753252       2

``` r
med_sil <- eclust(PC[,1:2], "pam", k = 4, hc_metric = "euclidean",
                 hc_method = "ward.D2", graph = FALSE)
fviz_silhouette(med_sil, palette = "jco",
                ggtheme = theme_classic())
```

    ##   cluster size ave.sil.width
    ## 1       1   10          0.46
    ## 2       2   13          0.33
    ## 3       3   10          0.52
    ## 4       4   16          0.46

![](final_report_files/figure-gfm/unnamed-chunk-32-1.png)<!-- --> The
kmedoid method does not “misclassify” any points and has an average
silhouette width of 0.44, which is the same as the width obtained
through kmeans.

## Cluster Assesment

``` r
library("clustertend")
```

    ## Package `clustertend` is deprecated.  Use package `hopkins` instead.

``` r
set.seed(123)
hopkins(PC, n = nrow(PC)-1)
```

    ## Warning: Package `clustertend` is deprecated.  Use package `hopkins` instead.

    ## $H
    ## [1] 0.4299043

The Hopkins statistic is pretty close to 0.5 for the first two principle
components, indicating likely no clusters exist in the data.

``` r
fviz_dist(dist(PC), show_labels = FALSE) +
  labs(title = "PCs")
```

![](final_report_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

VAT algorithm seems that there could be two clusters in the data.
Although that’s just my opinion

``` r
library("NbClust")
nb <- NbClust(PC, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")
```

![](final_report_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

    ## *** : The Hubert index is a graphical method of determining the number of clusters.
    ##                 In the plot of Hubert index, we seek a significant knee that corresponds to a 
    ##                 significant increase of the value of the measure i.e the significant peak in Hubert
    ##                 index second differences plot. 
    ## 

![](final_report_files/figure-gfm/unnamed-chunk-35-2.png)<!-- -->

    ## *** : The D index is a graphical method of determining the number of clusters. 
    ##                 In the plot of D index, we seek a significant knee (the significant peak in Dindex
    ##                 second differences plot) that corresponds to a significant increase of the value of
    ##                 the measure. 
    ##  
    ## ******************************************************************* 
    ## * Among all indices:                                                
    ## * 4 proposed 2 as the best number of clusters 
    ## * 9 proposed 3 as the best number of clusters 
    ## * 5 proposed 5 as the best number of clusters 
    ## * 1 proposed 7 as the best number of clusters 
    ## * 1 proposed 8 as the best number of clusters 
    ## * 3 proposed 10 as the best number of clusters 
    ## 
    ##                    ***** Conclusion *****                            
    ##  
    ## * According to the majority rule, the best number of clusters is  3 
    ##  
    ##  
    ## *******************************************************************

3 and 7 seems to be the optimal number of clusters.

``` r
library(clValid)
```

    ## Warning: package 'clValid' was built under R version 4.3.3

``` r
clmethods <- c("hierarchical", "kmeans", "pam")
intern <- clValid(PC, nClust = 2:6,
                  clMethods = clmethods, validation = "internal")

summary(intern)
```

    ## 
    ## Clustering Methods:
    ##  hierarchical kmeans pam 
    ## 
    ## Cluster sizes:
    ##  2 3 4 5 6 
    ## 
    ## Validation Measures:
    ##                                  2       3       4       5       6
    ##                                                                   
    ## hierarchical Connectivity  10.8762 16.1933 21.2274 24.0313 27.9492
    ##              Dunn           0.1609  0.1970  0.2596  0.2813  0.2813
    ##              Silhouette     0.2950  0.3356  0.3515  0.3175  0.3101
    ## kmeans       Connectivity  12.6040 21.8472 21.2274 24.0313 27.9492
    ##              Dunn           0.0542  0.1191  0.2596  0.2813  0.2813
    ##              Silhouette     0.2936  0.3374  0.3515  0.3175  0.3101
    ## pam          Connectivity  22.1060 18.9738 22.5385 34.5937 40.7294
    ##              Dunn           0.0886  0.1904  0.0874  0.0874  0.0874
    ##              Silhouette     0.2752  0.3392  0.3498  0.3001  0.2795
    ## 
    ## Optimal Scores:
    ## 
    ##              Score   Method       Clusters
    ## Connectivity 10.8762 hierarchical 2       
    ## Dunn          0.2813 hierarchical 5       
    ## Silhouette    0.3515 hierarchical 4

The 3 metrics do not agree on the best approach. According to
connectivity hierarchical with 2 clusters is the best, according to Dunn
kmeans with 6 clusters is the best and according to Silhouette kmeans
with 4 clusters is the best.

## Conclusion

The three statistical methods we used in this report had varying levels
of success in accurately representing patterns in the data. Multiple
linear regression resulted in a very small R^2, indicating that there is
very little linear relationship between hate crimes per 100k and the
rest of the variables. Principle component analysis was much more
successful. The first three principle components were able to explain
nearly 80 percent of the variance in the data, effectively reducing the
dimensions of the data. Clustering had little success. Clustering on the
original data yielded groupings with small average silhouette widths.
Clustering with respect to the first two principle components had larger
average silhouette width, but the results must be taken with a grain of
salt since the first two principle components only explain roughly 65
percent of the variability in the data. Overall principle component
analysis by itself seems to be the most effect method to represent the
data set.
