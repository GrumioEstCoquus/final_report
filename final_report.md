Final Report
================
Joseph Burks
2024-05-05

# Data Set

The data set is from 538 a news website owned by ABC News with a focus
on opinion poll analysis, economics and politics. The data set utilized
is a collection of various demographic and economic information per
state, such as percentage of no white citizens, percentage of non
citizens and median incomes, as wells as, hate crime statistics.

``` r
library(tidyverse)
library(factoextra)
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
<p>Generated by <a href='https://github.com/dcomtois/summarytools'>summarytools</a> 1.0.1 (<a href='https://www.r-project.org/'>R</a> version 4.3.1)<br/>2024-05-05</p>
</div>

## Distribution of Variables

``` r
#ggplot(hate_crimes) +
```

## Hate Crimes by State

``` r
ggplot(hate_crimes) +
  aes(x = reorder(state, avg_hatecrimes_per_100k_fbi), y = avg_hatecrimes_per_100k_fbi) +
  geom_bar(position="dodge",stat="identity",fill="blue",cex=0.75) + 
  coord_flip() +
  labs(title = "Average Hate Crime per 100K by State",
       x = "State", y = "Average Hate Crime per 100K")
```

![](final_report_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
District of Columbia has by far the largest amount of hate crimes per
100k, more than double the second largest of Massachusetts. The states
with smallest amount of hate crimes is Wyoming, Georgia, and
Pennsylvania.

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
we decided to drop the state.

``` r
no_na_crimes <- hate_crimes
no_na_crimes[20,6] <- 0.01
no_na_crimes[25,6] <- 0.01
no_na_crimes[42,6] <- 0.03

no_na_crimes[is.na(no_na_crimes)] <- 0
no_na_crimes[12,12] <- NA

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

# Principle Component Analysis

``` r
no_na_crimes <- drop_na(no_na_crimes)
scaled_crimes <- apply(no_na_crimes[,2:12], 2, scale)
crimes_cov <- cov(scaled_crimes)
crimes_eigen <- eigen(crimes_cov)
str(crimes_eigen)
```

    ## List of 2
    ##  $ values : num [1:11] 4.424 3.124 1.409 0.535 0.393 ...
    ##  $ vectors: num [1:11, 1:11] -0.1915 -0.2175 -0.3834 0.0992 -0.3996 ...
    ##  - attr(*, "class")= chr "eigen"

According to Kaiser’s rule the first three principle components are
enough since they are the only ones with eigenvalues greater than 1

``` r
phi <-crimes_eigen$vectors[,1:3]
phi <- -phi
colnames(phi) <- c("PC1","PC2","PC3")

pc1 <- scaled_crimes %*% phi[,1]
pc2 <- scaled_crimes %*% phi[,2]
pc3 <- scaled_crimes %*% phi[,3]

PC <- data.frame(State = no_na_crimes$state,pc1,pc2,pc3)
head(PC)
```

    ##        State        pc1       pc2        pc3
    ## 1    Alabama -1.1662995  2.758740  0.6278585
    ## 2     Alaska -0.2996980 -1.401397 -1.3415843
    ## 3    Arizona  1.6810905  1.497702 -0.5755910
    ## 4   Arkansas -1.4600824  2.096989 -0.1324033
    ## 5 California  3.4668429  1.604581 -1.4271620
    ## 6   Colorado  0.6986251 -1.505261 -0.2123170

``` r
results <- princomp(scaled_crimes)
fviz_pca_biplot(results)
```

![](final_report_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

 The biggest contributions for the first principle component are
positive contributions from the variables
share_population_with_highschool_degree and median_income, and negative
contributions from share_unemployed_seasonal and share_white_poverty.
This indicates that the first principle component is a kind of measure
of the states’ economy. The biggest contributions for the second
component share_non_citizen, share_population_in_metro_areas, and the
two hate crime variables. Together the first two principle components
explain the nearly 70 percent of the variability in the data.

``` r
ggplot(PC, aes(pc1, pc2)) + 
  modelr::geom_ref_line(h = 0) +
  modelr::geom_ref_line(v = 0) +
  geom_text(aes(label = State), size = 3) +
  xlab("First Principal Component") + 
  ylab("Second Principal Component") + 
  ggtitle("First Two Principal Components of USArrests Data")
```

![](final_report_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
Mississippi appears to have above average economy, at least according to
the metrics in the data set, and New Hampshire has the worst. California
has above average economy and above average hate_crimes, percentage of
non citizens and urbanization.

``` r
PVE <- crimes_eigen$values/sum(crimes_eigen$values)
round(PVE,3)
```

    ##  [1] 0.402 0.284 0.128 0.049 0.036 0.029 0.027 0.016 0.011 0.010 0.008

``` r
crimes_eigen$values
```

    ##  [1] 4.42444417 3.12397095 1.40879154 0.53524759 0.39330449 0.32259320
    ##  [7] 0.29694751 0.18136749 0.11675561 0.10999128 0.08658616

The first 3 Principle components explain roughly 80 percent of the
variability of the data and are the only pcs with eigenvalues greater
than 1, so 3 principle components are enough.
