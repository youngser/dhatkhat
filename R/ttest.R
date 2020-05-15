# paired sign test
#https://www.datanovia.com/en/lessons/sign-test-in-r/
library(rstatix)

labs <- c("tissue", "hemi", "Y", "region")

for (ll in labs) {
    df.long <- df.all %>% as_tibble() %>% select(i, index, contains(ll)) %>%
        gather(key = "group", value = "ari", -i, -index)

    df.long %>%
        group_by(group) %>%
        get_summary_stats(ari, type = "median_iqr")

    # stat.test <- df.long  %>%
    #     sign_test(ari ~ group) %>%
    #     add_significance()

    stat.test <- df.long  %>%
        t_test(ari ~ group, paired = TRUE, alternative = "greater") %>%
        filter(str_detect(group1, "SMS") | str_detect(group2, "SMS")) %>%
        print()

}



# paired t-test on sign test
#http://www.sthda.com/english/wiki/paired-samples-t-test-in-r
res <- t.test(before, after, paired = TRUE)
res

res <- t.test(weight ~ group, data = my_data, paired = TRUE)
res
res$p.value


res <- t.test(weight ~ group, data = my_data, paired = TRUE,
       alternative = "greater")
res
