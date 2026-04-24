library(ieugwasr)

cat("\n=== ieu-a-1098 Info ===\n")
info <- gwasinfo("ieu-a-1098")
print(info)

cat("\n=== Search: Graves ===\n")
options(width = 150)
print(as.data.frame(head(gwassearch("Graves"))))

cat("\n=== Search: hyperthyroidism ===\n")
print(as.data.frame(head(gwassearch("hyperthyroidism"))))

cat("\n=== Search: thyrotoxicosis ===\n")
print(as.data.frame(head(gwassearch("thyrotoxicosis"))))
