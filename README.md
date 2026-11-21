# SigBenchR

Bir gen imzasının (gen listesi) immünoterapi yanıtını tahmin etmede birden fazla veri setinde gerçekten güvenilir olup olmadığını istatistiksel olarak sınayan, tekrarlanabilir doğrulama aracı.

> **Önemli:** Bu araç klinik karar desteği değildir. Sonuçlar araştırma amaçlıdır; küçük örneklemlerde belirsizdir. “Kanıt yok” sonucu da netçe raporlanır. Eşikler ve varsayılanlar sonuçlara bakıldıktan sonra değiştirilmemelidir.

## Kurulum

```r
# Geliştirme (paket kökünden):
devtools::install(".")
# veya
# remotes::install_local(".")
```

HTML rapor indirme için Suggests bağımlılığı:

```r
install.packages("rmarkdown")
```

## Hızlı başlangıç

### CLI

```r
library(SigBenchR)

data(datasets_demo)
data(signatures_demo)

res <- run_benchmark(
  signatures = signatures_demo,
  datasets = datasets_demo,
  methods = c("mean", "zscore", "median", "rank"),
  seed = 42L
)

print(res)
plot_auc_heatmap(res, method = "mean")

sm <- summarize_auc(res, method = "mean")
plot_auc_forest(res, signature = "Cytokine-5", method = "mean", summary = sm)
jaccard_signatures(signatures_demo)
auc_correlation(res, method = "mean")

# Tekrarlanabilir R kodu / HTML rapor
cat(generate_r_code(signature_name = "Cytokine-5", methods = c("mean", "zscore"), seed = 42L))
# export_report(res, sm, check_sensitivity(res), params_list = list(...), output_file = "report.html")
```

### Shiny arayüzü

```r
SigBenchR::run_app()
```

Arayüzde sonuç tablosu, ısı haritası, ROC, DeLong, ağırlıklı özet, forest, Jaccard, AUC korelasyonu, **HTML rapor indirme** ve **“Bu analizin R kodunu kopyala”** bulunur.

## Özellik özeti

- Skorlama: `mean`, `zscore`, `median`, `rank`
- ROC-AUC + DeLong %95 GA / varyans; küçük kohort uyarısı
- PR-AUC + prevalans + tabakalı bootstrap GA (`PRROC`)
- Null permütasyon (`auc_mw` vektörize; p = (k+1)/(B+1)); `null_p_adj_bh`
- Eşleştirilmiş DeLong iki imza karşılaştırması (`compare_signatures`, BH)
- Ağırlıklı özet: `summarize_auc()` — n-ağırlıklı + ters-varyans, Cochran `Q`, `I²`
- İmza örtüşmesi: `jaccard_signatures()`; veri setleri boyunca `auc_correlation()`
- Forest grafik: `plot_auc_forest()`
- Gen kapsama denetimi + gri düşük-kapsama ısı haritası
- Skorlama dayanıklılığı: `check_sensitivity()` (eşik 0.10 keyfidir)
- HTML rapor: `export_report()` (`inst/rmd/benchmark_report.Rmd`)
- Tekrarlanabilir betik: `generate_r_code()`
- Sentetik demo veri + Shiny

## SigBenchR v1 kapsam dışı

- Aynı imzada gen bazlı `+` / `-` işaretli imzalar **yok**
- Takma ad (alias) eşleştirme **yok**
- DerSimonian–Laird rastgele etkiler özeti **yok**

## Veri yapısı

Bioconductor nesnesi kullanılmaz. Her veri seti sade bir listedir:

```r
list(
  expr = <gen x örnek matris>,
  sample_meta = <data.frame: sample_id, response_label, ...>,
  info = list(normalization = "TPM", scale = "log2")  # veya "linear"
)
```

## Yanıt kodlama

- Varsayılan: `CR` / `PR` = 1, `PD` = 0
- Stabil hastalık (`SD`): varsayılan **dışlanır** (`sd_handling = "exclude"`)

## Dağıtım (shinyapps.io)

Henüz otomatik deploy yoktur. `rsconnect` ile yayınlamak için:

1. Paketi kurun (`devtools::install(".")` veya uzak depodan).
2. Repo kökündeki ince sarmalayıcı `app.R` yalnızca `SigBenchR::run_app()` çağırır (`R CMD check` için `.Rbuildignore` içinde `^app\.R$`).
3. Suggests paketlerini (özellikle `shiny`, `DT`, `rmarkdown`, `knitr`) hedef ortamda kurun; aksi halde UI özellikleri sessizce düşer veya rapor indirme hata verir.
4. Örnek:

```r
# install.packages("rsconnect")
rsconnect::setAccountInfo(name = "...", token = "...", secret = "...")
rsconnect::deployApp(
  appDir = ".",
  appFiles = c("app.R"),
  appName = "SigBenchR"
)
```

Not: `rsconnect` paket listesini doğru çıkarsın diye önce `SigBenchR`’ın Imports/Suggests bağımlılıklarının kurulabilir olduğundan emin olun. Canlı uygulamada kimliklendirilebilir hasta verisi yüklemeyin.

## Bilimsel yorum sınırı

Sonuçlar araştırma bağlamındadır; klinik karar aracı değildir. Küçük kohortlarda güven aralıkları geniştir. Eşikler (kapsama, dayanıklılık) sonuçlara bakıldıktan sonra değiştirilmemelidir.

## Atıflar

- Kus M, et al. TCGEx: a powerful visual interface for exploring and analyzing cancer gene expression data. *EMBO Reports* (2025). doi:[10.1038/s44319-025-00407-7](https://doi.org/10.1038/s44319-025-00407-7). Kod: [github.com/atakanekiz/TCGEx](https://github.com/atakanekiz/TCGEx).
- Ayers M, et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. *J Clin Invest* (2017).
- Robin X, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics* (2011).
- Keilwagen J, et al. PRROC: computing and visualizing precision-recall and receiver operating characteristic curves in R. *Bioinformatics* (2014).

Cytokine-5 imzası TCGEx çalışmasında Pan-ICI verisinden türetilmiştir; yazarlar imzanın üstünlüğünü iddia etmediklerini ve ek doğrulama gerektiğini belirtmiştir. Bu proje bağımsız bir doğrulama aracıdır.

## Lisans

GNU GPL (>= 3). Ayrıntılar için `LICENSE.md`. TCGEx ile uyumludur.
