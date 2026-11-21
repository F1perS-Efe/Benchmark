# SigBenchR

Bir gen imzasının (gen listesi) immünoterapi yanıtını tahmin etmede birden fazla veri setinde gerçekten güvenilir olup olmadığını istatistiksel olarak sınayan, tekrarlanabilir doğrulama aracı.

> **Önemli:** Bu araç klinik karar desteği değildir. Sonuçlar araştırma amaçlıdır; küçük örneklemlerde belirsizdir. “Kanıt yok” sonucu da netçe raporlanır. Eşikler ve varsayılanlar sonuçlara bakıldıktan sonra değiştirilmemelidir.

## Kurulum (geliştirme)

```r
# Bağımlılıklar kurulduktan sonra:
devtools::install(".")
# veya
# remotes::install_local(".")
```

## Hızlı başlangıç

```r
library(SigBenchR)

data(datasets_demo)
data(signatures_demo)

res <- run_benchmark(
  signatures = signatures_demo,
  datasets = datasets_demo,
  methods = c("mean", "zscore", "median", "rank")
)

print(res)
plot_auc_heatmap(res, method = "mean")
```

Shiny arayüzü:

```r
SigBenchR::run_app()
```

## Kapsam (Faz 1 + 2A + 2B)

- Skorlama: `mean`, `zscore`, `median`, `rank`
- ROC-AUC + DeLong %95 GA / varyans; küçük kohort uyarısı
- PR-AUC + prevalans + tabakalı bootstrap GA (`PRROC`)
- Null permütasyon (`auc_mw` vektörize; p = (k+1)/(B+1))
- Eşleştirilmiş DeLong iki imza karşılaştırması (`compare_signatures`, BH)
- Gen kapsama denetimi + gri düşük-kapsama ısı haritası
- Skorlama dayanıklılığı: `check_sensitivity()` (eşik 0.10 keyfidir)
- Sentetik demo veri + Shiny

**v1 dışı:** Aynı imzada gen bazlı `+` / `-` işaretleri; takma ad eşleştirme.

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

## Atıflar

- Kus M, et al. TCGEx: a powerful visual interface for exploring and analyzing cancer gene expression data. *EMBO Reports* (2025). doi:[10.1038/s44319-025-00407-7](https://doi.org/10.1038/s44319-025-00407-7). Kod: [github.com/atakanekiz/TCGEx](https://github.com/atakanekiz/TCGEx).
- Ayers M, et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. *J Clin Invest* (2017).
- Robin X, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics* (2011).

Cytokine-5 imzası TCGEx çalışmasında Pan-ICI verisinden türetilmiştir; yazarlar imzanın üstünlüğünü iddia etmediklerini ve ek doğrulama gerektiğini belirtmiştir. Bu proje bağımsız bir doğrulama aracıdır.

## Lisans

GNU GPL (>= 3). Ayrıntılar için `LICENSE.md`. TCGEx ile uyumludur.
