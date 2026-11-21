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
  methods = c("mean", "zscore")
)

print(res)
plot_auc_heatmap(res, method = "mean")
```

Shiny arayüzü:

```r
SigBenchR::run_app()
```

## Faz 1 kapsamı

- Skorlama: `mean`, `zscore`
- ROC-AUC (sabit yön; `orient_score` sonrası her zaman `pROC` `direction = "<"`)
- Sentetik 2 demo veri seti + Cytokine-5 / IFNG6 / Random3
- Gen eşleştirme: büyük harfe çevrilmiş **tam eşleşme** (alias yok)
- Sade Shiny (`inst/shiny/app.R`)

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

`[LİSANS: benim seçimim; öneri GPL-3, TCGEx ile uyumlu]` — bkz. `LICENSE` dosyası. Kesin lisans Faz 3’te netleştirilecektir.
