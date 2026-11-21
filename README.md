# SigBenchR

[![R-CMD-check](https://github.com/F1perS-Efe/Benchmark/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/F1perS-Efe/Benchmark/actions/workflows/R-CMD-check.yaml)

SigBenchR, bir gen imzasının immünoterapi yanıt sınıflarını birden fazla kohortta ayırt edip etmediğini (discrimination) sınayan bir R paketidir. AUC burada klinik tahmin veya kalibrasyon modeli değildir.

Araç klinik karar desteği için değildir. Sonuçlar araştırma amaçlıdır; küçük örneklemlerde belirsizlik yüksektir. “Kanıt yok” da geçerli bir sonuçtur. Eşikleri sonuçlara baktıktan sonra değiştirmeyin.

Paket Bioconductor `Imports` kullanmaz. Veri setleri `list(expr, sample_meta, info)` sözleşmesini izler; isteğe bağlı GMT ve `SummarizedExperiment` adapter’ları vardır.

## Kurulum

R 4.2 veya üzeri gerekir. Depo kökünde:

```sh
Rscript setup.R
```

Betik `renv` proje kütüphanesini açar, `renv.lock` içindeki paketleri kurar, Shiny arayüzü, HTML rapor ve testler için Suggests paketlerini (`shiny`, `DT`, `rmarkdown`, `knitr`, `survival`, `testthat`) yükler ve `SigBenchR` paketini aynı kütüphaneye kurar. Ardından aynı klasörde açılan bir R oturumunda:

```r
library(SigBenchR)
SigBenchR::run_app()
```

Aynı adımlar elle:

```r
renv::restore()
renv::install(
  c("shiny", "DT", "rmarkdown", "knitr", "survival", "testthat"),
  lock = FALSE
)
renv::install(".", lock = FALSE)
```

## Hızlı başlangıç

### CLI

```r
library(SigBenchR)

data(datasets_ici)
data(signatures_demo)

res <- run_benchmark(
  signatures = signatures_demo,
  datasets = datasets_ici,
  methods = c("mean", "zscore", "median", "rank"),
  seed = 42L
)

print(res)
plot_auc_heatmap(res, method = "mean")

sm <- summarize_auc(res, method = "mean")
plot_auc_forest(res, signature = "IFNG6", method = "mean", summary = sm)
jaccard_signatures(signatures_demo)
auc_correlation(res, method = "mean")

cat(generate_r_code(
  signature_name = "IFNG6",
  methods = c("mean", "zscore"),
  seed = 42L,
  dataset_names = names(datasets_ici)
))
```

Sentetik `datasets_demo` (`demo_A`, `demo_B`) birim testleri ve hızlı denemeler içindir. Geniş literatür paneli için CLI’da `get_builtin_signatures()` kullanın.

### Shiny

```r
SigBenchR::run_app()
```

Shiny varsayılan kohortları `van_allen_2015`, `hugo_2016`, `riaz_2017`, `gide_2019`’dur. İmza seçimi varsayılan olarak `signatures_demo` (Cytokine-5, IFNG6, Random3) ile sınırlıdır. Null permütasyon kapalıdır; açıkça işaretlenmelidir.

Arayüzde sonuç tablosu, ısı haritası, ROC, DeLong, özet, forest, Jaccard, AUC korelasyonu, HTML rapor indirme ve “Bu analizin R kodunu kopyala” vardır.

## Özellikler

**Skorlama:** `mean`, `zscore` (ortalama gen z-skoru; GSVA değil), `median`, `rank` (kohort-içi sıralama; singscore değil), `ss_rank` (örnek-içi dikey sıralama; singscore normalizasyonuyla uyumlu, saf base R). İşaretli `+`/`-` imzalar `genes_up` / `genes_down` veya isimli sayısal `genes` (`c(IFNG = 1, VEGFA = -1)`) ile verilir.

**Değerlendirme:** ROC-AUC (DeLong %95 GA), PR-AUC (tabakalı bootstrap), isteğe bağlı null permütasyon (`compute_null = FALSE`; `B = 1000`; p = (k+1)/(B+1); BH: `null_p_adj_bh`), eşleştirilmiş DeLong (`compare_signatures`). Skorlama ve permütasyon yalnızca `enc$kept` örneklerinde çalışır.

**Özet ve yardımcıları:** `summarize_auc()` — n-ağırlıklı + ters-varyans, Cochran Q, I², logit-AUC REML + kesilmiş HKSJ; isteğe bağlı `stratify_by`. Forest: `plot_auc_forest()`. SD duyarlılığı: `check_sd_sensitivity()`. Dayanıklılık: `check_sensitivity()` (eşik 0.10 keyfidir). Hanley–McNeil güç: `estimate_auc_power()` / `power_auc70`. Biyolojik eksenler: `correlate_reference_axes()`. İsteğe bağlı Cox PH / Kaplan–Meier: `evaluate_survival()` (`survival` Suggests).

**Veri ve kapsama:** `read_user_dataset()`, `load_cohort_bundle()`, `read_gmt()`. Gen takma adı kurtarma (~53 immuno-onkoloji alias; `match_genes(use_aliases = TRUE)`). `coverage_status`: `ok` / `low_power` / `invalid`. İmza örtüşmesi: `jaccard_signatures()`; kohortlar arası: `auc_correlation()`.

**Rapor:** `export_report()`, `generate_r_code()`. Literatür paneli: `get_builtin_signatures()` (≥8 ICI imzası + 3 negatif kontrol).

## Veri yapısı

```r
list(
  expr = <gen x örnek matris>,
  sample_meta = <data.frame: sample_id, response_label, ...>,
  info = list(normalization = "TPM", scale = "log2")  # veya "linear"
)
```

`scale == "linear"` ise skorlamadan önce `log2(x + 1)` uygulanır.

## Yanıt kodlama

- Varsayılan: `CR` / `PR` = 1, `PD` = 0
- Stabil hastalık (`SD`): varsayılan olarak dışlanır (`sd_handling = "exclude"`)
- Alternatifler: `"as_0"` (ORR: CR/PR vs SD/PD), `"as_1"` (DCR: CR/PR/SD vs PD)

## ICI kohortları

| Kohort | n | n_pos | n_neg | n_SD | Normalizasyon | Ölçek | Kaynak |
|--------|--:|------:|------:|----:|---------------|-------|--------|
| `van_allen_2015` | 42 | 19 | 23 | 0 | TPM | linear | GitHub vanallenlab TPM |
| `hugo_2016` | 26 | 14 | 12 | 0 | FPKM | linear | GEO GSE78220 |
| `riaz_2017` | 49 | 10 | 23 | 16 | FPKM | linear | GEO GSE91061 (`PRCR`→`PR`, `UNK` dışlandı) |
| `gide_2019` | 41 | 19 | 16 | 6 | normalized counts | linear | GitHub Gide/Quek |

Liu 2019 (cBioPortal) S3/LFS 403 verdiği için dahil edilmedi.

**Klinik etiketler**

- `van_allen_2015` (`treatment = "anti-CTLA4"`): Etiketler IMS `mel_van_cli_data.csv` ikili `response` alanından gelir (`1`→PR/pozitif, `0`→PD/negatif; SD yok). Bu, Van Allen (Science 2015) klinik yarar tanımına dayanır (`n_pos = 19`, `n_neg = 23`, `n = 42`) ve Hugo/Riaz/Gide RECIST CR/PR–PD ekseniyle birebir aynı değildir.
- `hugo_2016` / `riaz_2017` / `gide_2019`: anti-PD-1; RECIST CR/PR vs PD. `riaz_2017` 16 SD içerir ve önceki anti-CTLA-4 maruziyeti taşır.
- `gide_2019` (`treatment = "anti-PD1"`): Yalnızca ENA `PD1_*_PRE` (anti-PD1 monoterapi; `n = 41`). `ipiPD1_*_PRE` kombinasyon dahil değildir (`n_pos = 19`, `n_neg = 16`, `n_SD = 6`).

**Ön işleme:** yalnızca pre-treatment; hasta başına tek örnek.

**Paket boyutu:** Her kohortta demo imza genleri + sıfır varyanslı olmayan genlerden `seed = 42` ile ~2.800 gen örneklenir (`compress = "xz"`; `data/` &lt; 5 MB). Kaynak URL’leri ve dışlama kuralları: [`data-raw/PROVENANCE.md`](data-raw/PROVENANCE.md). Yeniden üretim: `data-raw/prepare_ici_datasets.R`, `data-raw/ici_cohort_summary.csv`.

**Sınırlar:** ~2.800 genlik havuzda özel imzalarda kapsama tipik olarak ~%14–15’tir; bu sürüm keşifsel bir prototiptir, tam-genom ICI yeniden analizi değildir. `k ≤ 4` meta-analiz keşifseldir. Kohort normalizasyonları heterojendir; ham ifadeler birleştirilmez. Tam matris için `read_user_dataset()` / `load_cohort_bundle()` kullanın.

**Bağımsızlık:** `Cytokine-5$derived_from` bu kohort adlarını içerir (`independent = FALSE`). `IFNG6` / `Random3` için `derived_from = character(0)` (`independent = TRUE`). Kökeni belirtilmemiş özel imzalar `independent = NA` alır ve varsayılan bağımsız özete girmez.

## Dağıtım (shinyapps.io)

Otomatik deploy yoktur. `rsconnect` ile:

1. Paketi kurun (`devtools::install(".")` veya uzak depodan).
2. Repo kökündeki `app.R` yalnızca `SigBenchR::run_app()` çağırır (`.Rbuildignore` içinde `^app\.R$`).
3. Suggests paketlerini (`shiny`, `DT`, `rmarkdown`, `knitr`) hedef ortamda kurun; aksi halde UI özellikleri sessizce düşer veya rapor indirme hata verir.
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

`rsconnect` bağımlılık listesini doğru çıkarsın diye önce Imports/Suggests paketlerinin kurulabilir olduğundan emin olun. Uygulama dahili melanoma ICI kohortları üzerinde çalışır; yalnızca gen sembolleri girin.

## Atıflar

- Kus M, et al. TCGEx: a powerful visual interface for exploring and analyzing cancer gene expression data. *EMBO Reports* (2025). doi:[10.1038/s44319-025-00407-7](https://doi.org/10.1038/s44319-025-00407-7). Kod: [github.com/atakanekiz/TCGEx](https://github.com/atakanekiz/TCGEx).
- Ayers M, et al. IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade. *J Clin Invest* (2017).
- Robin X, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics* (2011).
- Keilwagen J, et al. PRROC: computing and visualizing precision-recall and receiver operating characteristic curves in R. *Bioinformatics* (2014).

Cytokine-5 imzası TCGEx çalışmasında Pan-ICI verisinden türetilmiştir; yazarlar imzanın üstünlüğünü iddia etmediklerini ve ek doğrulama gerektiğini belirtmiştir. Bu proje bağımsız bir doğrulama aracıdır.

## Lisans

GNU GPL (>= 3). Ayrıntılar için `LICENSE.md`. TCGEx ile uyumludur.
