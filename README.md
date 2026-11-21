# SigBenchR

Bir gen imzasının (gen listesi) immünoterapi yanıt sınıfları arasında **ayırt edici** (discrimination) olup olmadığını birden fazla veri setinde istatistiksel olarak sınayan, tekrarlanabilir doğrulama aracı. AUC klinik tahmin / kalibrasyon modeli değildir.

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

### CLI (gerçek ICI kohortları)

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

# Tekrarlanabilir R kodu / HTML rapor
cat(generate_r_code(
  signature_name = "IFNG6",
  methods = c("mean", "zscore"),
  seed = 42L,
  dataset_names = names(datasets_ici)
))
```

Sentetik `datasets_demo` (`demo_A`, `demo_B`) birim testleri ve hızlı denemeler için durur.

### Shiny arayüzü

```r
SigBenchR::run_app()
```

Shiny varsayılanı gerçek ICI kohortlarıdır (`van_allen_2015`, `hugo_2016`, `riaz_2017`, `gide_2019`). Arayüzde sonuç tablosu, ısı haritası, ROC, DeLong, ağırlıklı özet, forest, Jaccard, AUC korelasyonu, **HTML rapor indirme** ve **“Bu analizin R kodunu kopyala”** bulunur.

## Özellik özeti

- Skorlama: `mean`, `zscore` (basit ortalama gen z-skoru; GSVA değil), `median`, `rank` (kohort-içi göreli sıralama; singscore değil), **`ss_rank`** (kohorttan bağımsız dikey tek-örnek sıralama; singscore tek/çift yönlü normalizasyonla uyumlu saf base R)
- İşaretli (`+`/`-`) imzalar: `genes_up` / `genes_down` veya isimli sayısal `genes` (`c(IFNG = 1, VEGFA = -1)`)
- ROC-AUC (ayırt edicilik) + DeLong %95 GA / varyans; küçük kohort uyarısı
- PR-AUC + prevalans + tabakalı bootstrap GA (`PRROC`)
- Null permütasyon (`auc_mw` vektörize; p = (k+1)/(B+1); işaretli null için `n_genes_down`); `null_p_adj_bh`
- Eşleştirilmiş DeLong iki imza karşılaştırması (`compare_signatures`, BH)
- Ağırlıklı özet: `summarize_auc()` — n-ağırlıklı + ters-varyans, Cochran `Q`, `I²`; isteğe bağlı **`stratify_by`** (tedavi / klinik sonlanım alt grupları)
- İmza örtüşmesi: `jaccard_signatures()`; veri setleri boyunca `auc_correlation()`
- Forest grafik: `plot_auc_forest()`
- Gen kapsama denetimi + gri düşük-kapsama ısı haritası
- Skorlama dayanıklılığı: `check_sensitivity()` (eşik 0.10 keyfidir)
- HTML rapor: `export_report()` (`inst/rmd/benchmark_report.Rmd`)
- Tekrarlanabilir betik: `generate_r_code()`
- Gerçek melanoma ICI kohortları (`datasets_ici`) + sentetik demo + Shiny

## Faz 3B Paket A — gerçek kohortlar

| Kohort | n | n_pos | n_neg | n_SD | Normalizasyon | Ölçek | Kaynak |
|--------|--:|------:|------:|----:|---------------|-------|--------|
| `van_allen_2015` | 42 | 19 | 23 | 0 | TPM | linear | GitHub vanallenlab TPM |
| `hugo_2016` | 26 | 14 | 12 | 0 | FPKM | linear | GEO GSE78220 |
| `riaz_2017` | 49 | 10 | 23 | 16 | FPKM | linear | GEO GSE91061 (`PRCR`→`PR`, `UNK` dışlandı) |
| `gide_2019` | 41 | 19 | 16 | 6 | normalized counts | linear | GitHub Gide/Quek |

Liu 2019 (cBioPortal) S3/LFS **403** verdiği için dahil edilmedi.

**Klinik etiket şeffaflığı**

- `van_allen_2015` (`sample_meta$treatment = "anti-CTLA4"`): Yanıt etiketleri `xmuyulab/ims_gene_signature` `mel_van_cli_data.csv` ikili `response` alanından gelir (`1`→PR/pozitif, `0`→PD/negatif; bu kohortta SD/`-1` yok). IMS paketlemesi Van Allen (Science 2015) **klinik yarar** (clinical benefit: CR/PR veya uzun sağkalımlı SD) sınıflamasına dayanır → **`n_pos = 19`**, `n_neg = 23`, `n = 42`. Bu tanım Hugo/Riaz/Gide RECIST CR/PR–PD ekseni ile birebir değildir.
- `hugo_2016` / `riaz_2017` / `gide_2019`: anti-PD-1; RECIST CR/PR vs PD. `riaz_2017` **16 SD** içerir ve önceki anti-CTLA-4 maruziyeti heterojenliği taşır.
- `gide_2019` (`sample_meta$treatment = "anti-PD1"`): `n = 41` alt kümesi yalnızca ENA `PD1_*_PRE` örnekleridir (**anti-PD1 monoterapi**); normalizasyon **normalized counts**. `ipiPD1_*_PRE` (anti-PD1 + anti-CTLA4 kombinasyon) dahil değildir; IMS `mel_gide19_cli_data.csv` pre satırları ile birebir örtüşür (`n_pos = 19`, `n_neg = 16`, `n_SD = 6`).

**Ön işleme:** yalnızca pre-treatment; hasta başına tek örnek (pseudo-replication yok).  
**Tarafsız null havuzu / paket boyutu:** her kohortta demo imza genleri + sıfır varyanslı olmayan genlerden **`seed = 42`** ile rastgele **~2.800 gen**; `compress = "xz"`; `data/` toplamı **&lt; 5 MB**. Ayrıntı, kaynak URL’leri ve dışlama kuralları: [`data-raw/PROVENANCE.md`](data-raw/PROVENANCE.md); yeniden üretim: `data-raw/prepare_ici_datasets.R`, `data-raw/ici_cohort_summary.csv`.  
**Dürüst sınırlama:** ~2.800 genlik havuz nedeniyle özel (arbitrary) imzalarda gen kapsaması tipik olarak ~%14–15 ile sınırlıdır; bu sürüm **keşifsel bir prototiptir**, tam-genom ICI yeniden analizi değildir.  
**Bağımsızlık:** `Cytokine-5$derived_from` bu kohort adlarını içerir (`independent = FALSE`); `IFNG6` / `Random3` açıkça `derived_from = character(0)` taşır (`independent = TRUE`). Kökeni belirtilmemiş özel imzalar `independent = NA` alır ve varsayılan bağımsız özete girmez.  
**Yazar:** `DESCRIPTION` içinde `Efe Kaan Şengönül <fperskaan@gmail.com>` (`aut`, `cre`); placeholder yok.

## SigBenchR v1 kapsam dışı / sonraki sprintler

- Takma ad (alias) eşleştirme **yok**
- DerSimonian–Laird rastgele etkiler özeti **yok**
- Bioconductor bağımlılığı (`singscore` / `GSVA` / `SummarizedExperiment`) **yok** (`ss_rank` saf base R)
- Shiny UX yenilemesi sonraki sprintte

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

Not: `rsconnect` paket listesini doğru çıkarsın diye önce `SigBenchR`’ın Imports/Suggests bağımlılıklarının kurulabilir olduğundan emin olun. Uygulama dahili halka açık melanoma ICI kohortları üzerinde çalışır; yalnızca gen sembolleri girin.

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
