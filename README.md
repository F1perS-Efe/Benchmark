# SigBenchR

[![R-CMD-check](https://github.com/F1perS-Efe/Benchmark/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/F1perS-Efe/Benchmark/actions/workflows/R-CMD-check.yaml)

SigBenchR, bir gen imzasının immünoterapi yanıt sınıflarını birden fazla kohortta ayırt edip etmediğini (discrimination) ve temel T-hücresi infiltrasyonu ötesinde bağımsız bilgi katıp katmadığını tekrarlanabilir biçimde sınayan bir R/Bioconductor paketidir. AUC burada klinik tahmin veya kalibrasyon modeli değildir.

**Klinik Sorumluluk Reddi:** Araç klinik karar desteği için değildir. Sonuçlar yalnızca araştırma amaçlıdır; küçük örneklemlerde ($\text{EPV} < 5$ veya $k \le 4$) belirsizlik yüksektir. “Kanıt yok” da geçerli bir bilimsel sonuçtur. Eşikleri sonuçlara baktıktan sonra değiştirmeyin.

Paket zorunlu Bioconductor `Imports` bağımlılığı taşımaz (`LazyData: false`). Veri setleri `list(expr, sample_meta, info)` sözleşmesini izler; isteğe bağlı GMT (`read_gmt()`) ve `SummarizedExperiment` (`as_sigbenchr_dataset()`) adaptörleri `Suggests` üzerinden desteklenir.

## Kurulum

R 4.2 veya üzeri gerekir. Depo kökünde:

```sh
Rscript setup.R
```

Betik `renv` proje kütüphanesini açar, `renv.lock` içindeki paketleri kurar, Shiny arayüzü, HTML rapor ve testler için `Suggests` paketlerini (`shiny`, `DT`, `rmarkdown`, `knitr`, `survival`, `testthat`, `zip`) yükler ve `SigBenchR` paketini aynı kütüphaneye kurar. Ardından:

```r
library(SigBenchR)
SigBenchR::run_app()
```

## Hızlı Başlangıç

### CLI İş Akışı

```r
library(SigBenchR)

data(datasets_ici, package = "SigBenchR")
sigs <- get_builtin_signatures(include_controls = TRUE)

res <- run_benchmark(
  signatures = sigs[c("IFNG6", "T_Cell_Inflamed_GEP", "Cytolytic_Activity", "Myeloid_Suppression")],
  datasets = datasets_ici,
  methods = c("mean", "zscore", "ss_rank"),
  compute_pr = TRUE,
  n_boot = 200L,
  compute_null = TRUE,
  n_perm = 200L,
  match_expression = TRUE,
  seed = 42L
)

# Yöntem-içi (null_p_adj_bh) ve yöntemler-arası global (null_p_adj_bh_global) FDR
head(res[, c(
  "dataset", "signature", "method", "auc", "auc_ci_low", "auc_ci_high",
  "pr_auc", "ci_unreliable", "null_p", "null_p_adj_bh", "null_p_adj_bh_global",
  "background_type"
)])

# Logit-AUC REML + kesilmiş HKSJ meta-analiz özeti
sm <- summarize_auc(res, method = "zscore", stratify_by = c("treatment", "endpoint_type"))

# Temel T-hücresi infiltrasyonu ötesinde bağımsız değer testi (fold-içi CV + yönlü Stouffer Z)
bi <- test_beyond_infiltration(
  datasets = datasets_ici,
  signature = sigs[["T_Cell_Inflamed_GEP"]],
  infiltration_axis = "CD8_T_Cell",
  method = "zscore",
  seed = 42L
)
print(bi$summary)
print(bi$per_cohort[, c(
  "cohort", "n", "epv", "low_epv", "vif", "high_collinearity",
  "delta_auc_combined", "sd_delta_auc", "delta_auc_ci_low", "delta_auc_ci_high",
  "lrt_p", "lrt_p_adj_bh", "adds_independent_value", "analysis_status"
)])
```

### Shiny Arayüzü

```r
SigBenchR::run_app()
# Yerel CSV/TSV kohort yükleme özelliğini açmak için (100 MB yükleme sınırı):
SigBenchR::run_app(allow_upload = TRUE)
```

Shiny arayüzü varsayılan olarak tüm yerleşik literatür imzalarını ve gerçek negatif kontrolleri (`get_builtin_signatures(include_controls = TRUE)`), tek yönlü veya işaretli (`genes_up` / `genes_down`) özel imza girişini ve `van_allen_2015`, `hugo_2016`, `riaz_2017`, `gide_2019` kohortlarını sunar. Tüm sekmeler, grafikler, kopyalanabilir R kodu, HTML rapor (`download_report`), JSON manifest (`download_manifest`) ve Publication Bundle (`download_pub_bundle`) yalnızca son çalıştırılan **değişmez (immutable) `results()` snapshot** nesnesinden beslenir ve aynı deterministik `reproducibility_hash` (MD5) değerini paylaşır.

## Temel Bilimsel ve İstatistiksel Özellikler

1. **Skorlama Motoru ve Sürekli Ağırlıklı İmzalar (`score_signature()`):**
   - Yöntemler: `mean`, `zscore` (vektörleştirilmiş gen z-skoru; GSVA değil), `median`, `rank` (kohort-içi sıralama), `ss_rank` (örnek-içi dikey sıralama; `singscore` normalizasyonuyla uyumlu).
   - **Sürekli Regresyon Katsayıları (`weights`):** `T_Cell_Inflamed_GEP` (Ayers et al., *JCI* 2017; Cristescu et al., *Science* 2018) 18 genlik modelinde 17 pozitif IFN-$\gamma$/adaptif direnç geninin yanında negatif katsayılı **`CD276 = -0.023900`** ve **`STAT1 = 0.250229`** ağırlıkları orijinal doğrusal kombinasyon ($\sum w_i x_{i,j}$) olarak uygulanır. *Not:* RNA-seq üzerinde hesaplanan bu skorlar ticari NanoString nCounter klinik tahlili ile analitik eşdeğerlik iddia etmez.
   - **Biyolojik Eksen Ayrımı:** `Checkpoint_Axis` (9 genlik klasik adaptif checkpoint ekseni: `PDCD1`, `CD274`, `PDCD1LG2`, `CTLA4`, `LAG3`, `HAVCR2`, `TIGIT`, `IDO1`, `VSIR`) ile immün dışlanma (immune exclusion) belirteçleri olan 3 genlik `B7_Immune_Exclusion` (`CD276`, `VTCN1`, `HHLA2`) birbirinden ayrılmıştır. `Myeloid_Suppression` imzasından M1/TAM ayrımı yapmayan pan-makrofaj belirteci `CD68` ve IFN-$\gamma$ hedefi `IDO1` çıkarılmıştır (`direction = "negative"`).

2. **2D İfade-Dispersiyon Eşleşmeli Null Permutation ve Çoklu Test Aileleri (`null_permutation_test()`, `adjust_null_p()`):**
   - `match_expression = TRUE` seçildiğinde rastgele genler yalnızca ortalama ifadeye göre değil, **5 ortalama ifade quantile $\times$ 3 gen standart sapması/dispersiyon tertile** (15 hücreli 2D ızgara) içinden örneklenir.
   - Sürekli ağırlıklı (`weights`) ve işaretli (`genes_up`/`genes_down`) imzaların null çekilişleri, gözlenen imzanın katsayı büyüklük dağılımıyla skorlanır.
   - `background_type` alanı, kullanılan arka plan evrenini açıkça sınıflandırır: `"subsampled_panel_cohort"` (~2.800 gen), `"full_transcriptome"` ($\ge 5.000$ gen) veya `"partial_ortholog_subset"`.
   - `adjust_null_p()` hem yöntem-içi BH FDR (`null_p_adj_bh`) hem de tüm `method x signature x dataset` ailesi üzerinden global BH FDR (`null_p_adj_bh_global`) değerlerini ayrı sütunlarda raporlar.

3. **Fare-İnsan Ortolog Dönüşümü ve Arka Plan Bütünlüğü (`map_mouse_to_human()`):**
   - Polimorfik MHC paralogları, kemirgende bulunmayan genler ve kör `toupper()` dönüşümü varsayılan olarak kapalıdır.
   - `ortholog_policy = "strict_1to1"` (varsayılan) ve `"many_to_1"` politikalarını destekler; `1:many` belirsiz eşleşmeler reddedilir ve `attr(res, "ortholog_stats")` içinde raporlanır.
   - `keep_unmapped_context = TRUE` (varsayılan) iken tam transkriptom arka planı `attr(res, "transcriptome_bg_expr")` içinde korunarak `ss_rank` dikey sıralamasının ve `null_permutation_test()` örnekleme havuzunun eşleşen gen altkümesine daralması önlenir.

4. **Beyond Infiltration (`test_beyond_infiltration()`):**
   - **Fold-İçi Ön-İşleme:** Çapraz doğrulama (CV) sırasında `zscore` ve `rank` gen dönüşümleri ile lojistik regresyon girdisi standardizasyonu (`sig_score_z`, `inf_score_z`) yalnızca eğitim fold'u (`tr_idx`) üzerinde öğrenilip test fold'una (`te_idx`) uygulanır (test sızıntısı önlenir).
   - **Eşleştirilmiş Fold-İçi CV $\Delta\text{AUC}$:** Foldlar arası olasılık birleştirmesi (OOF pooling) yerine fold-içi eşleştirilmiş test AUC farkları üzerinden `delta_auc_combined`, `sd_delta_auc`, `delta_auc_ci_low` ve `delta_auc_ci_high` hesaplanır.
   - **Tanı Korumaları ve Yönlü Stouffer Z:** `adds_independent_value` kararı BH düzeltmeli **`lrt_p_adj_bh < 0.05`** şartına bağlıdır. Küçük örneklem ($\text{EPV} = \min(n_{\text{pos}}, n_{\text{neg}})/2 < 5 \implies \text{low\_epv}$) ve çoklu doğrusallık ($\text{VIF} = 1/(1-r^2) \ge 5 \implies \text{high\_collinearity}$; $\text{VIF} \ge 10 \implies \text{failed\_collinearity}$) açıkça işaretlenir. Kohortlar-arası p-değeri birleştirmesi örneklem büyüklüğü ağırlıklı **yönlü Stouffer Z** (`combined_lrt_p`, `stouffer_z`, `direction_conflict`) ile yapılır; zıt yönlü etkiler birbirini iptal eder.

5. **PR-AUC Bağlı Değer (Tie) ve Küçük Sınıf Korumaları (`pr_auc()`):**
   - Yoğun bağlı skorlarda Davis-Goadrich interpolasyonu sonlu değilse otomatik olarak `auc.integral` tahminine geri dönülür; `tie_fraction`, `small_positive_class` ($n_{\text{pos}} < 5$ veya $n_{\text{neg}} < 5$), `n_boot_valid` ve `ci_unreliable` alanları üretilir.

6. **Sağkalım Hipotez Aileleri (`run_survival_benchmark()`):**
   - Sürekli Cox PH regresyonu (`cox_p` $\to$ `p_adj_bh`) ve medyan-dikotomize Kaplan-Meier log-rank testi (`logrank_p` $\to$ `logrank_p_adj_bh`) ayrı hipotez aileleri olarak bağımsız BH düzeltmesine tabi tutulur.
   - Schoenfeld artıkları orantısal tehlike varsayımı kontrolü (`ph_p`, `ph_violated`) bir model tanı testi olduğundan FDR düzeltmesiyle gizlenmemesi için kasıtlı olarak ham bırakılır.

7. **Raporlama ve Publication Bundle (`export_report()`, `export_publication_bundle()`, `build_reproducibility_manifest()`):**
   - HTML raporu ve `.zip` Publication Bundle (`fig_auc_heatmap.{png,pdf}`, `fig_auc_forest.{png,pdf}`, `fig_roc_curves.{png,pdf}`, tam hassasiyetli CSV tabloları, `analysis_objects.rds`, `manifest.json`, `figure_legends.txt`) ham ve düzeltilmiş p-değerlerini ayrı başlıklar altında ve tüm Faz 4–7 tanı bayraklarıyla birlikte dışa aktarır.

## ICI Kohortları ve Kapsam Sınırlılıkları

| Kohort | n | n_pos | n_neg | n_SD | Normalizasyon | Ölçek | Kaynak |
|--------|--:|------:|------:|----:|---------------|-------|--------|
| `van_allen_2015` | 42 | 19 | 23 | 0 | TPM | linear | GitHub vanallenlab TPM |
| `hugo_2016` | 26 | 14 | 12 | 0 | FPKM | linear | GEO GSE78220 |
| `riaz_2017` | 49 | 10 | 23 | 16 | FPKM | linear | GEO GSE91061 (`PRCR`→`PR`, `UNK` dışlandı) |
| `gide_2019` | 41 | 19 | 16 | 6 | normalized counts | linear | GitHub Gide/Quek |

- **Alt-Örneklenmiş Paket Kohortları vs Tam Transkriptom:** Paket boyutunu <5 MB tutmak için `datasets_ici` kohortları ~2.800 gen içerir (`background_type = "subsampled_panel_cohort"`). Bu nedenle rastgele dış imzalarda gen kapsaması ~%14–15 olabilir. Tam transkriptom ($\ge 5.000$ gen, `background_type = "full_transcriptome"`) analizleri için `read_user_dataset()` veya `load_cohort_bundle()` kullanın.

## Lisans

GNU GPL (>= 3). Ayrıntılar için `LICENSE.md`.
