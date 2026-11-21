# SigBenchR `datasets_ici` veri provenansı

Bu belge yalnızca `data-raw/prepare_ici_datasets.R` ve
`data-raw/ici_cohort_summary.csv` içinde **kodla doğrulanabilen** kaynakları,
filtreleri ve sayıları özetler. Betikte yazılmayan ara örnek sayıları veya
dışlama gerekçeleri uydurulmaz; eksik bilgi açıkça işaretlenir.

Paketleme politikası (tüm kohortlar): `TARGET_N_GENES = 2800`,
`GENE_SEED = 42` (`subsample_neutral_pool`). `.Rbuildignore` içinde
`^data-raw$` aktiftir; bu dizin `R CMD check` tarball’ına girmez.

Klinik sonlanımda paket varsayılanı: `encode_response(..., sd_handling = "exclude")`
(CR/PR = 1, PD = 0, SD dışlanır). Duyarlılık: `as_0` (SD → 0), `as_1` (SD → 1).

---

## Ortak gen havuzu (~2.800 gen) ve kapsama sınırı

`prepare_ici_datasets.R`:

1. Demo imza genlerini zorunlu tutar: Cytokine-5 (`IFNG`, `CXCL9`, `CXCL13`,
   `IL21`, `CCL5`), IFNG6 (`IDO1`, `CXCL9`, `CXCL10`, `HLA-DRA`, `STAT1`,
   `IFNG`), Random3 (`GENA`, `GENB`, `GENC`).
2. Sıfır varyanslı satırları havuzdan çıkarır (`sd > 0`).
3. Kalanlardan `seed = 42` ile rastgele örnekleyerek hedefi ~2800 gene tamamlar.
4. `usethis::use_data(..., compress = "xz")` ile `data/datasets_ici.rda` yazar.

**Kapsama sınırı:** Betik yüzde kapsama hesaplamaz. ~2.800 genlik paket havuzu,
tipik insan transkriptomuna (~18–20k gen) göre yaklaşık **%14–15** gen
kapsamasına karşılık gelir; özel (arbitrary) imzalarda birçok gen paketlenmiş
matriste bulunmayabilir. Bu sürüm keşifsel bir prototiptir; tam-genom ICI
yeniden analizi değildir.

---

## `van_allen_2015`

### Orijinal kaynak / indirme (betikteki URL’ler)

| Rol | Adres |
|-----|--------|
| İfade (TPM) | `https://raw.githubusercontent.com/vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM/master/TPM_RSEM_VAScience2015.txt` |
| Klinik | `https://raw.githubusercontent.com/xmuyulab/ims_gene_signature/main/data/mel_van_cli_data.csv` |
| Sembol eşleme | `https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz` (Ensembl → sembol) |

Repo / kimlik (betik metni): GitHub `vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM`;
klinik `xmuyulab/ims_gene_signature`. Lisans notu: orijinal yayın / dbGaP
`phs000452` (betikteki `license` alanı). GEO/ENA erişim kimliği **bu kohort için
betikte yok**.

### Tedavi / normalizasyon / ölçek

- `treatment = "anti-CTLA4"` (ipilimumab monoterapi; betik sabiti)
- `normalization = "TPM"`, `scale = "linear"`
- `timing = "pre"` (betikte sabit atanır; ayrı on-treatment filtresi yok)

### Örnek seçimi ve sayılar

| Aşama | n | Kaynak |
|-------|--:|--------|
| Ham TPM sütun sayısı | — | kaynak betikte belirtilmemiş |
| Klinik CSV hasta kimliği ile kesişim (`MEL-IPI_` öneki kaldırılır) | — | ara n betikte yazılmamış |
| `response` NA olanlar düşülür | — | ara n betikte yazılmamış |
| **Analiz paketi** | **42** | `ici_cohort_summary.csv`: n_pos=19, n_neg=23, n_excluded_sd=0 |

`finalize_dataset`: `patient_id` yinelenirse ilk örnek tutulur (bu kohortta
tekrarlayan hasta düşümü betikte sayı olarak raporlanmaz).

### Klinik sonlanım

IMS `response` ikili alanı (`map_binary_response`): `1` → `PR`, `0` → `PD`,
`-1` → `SD`. Betik yorumu: IMS paketlemesi Van Allen (Science 2015) **klinik
yarar** (clinical benefit: CR/PR veya uzun sağkalımlı SD → 1; PD veya kısa
sağkalımlı SD → 0) sınıflamasına dayanır. Bu kohortta SD/`-1` yok
(`n_excluded_sd = 0`). RECIST CR/PR–PD ekseni ile birebir değildir.

`sd_handling = "exclude"` sonrası efektif N = 42 (SD yok).

---

## `hugo_2016`

### Orijinal kaynak / indirme

| Rol | Adres |
|-----|--------|
| İfade (FPKM) | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx` |
| Klinik / soft | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/soft/GSE78220_family.soft.gz` |

Erişim: **GEO GSE78220**.

### Tedavi / normalizasyon / ölçek

- `treatment = "anti-PD1"`, `timing = "pre"`
- `normalization = "FPKM"`, `scale = "linear"`

### Örnek seçimi ve sayılar

Kesin dışlama kuralları (betik):

1. Sütun adında `OnTx` veya `.on` geçen örnekler dışlanır (ör. on-treatment /
   `Pt16.OnTx`).
2. `Pt27A.baseline` ve `Pt27B.baseline` birlikte varsa `Pt27B.baseline` düşülür
   (hasta başına tek örnek).
3. Soft’tan `anti-pd-1 response` eşlemesi; COMPLETE→CR, PARTIAL→PR,
   PROGRESSIVE→PD, STABLE→SD; eşleşmeyen / NA etiketler düşülür.
4. `finalize_dataset` ile `patient_id` tekilleştirme.

| Aşama | n | Kaynak |
|-------|--:|--------|
| Ham FPKM sütun sayısı | — | kaynak betikte belirtilmemiş |
| OnTx / Pt27B sonrası | — | ara n betikte yazılmamış |
| **Analiz paketi** | **26** | summary: n_pos=14, n_neg=12, n_excluded_sd=0 |

### Klinik sonlanım

RECIST benzeri soft etiketleri → CR/PR/PD/SD. Bu pakette SD kalmamıştır
(`n_excluded_sd = 0`). `sd_handling = "exclude"` efektif N = 26.

---

## `riaz_2017`

### Orijinal kaynak / indirme

| Rol | Adres |
|-----|--------|
| İfade (FPKM) | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz` |
| Klinik / soft | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/soft/GSE91061_family.soft.gz` |

Erişim: **GEO GSE91061** (BMS038). Entrez → sembol: NCBI `Homo_sapiens.gene_info.gz`.

### Tedavi / normalizasyon / ölçek

- `treatment = "anti-PD1"`, `timing = "pre"`
- `normalization = "FPKM"`, `scale = "linear"`

Önceki anti-CTLA-4 maruziyeti betikte **filtre olarak uygulanmaz**; yalnızca
README/özet metninde heterojenlik notu vardır (kodda dışlama kuralı yok).

### Örnek seçimi ve sayılar

Kesin dışlama kuralları (betik):

1. Yalnız pre: sütun adında `_Pre_` **veya** soft `visit` / biopsy_time `PRE`.
2. Yanıt: `PRCR` → `PR`; `CR`/`PR`/`PD`/`SD` korunur; `UNK`, boş, `NA` → dışlanır.
3. `finalize_dataset` hasta tekilleştirme.

| Aşama | n | Kaynak |
|-------|--:|--------|
| Ham FPKM sütun sayısı | — | kaynak betikte belirtilmemiş |
| Pre alt kümesi | — | ara n betikte yazılmamış |
| UNK/NA düşümü sonrası | — | ara n betikte yazılmamış |
| **Analiz paketi** | **49** | summary: n_pos=10, n_neg=23, n_excluded_sd=16 |

### Klinik sonlanım

RECIST CR/PR vs PD; 16 SD pakette kalır. `sd_handling` politikaları:

| Politika | efektif n | n_pos | n_neg |
|----------|----------:|------:|------:|
| `exclude` (ana) | 33 | 10 | 23 |
| `as_0` | 49 | 10 | 39 |
| `as_1` | 49 | 26 | 23 |

---

## `gide_2019`

### Orijinal kaynak / indirme

| Rol | Adres |
|-----|--------|
| İfade (normalized counts) | `https://raw.githubusercontent.com/miabioinformatics/Gide_Quek_CancerCell2019/master/cancercell_normalized_counts_genenames.txt` |
| Klinik | `https://raw.githubusercontent.com/xmuyulab/ims_gene_signature/main/data/mel_gide19_cli_data.csv` |
| ENA eşleme | `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB23709&result=read_run&fields=run_accession,sample_alias,sample_title&format=tsv` |

Erişim: GitHub `miabioinformatics/Gide_Quek_CancerCell2019`; ENA **PRJEB23709**.

### Tedavi / normalizasyon / ölçek

- `normalization = "normalized_counts"`, `scale = "linear"`
- `sample_meta$treatment = "anti-PD1"` (monoterapi)
- Not: IMS CSV’deki `treatment` sütunu betikte **zamanlama** (`pre`/`on`)
  anlamında kullanılır; paket `sample_meta$treatment` alanı ilaç koludur.

### Örnek seçimi ve sayılar

Kesin dışlama kuralları (betik):

1. Yalnız sütun adı `_PRE$` ile biten örnekler.
2. ENA `sample_alias` → `run_accession`; yalnızca IMS `treatment == "pre"` satırlarıyla
   eşleşenler.
3. Güvenlik: `stopifnot(all(grepl("^PD1_", colnames(...))))` —
   `ipiPD1_*_PRE` (anti-PD1 + anti-CTLA4 kombinasyon) bu alt kümede yoktur /
   dahil edilmez.
4. IMS binary `response`: 1→PR, 0→PD, -1→SD; NA etiketler düşülür.

| Aşama | n | Kaynak |
|-------|--:|--------|
| Ham matris örnek sütunları | — | kaynak betikte belirtilmemiş |
| `_PRE` alt kümesi | — | ara n betikte yazılmamış |
| ENA ∩ IMS pre sonrası | — | ara n betikte yazılmamış; betik yorumu IMS pre n=41 |
| **Analiz paketi** | **41** | summary: n_pos=19, n_neg=16, n_excluded_sd=6 |

### Klinik sonlanım

RECIST CR/PR vs PD (IMS binary üzerinden PR/PD/SD etiketleri). `sd_handling`:

| Politika | efektif n | n_pos | n_neg |
|----------|----------:|------:|------:|
| `exclude` (ana) | 35 | 19 | 16 |
| `as_0` | 41 | 19 | 22 |
| `as_1` | 41 | 25 | 16 |

---

## Dahil edilmeyen kohort

**Liu 2019** (cBioPortal): betik açıkça S3 **403** / erişim engeli nedeniyle
atlanır; Paket A’da dört kohort kalır. URL indirme listesinde yoktur.

---

## Doğrulama özeti (koddan)

| Kohort | n | n_pos | n_neg | n_SD | exclude sonrası N | Norm | Ölçek |
|--------|--:|------:|------:|----:|------------------:|------|-------|
| van_allen_2015 | 42 | 19 | 23 | 0 | 42 | TPM | linear |
| hugo_2016 | 26 | 14 | 12 | 0 | 26 | FPKM | linear |
| riaz_2017 | 49 | 10 | 23 | 16 | 33 | FPKM | linear |
| gide_2019 | 41 | 19 | 16 | 6 | 35 | normalized_counts | linear |

### Bilinçli boşluklar (kaynak betikte belirtilmemiş)

- Her kohort için ham / ara örnek sayıları (filtre öncesi n) betikte
  `message` veya tablo olarak yazılmaz.
- Van Allen için GEO/ENA kimliği yok.
- Riaz’ta önceki CTLA-4 maruziyeti için kod filtresi yok.
- ~%14–15 kapsama oranı betikte hesaplanmaz; 2800 gen hedefinin yaklaşık
  sonucudur.
