# TODO — mod_diffanalysis: Run-on-button-click + mergeList staleness fix

## Steps

- [x] **Step 1 — Add `snap` reactiveValues + button-click observer**
  Inserted right after the `mtcoderDA <- eventReactive(...)` definition.
  - Snapshots `input$diff_factor`, `input$Cond1`, `input$Cond2`, `isNumFactor()` and `local_physeq()` at button-click time.
  - Explicitly calls `deseqDA()`, `mgSeqDA()`, `mtcoderDA()` from the observer to force evaluation regardless of which tab is currently visible.

- [x] **Step 2 — `output$deseqTab`: use `snap$physeq`**
  Replaced `local_physeq()` with `snap$physeq` so the tax_table / refseq joined to DESeq results matches the dataset that was actually analyzed.

- [x] **Step 3 — `mergeList()`: replace live contrast inputs with snapshots**
  - `input$diff_factor` → `snap$factor`
  - `input$Cond1` → `snap$cond1`
  - `input$Cond2` → `snap$cond2`
  - `isNumFactor()` → `snap$is_num`
  - `input$pval` left reactive (so users can re-threshold without re-running DA).

- [x] **Step 4 — `reacbarplot1()`: replace live contrast inputs with snapshots**
  - `input$diff_factor` → `snap$factor`
  - `input$Cond1` → `snap$cond1`
  - `input$Cond2` → `snap$cond2`
  - `input$Nmeth`, `input$minAb`, `input$Nfeat` left reactive (plot params).

- [ ] **Step 5 — Hand off to user for testing.**
