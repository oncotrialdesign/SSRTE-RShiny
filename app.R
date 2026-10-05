# =============================================================================
# Adaptive Sample Size Explorer — Totality-of-Evidence Approach for
# Rare Disease Clinical Trials (method: Shi, Lin, He, Shu; NEJSDS 2026)
#
# Tabs: Design -> Conditional Power -> Interim -> Final -> Simulation Study
#
# Setup / run instructions are in README.md.
# =============================================================================

library(shiny)
library(ggplot2)
library(bslib)
library(DT)

## ---- Dependency guard -------------------------------------------------------
## The app needs the SSRTE package. If it isn't installed we still let the UI
## load, but computations show a friendly message instead of crashing.
have_SSRTE <- requireNamespace("SSRTE", quietly = TRUE)
if (have_SSRTE) suppressMessages(library(SSRTE))

need_pkg_msg <- paste(
  "The analysis engine is not available in this R session.",
  "Please make sure the required packages are installed (see README),",
  "then restart the app.", sep = "\n"
)

## ---- Small helpers ----------------------------------------------------------

# Parse a user-typed numeric vector like "4,6,4,5,7,6" into c(4,6,4,5,7,6)
parse_vec <- function(txt) {
  v <- suppressWarnings(as.numeric(strsplit(gsub("\\s", "", txt), ",")[[1]]))
  v[!is.na(v)]
}

# Build an exchangeable correlation matrix (1 on diag, rho off-diag)
exch_rho <- function(K, rho) {
  m <- matrix(rho, K, K)
  diag(m) <- 1
  m
}

# The two Sigma structures used throughout the paper's simulations
sigma0_sd  <- function(K) rep(1, K)
sigma0_rho <- function(K) exch_rho(K, 0.3)

sigma1_sd  <- c(0.5, 0.5, 1, 1, 2, 2)
sigma1_rho <- matrix(c(
  1.0, 0.1, 0.3, 0.7, 0.1, 0.3,
  0.1, 1.0, 0.7, 0.1, 0.3, 0.7,
  0.3, 0.7, 1.0, 0.1, 0.3, 0.7,
  0.7, 0.1, 0.1, 1.0, 0.1, 0.3,
  0.1, 0.3, 0.3, 0.1, 1.0, 0.7,
  0.3, 0.7, 0.7, 0.3, 0.7, 1.0
), nrow = 6, byrow = TRUE)

# --- Plain-text report builders --------------------------------------------
# We build the reports directly from the result object's fields rather than
# calling the object's print() method. The print methods rely on cli, whose
# terminal-detection behaves unpredictably inside Shiny's capture.output()
# (it can emit ANSI escapes or almost nothing), which is what left the boxes
# blank. Reading the fields ourselves always yields clean, visible text.

fmt <- function(v, d = 3) {
  if (is.null(v) || length(v) == 0 || is.na(v)) return("NA")
  formatC(v, digits = d, format = "f")
}
rule <- function(title) paste0("\u2500\u2500 ", title, " ",
                               strrep("\u2500", max(2, 40 - nchar(title))))

interim_report <- function(x) {
  e  <- x$estimates; s <- x$ssr; ss <- x$sample_sizes
  s1 <- ss$stage1; s2 <- ss$stage2; sf <- ss$final
  lines <- c(
    rule("Interim analysis"),
    paste0("Decision: ", x$decision),
    paste0("Global test: ", x$method$global_test_type),
    paste0("SSR type: ", x$method$SSR_type),
    if (!is.null(x$method$n_permutations) && !is.na(x$method$n_permutations))
      paste0("Permutations: ", x$method$n_permutations),
    rule("Results"),
    paste0("Interim boundary (C1): ", fmt(x$boundary)),
    paste0("Test statistic (Z): ",   fmt(x$test_statistic_Z)),
    paste0("p-value: ",              fmt(x$p_value)),
    rule("Estimates"),
    paste0("Mean Cohen's d: ",                 fmt(e$mean_cohen_d)),
    paste0("A_obs (SE of bar_t_1): ",          fmt(e$A_obs)),
    paste0("Sum of pairwise correlations: ",   fmt(e$sum_rho)),
    rule("Sample-size re-estimation"),
    paste0("Estimated conditional power: ", fmt(s$estimated_CP)),
    paste0("Promising zone? ", if (isTRUE(s$promising_zone)) "yes" else "no"),
    if (!is.null(s$n_hat_ct) && !is.na(s$n_hat_ct))
      paste0("Estimated control N (before capping): ", s$n_hat_ct),
    paste0("Max allowed control N: ", s$N_ct_max),
    rule("Re-estimated sample sizes"),
    paste0("Stage 1 | Control: ", s1["control"], ", Treatment: ", s1["treatment"]),
    paste0("Stage 2 | Control: ", s2["control"], ", Treatment: ", s2["treatment"]),
    paste0("Final   | Control: ", sf["control"], ", Treatment: ", sf["treatment"])
  )
  paste(lines[!vapply(lines, is.null, logical(1))], collapse = "\n")
}

final_report <- function(x) {
  st1 <- x$stage1; st2 <- x$stage2
  lines <- c(
    rule("Final analysis"),
    paste0("Decision: ", x$decision),
    rule("Combined result"),
    paste0("Final Z: ",            fmt(x$test_statistic_final_Z)),
    paste0("Final p-value: ",      fmt(x$p_value_final)),
    paste0("Final boundary (C2): ",fmt(x$boundary_final)),
    rule("Interim component"),
    paste0("Interim Z (Z1): ",           fmt(st1$Z)),
    paste0("Interim boundary (C1): ",    fmt(x$interim_boundary)),
    paste0("Timing: ",                   fmt(st1$timing)),
    rule("Stage 2 component"),
    paste0("Stage 2 Z (Z2): ",           fmt(st2$Z)),
    paste0("Stage 2 p-value: ",          fmt(st2$p_value_t)),
    paste0("Stage 2 n: control = ", st2$m_ct, ", treatment = ", st2$m_trt),
    rule("Stage 2 estimates"),
    paste0("Mean Cohen's d (stage 2): ", fmt(st2$mean_cohen_d2)),
    paste0("A_obs2 (SE of bar_t_2): ",   fmt(st2$A_obs2)),
    paste0("Sum rho (stage 2): ",        fmt(st2$sum_rho2))
  )
  paste(lines, collapse = "\n")
}

## =============================================================================
## UI
## =============================================================================

ui <- page_navbar(
  title = "Adaptive SSR Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly",
                   primary = "#2C6E91", base_font = font_google("Inter")),
  fillable = FALSE,

  # ---------------------------------------------------------------- About ----
  nav_panel(
    "Overview",
    layout_columns(
      col_widths = c(8, 4),
      card(
        card_header("Adaptive Sample Size — Totality of Evidence"),
        markdown(
"An interactive tool for designing and analysing **rare-disease trials with
several continuous endpoints**, where a single primary endpoint is hard to
pre-specify. It aggregates the mean **Cohen's d** across endpoints into a
single global test and allows a one-time **sample-size re-estimation (SSR)**
at an interim look, within the **promising-zone** framework.

**Two global tests**

- *Exact OLS* — a small-sample t-test tuned for tiny samples; stricter type I
  error control when the sample size is small.
- *Permutation* — nonparametric, fewer assumptions, slightly higher power but
  noisier type I error.

**Two SSR rules**

- *SSR-Power* — re-solve the planning sample-size formula with interim
  estimates; aims at the target power.
- *SSR-CP* — re-solve so the interim **conditional power** hits the target
  1 − β; higher power but larger samples.

Work left-to-right through the tabs: **Design → Conditional Power → Interim →
Final → Simulation Study**."
        )
      ),
      card(
        card_header("Workflow"),
        markdown(
"1. **Design** — get the planned sample size.
2. **Conditional Power** — visualise the promising zone.
3. **Interim** — run the interim test, then SSR if promising.
4. **Final** — combine both stages (inverse-normal method).
5. **Simulation Study** — check type I error / power."
        ),
        if (!have_SSRTE)
          div(class = "alert alert-warning", style = "white-space:pre-wrap;",
              need_pkg_msg)
        else
          div(class = "alert alert-success",
              "Ready \u2014 all tabs are live.")
      )
    ),
    layout_columns(
      col_widths = c(8, 4),
      card(
        card_header("Reference"),
        markdown(
"Lan Shi, Yong Lin, Philip He, Di Shu. *Adaptive Sample Size Using a Totality of
Evidence Approach in Rare Disease Clinical Trials.* **N Engl J STAT DATA SCI**
(2026), 1\u201320. DOI:
[10.51387/26-NEJSDS103](https://doi.org/10.51387/26-NEJSDS103) \u00B7
[Full text](https://nejsds.nestat.org/journal/NEJSDS/article/108/text)"
        )
      ),
      card(
        card_header("Maintainer"),
        markdown(
"**Qin Huang**  
[doloreshqq@outlook.com](mailto:doloreshqq@outlook.com)"
        )
      )
    )
  ),

  # --------------------------------------------------------------- Design ----
  nav_panel(
    "1. Design",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = "Planning inputs (Eq. 2.11)",
        numericInput("d_theta", "Expected mean Cohen's d (\u03B8\u0304*)",
                     value = 0.40, min = 0.01, step = 0.01),
        numericInput("d_K", "Number of endpoints (K)",
                     value = 6, min = 2, step = 1),
        numericInput("d_alpha", "One-sided \u03B1", value = 0.025,
                     min = 0.001, max = 0.2, step = 0.005),
        sliderInput("d_power", "Target power (1 \u2212 \u03B2)",
                    value = 0.80, min = 0.5, max = 0.99, step = 0.01),
        numericInput("d_r", "Allocation ratio (trt:ctrl)",
                     value = 1, min = 0.1, step = 0.1),
        sliderInput("d_timing", "Interim timing (information fraction)",
                    value = 0.50, min = 0.10, max = 0.90, step = 0.05),
        radioButtons("d_rho_mode", "Correlation input",
                     c("Common value" = "common", "Full matrix (paste)" = "matrix"),
                     selected = "common"),
        conditionalPanel(
          "input.d_rho_mode == 'common'",
          numericInput("d_rho", "Common pairwise correlation \u03C1",
                       value = 0.30, min = -0.5, max = 0.99, step = 0.05)
        ),
        conditionalPanel(
          "input.d_rho_mode == 'matrix'",
          helpText("Paste K\u00D7K matrix, rows comma-separated, ",
                   "rows split by semicolons."),
          textAreaInput("d_rho_mat",
            "Correlation matrix",
            value = "1,0.3,0.3,0.3,0.3,0.3; 0.3,1,0.3,0.3,0.3,0.3; 0.3,0.3,1,0.3,0.3,0.3; 0.3,0.3,0.3,1,0.3,0.3; 0.3,0.3,0.3,0.3,1,0.3; 0.3,0.3,0.3,0.3,0.3,1",
            rows = 4)
        ),
        actionButton("d_go", "Compute planned N", class = "btn-primary w-100")
      ),
      layout_columns(
        col_widths = c(6, 6),
        value_box("Planned total N", textOutput("d_Ntot"),
                  showcase = bsicons::bs_icon("people-fill"),
                  theme = "primary"),
        value_box("Interim total N", textOutput("d_Nint"),
                  showcase = bsicons::bs_icon("hourglass-split"),
                  theme = "secondary")
      ),
      card(
        card_header("Full breakdown"),
        verbatimTextOutput("d_out")
      ),
      card(
        card_header("Sample-size sensitivity to the assumed effect size"),
        p("How the planned total N changes if the true mean Cohen's d differs ",
          "from your assumption (all other inputs fixed). The vertical line ",
          "marks your current \u03B8\u0304*."),
        plotOutput("d_plot", height = 320)
      )
    )
  ),

  # ------------------------------------------------------- Conditional Power --
  nav_panel(
    "2. Conditional Power",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = "Promising-zone explorer (Eq. 2.10)",
        p("Reproduces the zoning logic behind Figures 2 & 3: how conditional ",
          "power depends on the interim mean Cohen's d."),
        numericInput("cp_n", "Planned final control N (before SSR)",
                     value = 60, min = 4, step = 1),
        sliderInput("cp_timing", "Interim timing", value = 0.50,
                    min = 0.10, max = 0.90, step = 0.05),
        numericInput("cp_K", "Number of endpoints (K)",
                     value = 6, min = 2, step = 1),
        numericInput("cp_sumrho", "Sum of off-diagonal correlations (\u03A3\u03C1)",
                     value = 9, step = 0.5),
        helpText("For K=6, equal \u03C1=0.3 gives \u03A3\u03C1 = 2\u00D715\u00D70.3 = 9."),
        numericInput("cp_r", "Allocation ratio", value = 1, min = 0.1, step = 0.1),
        sliderInput("cp_power", "Target power (upper zone edge)",
                    value = 0.80, min = 0.5, max = 0.99, step = 0.01),
        sliderInput("cp_min", "CP_min (lower zone edge)",
                    value = 0.20, min = 0.01, max = 0.6, step = 0.01),
        numericInput("cp_mark", "Mark an interim Cohen's d",
                     value = 0.30, step = 0.01)
      ),
      card(
        card_header("Conditional power vs. interim mean Cohen's d"),
        plotOutput("cp_plot", height = 420),
        verbatimTextOutput("cp_point")
      )
    )
  ),

  # -------------------------------------------------------------- Interim ----
  nav_panel(
    "3. Interim",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        title = "Interim analysis + SSR",
        h6("Planned design"),
        numericInput("ia_Nct", "Planned control N (N_ct)", value = 58, min = 4),
        numericInput("ia_Ntrt", "Planned treatment N (N_trt)", value = 58, min = 4),
        numericInput("ia_maxidx", "Max inflation factor", value = 2,
                     min = 1, max = 5, step = 0.5),
        sliderInput("ia_timing", "Interim timing", value = 0.50,
                    min = 0.10, max = 0.90, step = 0.05),
        numericInput("ia_alpha", "One-sided \u03B1", value = 0.025, step = 0.005),
        sliderInput("ia_power", "Target power", value = 0.80,
                    min = 0.5, max = 0.99, step = 0.01),
        numericInput("ia_r", "Allocation ratio", value = 1, min = 0.1, step = 0.1),
        numericInput("ia_K", "Endpoints (K)", value = 6, min = 2, step = 1),
        selectInput("ia_ssr", "SSR type",
                    c("SSR-Power", "SSR-CP", "No SSR")),
        selectInput("ia_test", "Global test",
                    c("exact OLS", "Permutation")),
        conditionalPanel("input.ia_test == 'Permutation'",
          numericInput("ia_nPM", "# permutations", value = 500, min = 50, step = 50),
          numericInput("ia_seedPM", "Permutation seed", value = 1, step = 1)),
        sliderInput("ia_promLL", "Promising-zone lower bound",
                    value = 0.20, min = 0.01, max = 0.6, step = 0.01),
        hr(),
        h6("Stage-1 data (simulated)"),
        p(style="font-size:0.85em;color:#666;",
          "Stage-1 data are simulated with a true mean Cohen's d that you set ",
          "below (e.g. lower than the planning value to mimic an ",
          "over-optimistic design)."),
        numericInput("ia_trueD", "True mean Cohen's d for the data",
                     value = 0.30, step = 0.01),
        radioButtons("ia_sigma", "Covariance structure",
                     c("\u03A30 (equal var/corr)" = "S0",
                       "\u03A31 (6-endpt, unequal)" = "S1"),
                     selected = "S0"),
        conditionalPanel("input.ia_sigma == 'S0'",
          numericInput("ia_rho0", "Common \u03C1 (for \u03A30)", value = 0.3, step = 0.05)),
        textInput("ia_muct_txt", "Control means (comma-separated)",
                  value = "4,6,4,5,7,6"),
        numericInput("ia_seed", "Data seed", value = 42, step = 1),
        actionButton("ia_go", "Run interim analysis", class = "btn-primary w-100")
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box("Decision", textOutput("ia_decision"),
                  showcase = bsicons::bs_icon("signpost-split"), theme = "info"),
        value_box("Interim Z (C1 boundary)", textOutput("ia_z"),
                  showcase = bsicons::bs_icon("graph-up"), theme = "secondary"),
        value_box("Conditional power / zone", textOutput("ia_cp"),
                  showcase = bsicons::bs_icon("bullseye"), theme = "warning")
      ),
      card(
        card_header("Interim report"),
        verbatimTextOutput("ia_out")
      ),
      card(
        card_header("Re-estimated sample sizes"),
        tableOutput("ia_ss")
      )
    )
  ),

  # ---------------------------------------------------------------- Final ----
  nav_panel(
    "4. Final",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        title = "Final analysis",
        p("Run the interim tab first. Stage-2 data are simulated at the ",
          "re-estimated size, then combined via the inverse-normal method ",
          "(Eq. 2.17)."),
        numericInput("fa_trueD", "True mean Cohen's d for stage-2 data",
                     value = 0.30, step = 0.01),
        numericInput("fa_seed", "Stage-2 data seed", value = 58, step = 1),
        actionButton("fa_go", "Run final analysis", class = "btn-primary w-100"),
        hr(),
        uiOutput("fa_status")
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box("Decision", textOutput("fa_decision"),
                  showcase = bsicons::bs_icon("flag-fill"), theme = "success"),
        value_box("Final Z (C2 boundary)", textOutput("fa_z"),
                  showcase = bsicons::bs_icon("graph-up-arrow"), theme = "secondary"),
        value_box("Final p-value", textOutput("fa_p"),
                  showcase = bsicons::bs_icon("percent"), theme = "info")
      ),
      card(
        card_header("Final report"),
        verbatimTextOutput("fa_out")
      )
    )
  ),

  # ------------------------------------------------------- Simulation study --
  nav_panel(
    "5. Simulation Study",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        title = "Type I error / Power",
        p("Type I error \u2192 set true d = 0.   Power \u2192 set true d > 0."),
        numericInput("ss_nsim", "# simulations (nsim)", value = 200,
                     min = 20, max = 5000, step = 20),
        div(class="alert alert-info", style="font-size:0.85em;",
            "Each replicate runs a full two-stage trial. Permutation + large ",
            "nsim can be slow; start small (e.g. 200)."),
        numericInput("ss_trueD", "True mean Cohen's d", value = 0.30, step = 0.01),
        numericInput("ss_planD", "Planning mean Cohen's d",
                     value = 0.386, step = 0.001),
        helpText("If true d < planning d, the trial is 'underestimated'."),
        numericInput("ss_alpha", "One-sided \u03B1", value = 0.025, step = 0.005),
        sliderInput("ss_power", "Target power", value = 0.80,
                    min = 0.5, max = 0.99, step = 0.01),
        sliderInput("ss_timing", "Interim timing", value = 0.50,
                    min = 0.10, max = 0.90, step = 0.05),
        numericInput("ss_r", "Allocation ratio", value = 1, min = 0.1, step = 0.1),
        numericInput("ss_maxidx", "Max inflation factor", value = 2,
                     min = 1, max = 5, step = 0.5),
        selectInput("ss_ssr", "SSR type", c("SSR-Power", "SSR-CP", "No SSR")),
        selectInput("ss_test", "Global test", c("exact OLS", "Permutation")),
        conditionalPanel("input.ss_test == 'Permutation'",
          numericInput("ss_nPM", "# permutations", value = 200, min = 50, step = 50)),
        radioButtons("ss_sigma", "Covariance structure",
                     c("\u03A30 (equal var/corr)" = "S0",
                       "\u03A31 (6-endpt, unequal)" = "S1"),
                     selected = "S0"),
        conditionalPanel("input.ss_sigma == 'S0'",
          numericInput("ss_rho0", "Common \u03C1 (for \u03A30)", value = 0.3, step = 0.05)),
        numericInput("ss_seed1", "Stage-1 base seed", value = 51, step = 1),
        numericInput("ss_seed2", "Stage-2 base seed", value = 82, step = 1),
        actionButton("ss_go", "Run simulation study", class = "btn-primary w-100")
      ),
      layout_columns(
        col_widths = c(6, 6),
        value_box("Rejection proportion",
                  textOutput("ss_rej"),
                  showcase = bsicons::bs_icon("clipboard-data"),
                  theme = "primary"),
        value_box("Interpretation", textOutput("ss_interp"),
                  showcase = bsicons::bs_icon("info-circle"), theme = "secondary")
      ),
      card(
        card_header("Interim-decision breakdown"),
        plotOutput("ss_zoneplot", height = 300)
      ),
      card(
        card_header("Per-replicate interim summary (first 200 rows)"),
        DT::DTOutput("ss_table")
      )
    )
  ),

  nav_spacer(),
  nav_item(tags$span(style = "color:#888;font-size:0.85em;",
                     "Method: Shi, Lin, He & Shu (NEJSDS, 2026)"))
)

## =============================================================================
## SERVER
## =============================================================================

server <- function(input, output, session) {

  # Guard used by every compute path
  require_pkg <- function() {
    validate(need(have_SSRTE, need_pkg_msg))
  }

  # ---- resolve a correlation matrix from the Design tab ---------------------
  design_rho <- reactive({
    K <- input$d_K
    if (input$d_rho_mode == "common") {
      exch_rho(K, input$d_rho)
    } else {
      rows <- strsplit(input$d_rho_mat, ";")[[1]]
      m <- do.call(rbind, lapply(rows, parse_vec))
      validate(need(is.matrix(m) && nrow(m) == K && ncol(m) == K,
                    sprintf("Matrix must be %d\u00D7%d.", K, K)))
      m
    }
  })

  # ----------------------------- 1. DESIGN -----------------------------------
  design_res <- eventReactive(input$d_go, {
    require_pkg()
    rho <- design_rho()
    SSRTE::get_initial_tot_ss(
      beta0    = 1 - input$d_power,
      r        = input$d_r,
      theta_k  = input$d_theta,
      n_endpts = input$d_K,
      alpha0   = input$d_alpha,
      rho      = rho,
      timing0  = input$d_timing
    )
  })

  output$d_Ntot <- renderText({
    r <- design_res(); paste0(r$N_total_actual, " (", r$N_ct, " ctrl / ", r$N_trt, " trt)")
  })
  output$d_Nint <- renderText({
    r <- design_res(); paste0(r$N_interim_actual, "  @ IF=", round(r$timing_actual, 3))
  })
  output$d_out <- renderPrint({
    r <- design_res()
    cat("Planned total (raw formula) N_tot0 :", round(r$N_tot0, 3), "\n")
    cat("Control total N_ct                 :", r$N_ct, "\n")
    cat("Treatment total N_trt              :", r$N_trt, "\n")
    cat("Control @ interim n_ct             :", r$n_ct, "\n")
    cat("Treatment @ interim n_trt          :", r$n_trt, "\n")
    cat("Interim total  N_interim_actual    :", r$N_interim_actual, "\n")
    cat("Final total    N_total_actual      :", r$N_total_actual, "\n")
    cat("Actual timing  timing_actual       :", round(r$timing_actual, 4), "\n")
  })

  output$d_plot <- renderPlot({
    require_pkg()
    rho <- design_rho()
    ds  <- seq(max(0.05, input$d_theta - 0.3), input$d_theta + 0.3, by = 0.01)
    Ntot <- vapply(ds, function(dd) {
      SSRTE::get_initial_tot_ss(
        beta0 = 1 - input$d_power, r = input$d_r, theta_k = dd,
        n_endpts = input$d_K, alpha0 = input$d_alpha,
        rho = rho, timing0 = input$d_timing)$N_total_actual
    }, numeric(1))
    df <- data.frame(d = ds, N = Ntot)
    ggplot(df, aes(d, N)) +
      geom_line(linewidth = 1.1, colour = "#2C6E91") +
      geom_vline(xintercept = input$d_theta, linetype = 2, colour = "#C0392B") +
      annotate("text", x = input$d_theta, y = max(df$N),
               label = "your \u03B8\u0304*", hjust = -0.1, colour = "#C0392B") +
      labs(x = "True mean Cohen's d", y = "Planned total N") +
      theme_minimal(base_size = 14)
  })

  # ------------------------- 2. CONDITIONAL POWER ----------------------------
  output$cp_plot <- renderPlot({
    require_pkg()
    dgrid <- seq(0.0, 0.6, by = 0.005)
    cpv <- vapply(dgrid, function(dd) {
      tryCatch(
        SSRTE::get_CP(timing = input$cp_timing, n = input$cp_n,
                      mean_cohen_d = max(dd, 1e-4), K = input$cp_K,
                      sum_rho = input$cp_sumrho, alloc_rate = input$cp_r),
        error = function(e) NA_real_)
    }, numeric(1))
    df <- data.frame(d = dgrid, cp = cpv)

    lo <- input$cp_min; hi <- input$cp_power
    ggplot(df, aes(d, cp)) +
      # zone shading
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = lo,
               fill = "#BDC3C7", alpha = 0.35) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi,
               fill = "#F7DC6F", alpha = 0.35) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = hi, ymax = Inf,
               fill = "#AED6F1", alpha = 0.35) +
      geom_hline(yintercept = c(lo, hi), linetype = 2, colour = "#C0392B") +
      geom_line(linewidth = 1.2, colour = "#1B4F72") +
      geom_vline(xintercept = input$cp_mark, linetype = 3) +
      annotate("text", x = 0.02, y = lo/2, label = "Unfavorable",
               hjust = 0, size = 4) +
      annotate("text", x = 0.02, y = (lo+hi)/2, label = "Promising \u2192 SSR",
               hjust = 0, size = 4) +
      annotate("text", x = 0.02, y = (hi+1)/2, label = "Favorable",
               hjust = 0, size = 4) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(x = "Interim mean Cohen's d", y = "Conditional power") +
      theme_minimal(base_size = 14)
  })

  output$cp_point <- renderPrint({
    require_pkg()
    cp <- tryCatch(
      SSRTE::get_CP(timing = input$cp_timing, n = input$cp_n,
                    mean_cohen_d = max(input$cp_mark, 1e-4), K = input$cp_K,
                    sum_rho = input$cp_sumrho, alloc_rate = input$cp_r),
      error = function(e) NA_real_)
    bnd <- SSRTE::get_efficacy_boundaries(timing = input$cp_timing,
                                          alpha0 = input$d_alpha %||% 0.025)
    zone <- if (is.na(cp)) "NA"
            else if (cp <= input$cp_min) "Unfavorable (no SSR)"
            else if (cp <  input$cp_power) "Promising (SSR triggered)"
            else "Favorable (no SSR)"
    cat("At interim mean Cohen's d =", input$cp_mark, "\n")
    cat("  Conditional power :", ifelse(is.na(cp), "NA", sprintf("%.1f%%", 100*cp)), "\n")
    cat("  Zone             :", zone, "\n")
    cat("  Efficacy bounds  : C1 =", round(bnd[1],3), " C2 =", round(bnd[2],3), "\n")
  })

  # ------------------------------ shared covariance --------------------------
  make_cov_inputs <- function(sigma, rho0, K, muct_txt) {
    mu <- parse_vec(muct_txt)
    if (sigma == "S1") {
      list(sd = sigma1_sd, rho = sigma1_rho, K = 6,
           mu = if (length(mu) == 6) mu else c(4,6,4,5,7,6))
    } else {
      list(sd = sigma0_sd(K), rho = exch_rho(K, rho0), K = K,
           mu = if (length(mu) == K) mu else rep(5, K))
    }
  }

  # ------------------------------ 3. INTERIM ---------------------------------
  # store interim object + the inputs needed to rebuild stage-2
  ia_store <- reactiveVal(NULL)

  observeEvent(input$ia_go, {
    require_pkg()
    cfg <- if (input$ia_sigma == "S1")
      make_cov_inputs("S1", NA, 6, input$ia_muct_txt)
    else
      make_cov_inputs("S0", input$ia_rho0, input$ia_K, input$ia_muct_txt)

    K <- cfg$K
    n_ct1  <- round(input$ia_Nct  * input$ia_timing)
    n_trt1 <- round(input$ia_Ntrt * input$ia_timing)

    sim1 <- SSRTE::simulate_example_one_stage_data(
      n_trt = n_trt1, n_ct = n_ct1, n_endpts = K,
      exp_mean_cohen_d = input$ia_trueD,
      mu_ct = cfg$mu, sd = cfg$sd, rho = cfg$rho, seed = input$ia_seed)

    if (input$ia_test == "Permutation") set.seed(input$ia_seedPM)

    ia <- tryCatch(
      SSRTE::interim_analysis(
        y_trt1 = sim1$y_trt, y_ct1 = sim1$y_ct,
        N_trt = input$ia_Ntrt, N_ct = input$ia_Nct,
        N_ct_max = ceiling(input$ia_maxidx * input$ia_Nct),
        timing_actual = input$ia_timing,
        alpha0 = input$ia_alpha, beta0 = 1 - input$ia_power,
        alloc_rate = input$ia_r, n_endpts = K,
        SSR_type = input$ia_ssr, global_test_type = input$ia_test,
        nPM = if (input$ia_test == "Permutation") input$ia_nPM else 100,
        promising_LL = input$ia_promLL),
      error = function(e) { showNotification(conditionMessage(e), type="error"); NULL })

    if (is.null(ia)) return()
    ia_store(list(ia = ia, cfg = cfg, K = K, r = input$ia_r))
  })

  output$ia_decision <- renderText({
    s <- ia_store(); if (is.null(s)) return("\u2014")
    if (grepl("Reject", s$ia$decision)) "Reject H0 (stop)" else "Continue"
  })
  output$ia_z <- renderText({
    s <- ia_store(); if (is.null(s)) return("\u2014")
    sprintf("%.3f  (C1 = %.3f)", s$ia$test_statistic_Z, s$ia$boundary)
  })
  output$ia_cp <- renderText({
    s <- ia_store(); if (is.null(s)) return("\u2014")
    cp <- s$ia$ssr$estimated_CP
    zone <- if (isTRUE(s$ia$ssr$promising_zone)) "Promising" else "not promising"
    if (is.na(cp)) "\u2014 (stopped)" else sprintf("%.1f%% (%s)", 100*cp, zone)
  })
  output$ia_out <- renderPrint({
    s <- ia_store()
    validate(need(!is.null(s), "Run the interim analysis to see the report."))
    cat(interim_report(s$ia))
  })
  output$ia_ss <- renderTable({
    s <- ia_store(); validate(need(!is.null(s), ""))
    ss <- s$ia$sample_sizes
    data.frame(
      Stage = c("Stage 1", "Stage 2 (new)", "Final total"),
      Control   = c(ss$stage1["control"],  ss$stage2["control"],  ss$final["control"]),
      Treatment = c(ss$stage1["treatment"],ss$stage2["treatment"],ss$final["treatment"]),
      check.names = FALSE)
  }, striped = TRUE, bordered = TRUE)

  # ------------------------------- 4. FINAL ----------------------------------
  fa_store <- reactiveVal(NULL)

  output$fa_status <- renderUI({
    s <- ia_store()
    if (is.null(s))
      div(class="alert alert-warning", "No interim result yet \u2014 run tab 3 first.")
    else if (grepl("Reject", s$ia$decision))
      div(class="alert alert-danger",
          "Interim stopped for efficacy \u2014 there is no stage 2 / final step.")
    else {
      s2 <- s$ia$sample_sizes$stage2
      div(class="alert alert-info",
          sprintf("Interim continued. Stage-2 to simulate: %d ctrl / %d trt.",
                  s2["control"], s2["treatment"]))
    }
  })

  observeEvent(input$fa_go, {
    require_pkg()
    s <- ia_store()
    validate(need(!is.null(s), "Run the interim analysis first."))
    if (grepl("Reject", s$ia$decision)) {
      showNotification("Trial stopped at interim; no final analysis needed.",
                       type = "warning"); return()
    }
    s2 <- s$ia$sample_sizes$stage2
    n_ct2 <- as.integer(s2["control"]); n_trt2 <- as.integer(s2["treatment"])
    validate(need(n_ct2 > 1 && n_trt2 > 1, "Stage-2 size is too small."))

    sim2 <- SSRTE::simulate_example_one_stage_data(
      n_trt = n_trt2, n_ct = n_ct2, n_endpts = s$K,
      exp_mean_cohen_d = input$fa_trueD,
      mu_ct = s$cfg$mu, sd = s$cfg$sd, rho = s$cfg$rho, seed = input$fa_seed)

    fa <- tryCatch(
      SSRTE::final_analysis(interim = s$ia, y_trt2 = sim2$y_trt,
                            y_ct2 = sim2$y_ct, alloc_rate = s$r),
      error = function(e){ showNotification(conditionMessage(e),type="error"); NULL })
    if (is.null(fa)) return()
    fa_store(fa)
  })

  output$fa_decision <- renderText({
    fa <- fa_store(); if (is.null(fa)) return("\u2014")
    if (grepl("Reject", fa$decision)) "Reject H0" else "Fail to reject"
  })
  output$fa_z <- renderText({
    fa <- fa_store(); if (is.null(fa)) return("\u2014")
    sprintf("%.3f  (C2 = %.3f)", fa$test_statistic_final_Z, fa$boundary_final)
  })
  output$fa_p <- renderText({
    fa <- fa_store(); if (is.null(fa)) return("\u2014")
    sprintf("%.4f", fa$p_value_final)
  })
  output$fa_out <- renderPrint({
    fa <- fa_store()
    validate(need(!is.null(fa), "Run the final analysis to see the report."))
    cat(final_report(fa))
  })

  # -------------------------- 5. SIMULATION STUDY ----------------------------
  ss_store <- reactiveVal(NULL)

  observeEvent(input$ss_go, {
    require_pkg()
    cfg <- if (input$ss_sigma == "S1")
      list(sd = sigma1_sd, rho = sigma1_rho, K = 6, mu = c(4,6,4,5,7,6))
    else
      list(sd = sigma0_sd(6), rho = exch_rho(6, input$ss_rho0), K = 6,
           mu = c(4,6,4,5,7,6))

    withProgress(message = "Running simulation study...", value = 0.3, {
      res <- tryCatch(
        SSRTE::SSRTE_simstudy(
          nsim = input$ss_nsim, alpha = input$ss_alpha,
          beta = 1 - input$ss_power, n_endpts = cfg$K,
          allocate_rate = input$ss_r,
          mean_cohen_d_truth = input$ss_trueD,
          mu_ct = cfg$mu, sd = cfg$sd, rho = cfg$rho,
          timing_interim = input$ss_timing,
          initial_tot_ss = TRUE,
          exp_mean_cohen_d_at_planning = input$ss_planD,
          SSR_type = input$ss_ssr, global_test_type = input$ss_test,
          max_ss_index = input$ss_maxidx,
          seed1 = input$ss_seed1, seed2 = input$ss_seed2,
          nPM = if (input$ss_test == "Permutation") input$ss_nPM else 100),
        error = function(e){ showNotification(conditionMessage(e),type="error"); NULL })
      incProgress(0.7)
      if (!is.null(res)) ss_store(list(res = res, trueD = input$ss_trueD))
    })
  })

  output$ss_rej <- renderText({
    s <- ss_store(); if (is.null(s)) return("\u2014")
    sprintf("%.2f%%", 100 * s$res$proportion_of_rejection)
  })
  output$ss_interp <- renderText({
    s <- ss_store(); if (is.null(s)) return("\u2014")
    if (s$trueD == 0) "Type I error (target \u2264 2.5%)" else "Power (target ~80%)"
  })

  output$ss_zoneplot <- renderPlot({
    s <- ss_store()
    validate(need(!is.null(s), "Run the simulation study to see results."))
    ia <- s$res$interim_analysis
    validate(need(nrow(ia) > 0, "No interim rows returned."))
    # classify each replicate's interim outcome
    stopped <- grepl("Reject", ia$decision)
    zone <- ifelse(stopped, "Reject / Stop",
             ifelse(isTRUE_vec(ia$promising_zone), "Cont. Promising",
              ifelse(!is.na(ia$est_CP) & ia$est_CP >= (1 - (1-input$ss_power)),
                     "Cont. Favorable", "Cont. Unfavorable")))
    df <- as.data.frame(table(Zone = zone))
    df$Pct <- 100 * df$Freq / sum(df$Freq)
    ggplot(df, aes(reorder(Zone, -Pct), Pct, fill = Zone)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = sprintf("%.1f%%", Pct)), vjust = -0.3) +
      labs(x = NULL, y = "Percent of replicates") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  })

  output$ss_table <- DT::renderDT({
    s <- ss_store()
    validate(need(!is.null(s), "Run the simulation study to see the table."))
    ia <- s$res$interim_analysis
    keep <- intersect(c("sim_id","decision","test_statistic_Z","p_value",
                        "mean_cohen_d","sum_rho","est_CP","promising_zone",
                        "n_hat_ct","n_ct_final","n_trt_final"), names(ia))
    DT::datatable(head(ia[keep], 200), rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE)) |>
      DT::formatRound(intersect(c("test_statistic_Z","p_value","mean_cohen_d",
                                  "sum_rho","est_CP"), keep), 3)
  })
}

# helper for vectorized isTRUE over a logical column that may contain NA
isTRUE_vec <- function(x) !is.na(x) & x == TRUE
`%||%` <- function(a, b) if (is.null(a)) b else a

shinyApp(ui, server)
