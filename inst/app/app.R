library(shiny)
library(survival)
library(splines)
library(ggplot2)

# =============================================================================
# 1. DATA LOADING (relative to app directory)
# =============================================================================

# When the app runs via runApp(), the working directory is inst/app/ itself,
# so we load models from the current directory.
mod_pfs  <- readRDS("FINAL_PFS_RUBIN.rds")
mod_os   <- readRDS("FINAL_OS_RUBIN.rds")
mod_resp <- readRDS("FINAL_RESP_RUBIN_NS13.rds")

# =============================================================================
# 2. PREDICTION ENGINE (WITH ERROR BYPASS)
# =============================================================================

predict_final <- function(model, newdata) {

  # A. PREPARATION
  knots_logsiri <- model$knots
  trms <- delete.response(model$terms)
  environment(trms) <- environment()

  # B. DESIGN MATRIX
  mf <- model.frame(trms, newdata, xlev = model$xlevels, na.action = na.pass)
  X <- model.matrix(trms, mf)

  # C. COEFFICIENT MAPPING
  beta_final <- numeric(ncol(X))
  names(beta_final) <- colnames(X)

  coefs_orig <- model$coef
  src_coef_names <- names(coefs_orig)

  # Function to normalize names (ignores order A:B vs B:A)
  norm_str <- function(s) {
    if (!grepl(":", s)) return(trimws(s))
    paste(sort(trimws(unlist(strsplit(s, ":")))), collapse = ":")
  }

  # Create search map
  coef_map <- list()
  for (nm in src_coef_names) coef_map[[norm_str(nm)]] <- nm

  # Loop through matrix columns
  for (col_nm in colnames(X)) {
    val <- NA

    # 1. Exact search
    if (col_nm %in% src_coef_names) {
      val <- coefs_orig[[col_nm]]
    }
    # 2. Normalized search
    else {
      norm_nm <- norm_str(col_nm)
      if (norm_nm %in% names(coef_map)) {
        real_nm <- coef_map[[norm_nm]]
        val <- coefs_orig[[real_nm]]
      }
    }
    beta_final[col_nm] <- val
  }

  # D. SPLINE FILLING (By position)
  idx_nas <- which(is.na(beta_final))
  if (length(idx_nas) > 0) {
    missing_names <- names(beta_final)[idx_nas]

    # If they are splines, fill by order
    if (all(grepl("ns\\(", missing_names))) {
      idx_spline_coef <- grep("ns\\(", src_coef_names)
      idx_spline_X    <- grep("ns\\(", colnames(X))

      if (length(idx_spline_coef) == length(idx_spline_X)) {
        beta_final[idx_spline_X] <- coefs_orig[idx_spline_coef]
      }
    }
  }

  # --- E. PROBLEM SOLUTION (BYPASS) ---
  if (any(is.na(beta_final))) {
    beta_final[is.na(beta_final)] <- 0
  }

  # F. FINAL CALCULATION
  lp <- as.vector(X %*% beta_final)

  if (!is.null(model$dist) && model$dist == "weibull") {
    # Weibull model
    scale <- model$scale
    log_log_2 <- log(-log(0.5))
    median_val <- exp(lp + scale * log_log_2)

    # Calculate real SE using variance-covariance matrix
    if (!is.null(model$vcov)) {
      var_mat <- model$vcov
      n_coef <- length(model$coef)

      if (nrow(var_mat) > n_coef) {
        var_mat <- var_mat[1:n_coef, 1:n_coef]
      }

      var_mapped <- matrix(0, ncol(X), ncol(X))
      rownames(var_mapped) <- colnames(X)
      colnames(var_mapped) <- colnames(X)

      coef_names <- names(model$coef)
      for (i in seq_along(coef_names)) {
        for (j in seq_along(coef_names)) {
          ni <- coef_names[i]
          nj <- coef_names[j]
          if (ni %in% colnames(X) && nj %in% colnames(X)) {
            var_mapped[ni, nj] <- var_mat[i, j]
          }
        }
      }

      var_mapped[is.na(var_mapped)] <- 0
      se_lp <- sqrt(as.vector(X %*% var_mapped %*% t(X)))
    } else {
      se_lp <- 0.1 * abs(lp)
    }

    ic_low <- exp((lp - 1.96 * se_lp) + scale * log_log_2)
    ic_up  <- exp((lp + 1.96 * se_lp) + scale * log_log_2)

    return(list(val = median_val, lower = ic_low, upper = ic_up, lp = lp, scale = scale, se_lp = se_lp))

  } else {
    # Logistic model
    prob <- 1 / (1 + exp(-lp))
    return(list(val = prob))
  }
}

# =============================================================================
# 2B. PREDICTION ENGINE FOR RESPONSE MODEL WITH SPLINES
# =============================================================================

predict_response <- function(model, newdata) {

  # Extract model components
  coefs <- model$coef

  # Interior and boundary knots (fixed from model specification)
  knots_interior <- c(1, 3)
  boundary_knots <- c(-2.40794560865187, 4.91190045341057)

  # Get logsiri value
  logsiri_val <- newdata$logsiri

  # Create spline basis manually with correct knots
  ns_basis <- splines::ns(logsiri_val, knots = knots_interior, Boundary.knots = boundary_knots)

  # Start linear predictor with intercept
  lp <- coefs["(Intercept)"]

  # Add spline terms (match coefficient names dynamically)
  spline_coef_idx <- sort(grep("^ns\\(logsiri", names(coefs)))
  for (i in seq_along(spline_coef_idx)) {
    lp <- lp + ns_basis[1, i] * coefs[spline_coef_idx[i]]
  }

  # Add diam_low effect
  if (as.character(newdata$diam_low) == "Low") {
    lp <- lp + coefs["diam_lowLow"]
  }

  # Add regimen effect
  regimen <- as.character(newdata$regimen_cat)
  if (regimen == "Gem-Abraxane") {
    lp <- lp + coefs["regimen_catGem-Abraxane"]
    lp <- lp + logsiri_val * coefs["regimen_catGem-Abraxane:logsiri"]
  } else if (regimen == "Gemcitabine mono") {
    lp <- lp + coefs["regimen_catGemcitabine mono"]
    lp <- lp + logsiri_val * coefs["regimen_catGemcitabine mono:logsiri"]
  } else if (regimen == "Other") {
    lp <- lp + coefs["regimen_catOther"]
  } else {
    # FOLFIRINOX is reference, add interaction only
    lp <- lp + logsiri_val * coefs["regimen_catFOLFIRINOX:logsiri"]
  }

  # Add ECOG effect
  ecog <- as.character(newdata$ecog_cat_3)
  if (ecog == "1") {
    lp <- lp + coefs["ecog_cat_31"]
  } else if (ecog == "2+") {
    lp <- lp + coefs["ecog_cat_32+"]
  }

  # Add CACS effect
  if (as.character(newdata$CACS) == "Yes") {
    lp <- lp + coefs["CACSYes"]
  }

  # Convert to probability
  prob <- 1 / (1 + exp(-lp))

  return(list(val = prob))
}

# =============================================================================
# 3. USER INTERFACE (UI)
# =============================================================================

lvls <- mod_pfs$xlevels

# FILTER: Remove "Other" from regimen options
regimen_opts <- lvls$regimen_cat[lvls$regimen_cat != "Other"]

ui <- fluidPage(
  tags$head(tags$style(HTML("
    /* --- GLOBAL BACKGROUND AND TYPOGRAPHY --- */
    body {
      background: linear-gradient(135deg, #e0e4e8 0%, #f3f5f8 100%);
      font-family: 'Segoe UI', system-ui, sans-serif;
      color: #2c3e50;
    }
    /* --- CONTAINERS (SIDEBAR AND PANELS) --- */
    .well {
      background: #ffffff;
      border: 1px solid #d1d9e6;
      box-shadow: 0 10px 20px rgba(0,0,0,0.08);
      border-radius: 12px;
    }
    .shiny-input-container {
      color: #34495e;
    }
    /* --- MAIN BUTTON (BLUE METAL EFFECT) --- */
    .btn-primary {
      background: linear-gradient(to bottom, #3a6186, #26425f);
      border: none;
      box-shadow: 0 4px 6px rgba(0,0,0,0.2);
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 1px;
      transition: all 0.3s ease;
    }
    .btn-primary:hover {
      background: linear-gradient(to bottom, #4b7aa8, #2c4e70);
      box-shadow: 0 6px 10px rgba(0,0,0,0.3);
      transform: translateY(-1px);
    }
    /* --- TITLES AND HEADERS --- */
    h2 {
      color: #2c3e50;
      font-weight: 800;
      text-shadow: 1px 1px 2px #bdc3c7;
      margin-top: 0;
    }
    h3 {
      color: #2c3e50;
      font-weight: 700;
      text-shadow: 1px 1px 2px #bdc3c7;
    }
    h4 {
      color: #5d6d7e;
      font-weight: 600;
      border-bottom: 2px solid #3a6186;
      padding-bottom: 10px;
      margin-bottom: 20px;
    }
    h5 {
      color: #34495e;
      font-weight: bold;
      margin-top: 20px;
    }
    hr {
      border-top: 1px solid #bdc3c7;
    }
    /* --- RESULT CARDS (NEUMORPHISM) --- */
    .result-card {
      background: linear-gradient(145deg, #ffffff, #f0f0f0);
      padding: 20px;
      border-radius: 15px;
      text-align: center;
      box-shadow: 5px 5px 10px #bebebe,
                   -5px -5px 10px #ffffff;
      margin-bottom: 20px;
    }
    .result-value {
      font-size: 32px;
      font-weight: 800;
      color: #26425f;
    }
    .result-label {
      font-size: 14px;
      text-transform: uppercase;
      color: #7f8c8d;
      letter-spacing: 1px;
    }
    .result-ci {
      font-size: 13px;
      color: #95a5a6;
      margin-top: 5px;
    }
    /* --- PLOT CONTAINER --- */
    .plot-container {
      background: linear-gradient(145deg, #ffffff, #f8f9fa);
      border-radius: 15px;
      padding: 15px;
      box-shadow: 5px 5px 10px #bebebe,
                   -5px -5px 10px #ffffff;
      margin-top: 20px;
    }
    .plot-title {
      color: #2c3e50;
      font-weight: 700;
      text-align: center;
      margin-bottom: 10px;
      font-size: 16px;
    }
    .btn-calc { font-weight: bold; font-size: 16px; padding: 12px; }
    /* --- MAIN TITLE STYLING --- */
    .main-title {
      color: #2c3e50;
      font-weight: 800;
      text-shadow: 1px 1px 2px #bdc3c7;
      margin-bottom: 5px;
    }
    .subtitle {
      color: #5d6d7e;
      font-size: 14px;
      font-style: italic;
      margin-bottom: 20px;
    }
    /* --- TOOLTIP ICON STYLING --- */
    .info-icon {
      display: inline-block;
      width: 16px;
      height: 16px;
      background: linear-gradient(to bottom, #3a6186, #26425f);
      color: white;
      border-radius: 50%;
      text-align: center;
      font-size: 11px;
      font-weight: bold;
      line-height: 16px;
      cursor: help;
      margin-left: 5px;
      vertical-align: middle;
    }
    .info-icon:hover {
      background: linear-gradient(to bottom, #4b7aa8, #2c4e70);
    }
    /* --- FOOTER STYLING --- */
    .footer-section {
      margin-top: 30px;
      padding-top: 20px;
      border-top: 1px solid #bdc3c7;
      text-align: center;
    }
    .footer-credit {
      color: #5d6d7e;
      font-size: 14px;
      font-style: italic;
      margin-bottom: 10px;
    }
    .footer-disclaimer {
      color: #95a5a6;
      font-size: 12px;
      font-style: italic;
    }
  "))),

  # Main title with subtitle
  div(class = "main-title", style = "font-size: 24px; text-align: center;",
      "Pancreatic Cancer Prognostic Calculator (Pantheia Group)"),
  div(class = "subtitle", style = "text-align: center;",
      "Illustrating the prognostic impact of SIRI on therapeutic effects"),

  sidebarLayout(
    sidebarPanel(
      h4("Clinical Data"),

      # SIRI with tooltip
      div(
        tags$label(
          "SIRI (absolute value):",
          tags$span(class = "info-icon", title = "Systemic Inflammation Response Index. Formula: SIRI = (Neutrophils x Monocytes) / Lymphocytes. Calculated from absolute blood cell counts.", "?")
        ),
        numericInput("siri", label = NULL, value = 2.1, min = 0.01, step = 0.1)
      ),

      # Tumour burden (sum of RECIST target lesions), three levels
      selectInput("diam", "Tumour burden (sum of RECIST target lesions):",
                  choices = c(">5 cm" = "GT5", "<=5 cm" = "LE5", "Non-measurable disease" = "NonMeasurable"),
                  selected = "GT5"),

      # Regimen options (excluding 'Other')
      selectInput("regimen", "Regimen:", choices = regimen_opts, selected = "FOLFIRINOX"),

      selectInput("ecog", "ECOG PS:", choices = lvls$ecog_cat_3, selected = "1"),

      # CACS with tooltip
      div(
        tags$label(
          "CACS:",
          tags$span(class = "info-icon", title = "Cancer anorexia-cachexia syndrome. In this model, a baseline symptom composite: presence of anorexia, cachexia, asthenia, or weight loss >5%.", "?")
        ),
        selectInput("cacs", label = NULL, choices = lvls$CACS, selected = "Yes")
      ),

      br(),
      actionButton("btn_calc", "CALCULATE NOW", class = "btn-primary btn-block btn-calc")
    ),

    mainPanel(
      conditionalPanel(
        condition = "input.btn_calc > 0",
        h3("Model Estimates"),
        hr(),

        fluidRow(
          # PFS
          column(4,
                 div(class = "result-card",
                     div(class = "result-label", "Median PFS"),
                     div(class = "result-value", textOutput("pfs_val")),
                     div(class = "result-ci", textOutput("pfs_ci"))
                 )
          ),
          # OS
          column(4,
                 div(class = "result-card",
                     div(class = "result-label", "Median OS"),
                     div(class = "result-value", textOutput("os_val")),
                     div(class = "result-ci", textOutput("os_ci"))
                 )
          ),
          # Response
          column(4,
                 div(class = "result-card",
                     div(class = "result-label", "Response Prob."),
                     div(class = "result-value", style = "color:#27ae60", textOutput("resp_val")),
                     div(class = "result-ci", "Clinical Probability")
                 )
          )
        ),

        hr(),

        # --- SURVIVAL PLOTS ---
        h4("Predicted Survival Curves"),
        fluidRow(
          column(6,
                 div(class = "plot-container",
                     div(class = "plot-title", "Progression-Free Survival (PFS)"),
                     plotOutput("plot_pfs", height = "320px")
                 )
          ),
          column(6,
                 div(class = "plot-container",
                     div(class = "plot-title", "Overall Survival (OS)"),
                     plotOutput("plot_os", height = "320px")
                 )
          )
        ),

        hr(),
        div(style = "color:gray; font-size:11px; text-align:right;",
            "Experimental calculator for illustrative purposes only: not for clinical use."),

        # --- FOOTER ---
        div(class = "footer-section",
            div(class = "footer-credit", "Developed by Pantheia/SEOM project researchers"),
            div(class = "footer-disclaimer",
                "Disclaimer: This is an experimental tool designed solely to illustrate the prognostic impact of SIRI on all endpoints. Not intended for clinical decision-making.")
        )
      )
    )
  )
)

# =============================================================================
# 4. SERVER
# =============================================================================

server <- function(input, output) {

  vals <- reactiveValues(pfs = NULL, os = NULL, resp = NULL, base_df = NULL)

  observeEvent(input$btn_calc, {

    # Patient dataframe
    df <- data.frame(
      logsiri     = log(input$siri),
      diam3       = factor(input$diam, levels = c("GT5", "LE5", "NonMeasurable")),
      diam_low    = factor(ifelse(input$diam == "LE5", "Low", "High"), levels = mod_pfs$xlevels$diam_low),
      regimen_cat = factor(input$regimen, levels = mod_pfs$xlevels$regimen_cat),
      ecog_cat_3  = factor(input$ecog, levels = mod_pfs$xlevels$ecog_cat_3),
      CACS        = factor(input$cacs, levels = mod_pfs$xlevels$CACS)
    )

    # Store base df for HR calculations
    vals$base_df <- df

    # Safe calculations
    vals$pfs  <- predict_final(mod_pfs, df)
    vals$os   <- predict_final(mod_os, df)
    vals$resp <- predict_final(mod_resp, df)
  })

  # Text outputs
  output$pfs_val <- renderText({ if (is.null(vals$pfs)) "-" else paste(round(vals$pfs$val, 2), "months") })
  output$pfs_ci  <- renderText({ if (is.null(vals$pfs)) "" else paste0("CI: ", round(vals$pfs$lower, 1), " - ", round(vals$pfs$upper, 1)) })

  output$os_val <- renderText({ if (is.null(vals$os)) "-" else paste(round(vals$os$val, 2), "months") })
  output$os_ci  <- renderText({ if (is.null(vals$os)) "" else paste0("CI: ", round(vals$os$lower, 1), " - ", round(vals$os$upper, 1)) })

  output$resp_val <- renderText({ if (is.null(vals$resp)) "-" else paste0(round(vals$resp$val * 100, 1), "%") })

  # --- PFS PLOT ---
  output$plot_pfs <- renderPlot({
    req(vals$pfs)

    lp <- vals$pfs$lp
    scale <- vals$pfs$scale
    se_lp <- vals$pfs$se_lp

    lp_lo <- lp - 1.96 * se_lp
    lp_hi <- lp + 1.96 * se_lp

    times <- seq(0, 24, by = 0.1)

    # Weibull AFT survival function
    surv_est <- exp(-(times / exp(lp))^(1 / scale))
    surv_lo  <- exp(-(times / exp(lp_hi))^(1 / scale))
    surv_hi  <- exp(-(times / exp(lp_lo))^(1 / scale))

    df_plot <- data.frame(time = times, surv = surv_est, lower = surv_lo, upper = surv_hi)

    # Median text
    med_txt <- paste0("Median: ", round(vals$pfs$val, 1), " months\n(95% CI: ",
                      round(vals$pfs$lower, 1), " - ", round(vals$pfs$upper, 1), ")")

    ggplot(df_plot, aes(x = time)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2E86C1", alpha = 0.2) +
      geom_line(aes(y = surv), color = "#2E86C1", linewidth = 1.5) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "#7f8c8d", linewidth = 0.8) +
      annotate("text", x = 18, y = 0.85, label = med_txt, size = 4.5, color = "#2c3e50", fontface = "bold") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent) +
      scale_x_continuous(breaks = seq(0, 24, 3)) +
      labs(x = "Time (Months)", y = "Survival Probability") +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ecf0f1"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(color = "#34495e", face = "bold", size = 14),
        axis.text = element_text(color = "#5d6d7e", size = 13)
      )
  }, bg = "transparent")

  # --- OS PLOT ---
  output$plot_os <- renderPlot({
    req(vals$os)

    lp <- vals$os$lp
    scale <- vals$os$scale
    se_lp <- vals$os$se_lp

    lp_lo <- lp - 1.96 * se_lp
    lp_hi <- lp + 1.96 * se_lp

    times <- seq(0, 36, by = 0.1)

    # Weibull AFT survival function
    surv_est <- exp(-(times / exp(lp))^(1 / scale))
    surv_lo  <- exp(-(times / exp(lp_hi))^(1 / scale))
    surv_hi  <- exp(-(times / exp(lp_lo))^(1 / scale))

    df_plot <- data.frame(time = times, surv = surv_est, lower = surv_lo, upper = surv_hi)

    # Median text
    med_txt <- paste0("Median: ", round(vals$os$val, 1), " months\n(95% CI: ",
                      round(vals$os$lower, 1), " - ", round(vals$os$upper, 1), ")")

    ggplot(df_plot, aes(x = time)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#C0392B", alpha = 0.2) +
      geom_line(aes(y = surv), color = "#C0392B", linewidth = 1.5) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "#7f8c8d", linewidth = 0.8) +
      annotate("text", x = 27, y = 0.85, label = med_txt, size = 4.5, color = "#2c3e50", fontface = "bold") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent) +
      scale_x_continuous(breaks = seq(0, 36, 3)) +
      labs(x = "Time (Months)", y = "Survival Probability") +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ecf0f1"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(color = "#34495e", face = "bold", size = 14),
        axis.text = element_text(color = "#5d6d7e", size = 13)
      )
  }, bg = "transparent")
}

shinyApp(ui, server)
