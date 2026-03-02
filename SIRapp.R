# ============================================================
#  SIR / SEIR / SIRS Infectious Disease Model — Shiny App
# ============================================================
# install.packages(c("shiny", "deSolve", "ggplot2", "patchwork", "scales"))
# ============================================================

library(shiny)
library(deSolve)
library(ggplot2)
library(patchwork)
library(scales)

# ================================================================
#  Simulation helpers
# ================================================================

# ---- Continuous ODE solver ----------------------------------
run_continuous <- function(model, N, S0, I0, E0, beta, gamma, sigma,
                           omega, mu, nu, tmax) {
  I0 <- min(I0, N); E0 <- min(E0, N - I0)
  R0_init <- N * nu          # initially vaccinated → R
  S0_eff  <- max(N - I0 - E0 - R0_init, 1)
  
  state <- switch(model,
                  "SIR"  = c(S=S0_eff, I=I0, R=R0_init),
                  "SEIR" = c(S=S0_eff, E=E0, I=I0, R=R0_init),
                  "SIRS" = c(S=S0_eff, I=I0, R=R0_init)
  )
  
  times <- seq(0, tmax, by = 0.2)
  
  ode_fn <- function(t, state, parms) {
    with(as.list(c(state, parms)), {
      Stot <- if ("E" %in% names(state)) S + E + I + R else S + I + R
      switch(model,
             "SIR" = {
               dS <- -beta * S * I / Stot
               dI <-  beta * S * I / Stot - gamma * I
               dR <-  gamma * I
               list(c(dS, dI, dR))
             },
             "SEIR" = {
               dS <- -beta * S * I / Stot
               dE <-  beta * S * I / Stot - sigma * E
               dI <-  sigma * E            - gamma * I
               dR <-  gamma * I
               list(c(dS, dE, dI, dR))
             },
             "SIRS" = {
               dS <- -beta * S * I / Stot + omega * R
               dI <-  beta * S * I / Stot - gamma * I
               dR <-  gamma * I           - omega * R
               list(c(dS, dI, dR))
             }
      )
    })
  }
  
  tryCatch(
    as.data.frame(ode(y=state, times=times, func=ode_fn,
                      parms=list(), method="lsoda")),
    error = function(e) NULL
  )
}

# ---- Discrete difference equations -------------------------
run_discrete <- function(model, N, S0, I0, E0, beta, gamma, sigma,
                         omega, mu, nu, tmax) {
  steps   <- as.integer(tmax)
  R0_init <- round(N * nu)
  S_init  <- max(N - I0 - round(E0) - R0_init, 1)
  
  results <- switch(model,
                    "SIR"  = { mat <- matrix(0, steps+1, 3)
                    colnames(mat) <- c("S","I","R")
                    mat[1,] <- c(S_init, I0, R0_init); mat },
                    "SEIR" = { mat <- matrix(0, steps+1, 4)
                    colnames(mat) <- c("S","E","I","R")
                    mat[1,] <- c(S_init, E0, I0, R0_init); mat },
                    "SIRS" = { mat <- matrix(0, steps+1, 3)
                    colnames(mat) <- c("S","I","R")
                    mat[1,] <- c(S_init, I0, R0_init); mat }
  )
  
  for (t in seq_len(steps)) {
    r <- results[t, ]
    Stot <- sum(r)
    S <- max(r["S"], 0); I <- max(r["I"], 0)
    E <- if ("E" %in% colnames(results)) max(r["E"], 0) else 0
    R <- max(r["R"], 0)
    
    new_inf  <- beta  * S * I / max(Stot, 1)
    new_rec  <- gamma * I
    new_loss <- if (model == "SIRS") omega * R else 0
    new_exp  <- if (model == "SEIR") sigma * E else 0
    
    switch(model,
           "SIR" = {
             results[t+1, "S"] <- max(S - new_inf, 0)
             results[t+1, "I"] <- max(I + new_inf - new_rec, 0)
             results[t+1, "R"] <- R + new_rec
           },
           "SEIR" = {
             results[t+1, "S"] <- max(S - new_inf, 0)
             results[t+1, "E"] <- max(E + new_inf - new_exp, 0)
             results[t+1, "I"] <- max(I + new_exp - new_rec, 0)
             results[t+1, "R"] <- R + new_rec
           },
           "SIRS" = {
             results[t+1, "S"] <- max(S - new_inf + new_loss, 0)
             results[t+1, "I"] <- max(I + new_inf - new_rec,  0)
             results[t+1, "R"] <- max(R + new_rec - new_loss, 0)
           }
    )
  }
  
  df      <- as.data.frame(results)
  df$time <- 0:steps
  df
}

# ================================================================
#  UI
# ================================================================
ui <- fluidPage(
  titlePanel("🦠 SIR Infectious Disease Model"),
  
  tags$head(tags$style(HTML("
    body { background:#f4f6f8; font-family:'Segoe UI', sans-serif; }
    .well { background:#fff; border:1px solid #dee2e6; border-radius:8px; }
    .section-header {
      color:#444; font-size:11px; font-weight:bold;
      text-transform:uppercase; letter-spacing:1px;
      margin-top:14px; margin-bottom:3px;
      border-bottom:1px solid #e0e0e0; padding-bottom:3px;
    }
    .info-box {
      padding:9px 12px; border-radius:5px; font-size:12.5px;
      line-height:1.75; margin-top:8px;
    }
    .r0-high   { background:#fdecea; border-left:4px solid #f44336; color:#b71c1c; }
    .r0-medium { background:#fff8e1; border-left:4px solid #ffc107; color:#7b5800; }
    .r0-low    { background:#e8f5e9; border-left:4px solid #4caf50; color:#1b5e20; }
    .mode-badge {
      display:inline-block; padding:2px 10px; border-radius:12px;
      font-size:11.5px; font-weight:700; margin-left:6px; vertical-align:middle;
    }
    .badge-cont { background:#d0f0ff; color:#0077aa; }
    .badge-disc { background:#ffe0b2; color:#bf360c; }
    .stat-grid {
      display:grid; grid-template-columns:1fr 1fr; gap:6px; margin-top:8px;
    }
    .stat-card {
      background:#f9f9f9; border:1px solid #e0e0e0; border-radius:6px;
      padding:7px 10px; font-size:12px; text-align:center;
    }
    .stat-card b { display:block; font-size:15px; color:#2c3e50; }
  "))),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # ---- Model selection ----------------------------------
      div(class="section-header", "🔬 Model"),
      selectInput("model", NULL,
                  choices  = c("SIR", "SEIR", "SIRS"),
                  selected = "SIR"),
      uiOutput("model_desc"),
      
      # ---- Time formulation ---------------------------------
      div(class="section-header", "⏱️ Time Formulation"),
      radioButtons("time_mode", NULL,
                   choices  = c("Continuous (ODE)"="continuous",
                                "Discrete (Difference Eq.)"="discrete"),
                   selected = "continuous", inline=FALSE),
      
      hr(),
      # ---- Population ---------------------------------------
      div(class="section-header", "👥 Population"),
      numericInput("N",  "Total population (N)",
                   value=10000, min=100, max=1e7, step=100),
      sliderInput("I0", "Initially infected (I₀)",
                  min=1, max=500, value=10, step=1),
      
      # SEIR: exposed compartment
      conditionalPanel("input.model == 'SEIR'",
                       sliderInput("E0", "Initially exposed (E₀)",
                                   min=0, max=500, value=0, step=1)
      ),
      
      # Vaccination
      sliderInput("nu", "Vaccination coverage (ν, fraction immune at t=0)",
                  min=0, max=0.99, value=0, step=0.01),
      
      hr(),
      # ---- Transmission parameters --------------------------
      div(class="section-header", "⚙️ Transmission Parameters"),
      sliderInput("beta",  "β — Transmission rate (contacts × prob.)",
                  min=0.01, max=2.0, value=0.3, step=0.01),
      sliderInput("gamma", "γ — Recovery rate (1/γ = infectious period)",
                  min=0.01, max=1.0, value=0.1, step=0.01),
      
      # SEIR: incubation rate
      conditionalPanel("input.model == 'SEIR'",
                       sliderInput("sigma", "σ — Incubation rate (1/σ = incubation period)",
                                   min=0.01, max=1.0, value=0.2, step=0.01)
      ),
      
      # SIRS: waning immunity
      conditionalPanel("input.model == 'SIRS'",
                       sliderInput("omega", "ω — Waning immunity rate (1/ω = immunity duration)",
                                   min=0.001, max=0.1, value=0.01, step=0.001)
      ),
      
      hr(),
      # ---- Simulation ---------------------------------------
      div(class="section-header", "📅 Simulation"),
      sliderInput("tmax", "Duration (days)",
                  min=30, max=1000, value=200, step=10),
      
      hr(),
      # ---- Key metrics --------------------------------------
      div(class="section-header", "📊 Key Metrics"),
      uiOutput("r0_box"),
      uiOutput("stat_cards")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("📈 Epidemic Curves",
                 br(),
                 plotOutput("epi_curve_plot", height="440px")
        ),
        tabPanel("🥧 Compartment Dynamics",
                 br(),
                 plotOutput("area_plot", height="440px")
        ),
        tabPanel("📉 Log Scale",
                 br(),
                 plotOutput("log_plot", height="440px")
        ),
        tabPanel("🔄 Side-by-Side (Cont. vs Disc.)",
                 br(),
                 plotOutput("comparison_plot", height="540px")
        ),
        tabPanel("ℹ️ Theory",
                 br(),
                 uiOutput("theory_panel")
        )
      )
    )
  )
)

# ================================================================
#  Server
# ================================================================
server <- function(input, output, session) {
  
  # ---- Model description ------------------------------------
  output$model_desc <- renderUI({
    desc <- switch(input$model,
                   "SIR"  = "<b>SIR:</b> Susceptible → Infectious → Recovered. Classic model for immunising infections (measles, COVID-19).",
                   "SEIR" = "<b>SEIR:</b> Adds an <i>Exposed</i> (latent) compartment. Individuals are infected but not yet infectious — models incubation periods (flu, Ebola).",
                   "SIRS" = "<b>SIRS:</b> Recovered individuals lose immunity and return to Susceptible. Models endemic diseases with waning immunity (RSV, seasonal flu)."
    )
    div(class="info-box", style="background:#f0f4ff; border-left:4px solid #3f51b5; color:#1a237e;",
        HTML(desc))
  })
  
  # ---- Run simulation (active mode) -------------------------
  sim <- reactive({
    run_sim(input$model, input$time_mode)
  })
  
  run_sim <- function(model, mode) {
    E0 <- if (!is.null(input$E0)) input$E0 else 0
    sigma <- if (!is.null(input$sigma)) input$sigma else 0.2
    omega <- if (!is.null(input$omega)) input$omega else 0.01
    if (mode == "continuous")
      run_continuous(model, input$N, NULL, input$I0, E0,
                     input$beta, input$gamma, sigma, omega, 0, input$nu, input$tmax)
    else
      run_discrete(model, input$N, NULL, input$I0, E0,
                   input$beta, input$gamma, sigma, omega, 0, input$nu, input$tmax)
  }
  
  # ---- Both simulations for comparison ----------------------
  sim_both <- reactive({
    list(
      cont = run_sim(input$model, "continuous"),
      disc = run_sim(input$model, "discrete")
    )
  })
  
  # ---- R₀ and derived metrics --------------------------------
  R0 <- reactive({ input$beta / input$gamma })
  
  output$r0_box <- renderUI({
    r0  <- R0()
    herd <- paste0(round((1 - 1/r0) * 100, 1), "%")
    css  <- if (r0 > 3) "r0-high" else if (r0 > 1) "r0-medium" else "r0-low"
    icon <- if (r0 > 3) "🔴" else if (r0 > 1) "🟡" else "🟢"
    div(class=paste("info-box", css),
        HTML(paste0(icon, " <b>R₀ = β/γ = ", round(r0, 2), "</b><br>",
                    "Herd immunity threshold: <b>", herd, "</b><br>",
                    "Infectious period: <b>", round(1/input$gamma, 1), " days</b>",
                    if (input$model == "SEIR")
                      paste0("<br>Incubation period: <b>",
                             round(1/input$sigma, 1), " days</b>") else "",
                    if (input$model == "SIRS")
                      paste0("<br>Immunity duration: <b>",
                             round(1/input$omega, 1), " days</b>") else ""
        )))
  })
  
  output$stat_cards <- renderUI({
    df <- sim(); req(df)
    peak_I   <- round(max(df$I))
    peak_day <- df$time[which.max(df$I)]
    final_R  <- round(tail(df$R, 1))
    attack   <- paste0(round(final_R / input$N * 100, 1), "%")
    
    div(class="stat-grid",
        div(class="stat-card", "Peak infected", tags$b(comma(peak_I))),
        div(class="stat-card", "Peak day",      tags$b(peak_day)),
        div(class="stat-card", "Final recovered", tags$b(comma(final_R))),
        div(class="stat-card", "Attack rate",   tags$b(attack))
    )
  })
  
  # ---- Compartment colors ------------------------------------
  comp_colors <- c(
    S = "#2196F3", E = "#FF9800", I = "#F44336", R = "#4CAF50"
  )
  comp_labels <- c(S="Susceptible", E="Exposed", I="Infectious", R="Recovered")
  
  # ---- Wide → long helper -----------------------------------
  to_long <- function(df) {
    cols <- intersect(c("S","E","I","R"), names(df))
    do.call(rbind, lapply(cols, function(col) {
      data.frame(time=df$time, value=df[[col]], compartment=col)
    }))
  }
  
  # ---- Mode badge -------------------------------------------
  mode_label <- reactive({
    if (input$time_mode == "continuous") "Continuous" else "Discrete"
  })
  
  base_theme <- function(sz=14) {
    theme_minimal(base_size=sz) +
      theme(panel.grid.minor=element_blank(),
            plot.title=element_text(face="bold", color="#2c3e50"),
            plot.subtitle=element_text(color="#888", size=sz-3),
            legend.position="top")
  }
  
  # ---- Epidemic curves tab ----------------------------------
  output$epi_curve_plot <- renderPlot({
    df <- sim(); req(df)
    dl <- to_long(df)
    cols_used <- intersect(c("S","E","I","R"), unique(dl$compartment))
    dl$compartment <- factor(dl$compartment, levels=c("S","E","I","R"))
    
    g_type <- if (input$time_mode == "discrete") geom_step else geom_line
    
    ggplot(dl, aes(x=time, y=value, color=compartment)) +
      g_type(linewidth=1.2) +
      scale_color_manual(values=comp_colors[cols_used],
                         labels=comp_labels[cols_used]) +
      scale_y_continuous(labels=comma) +
      # Peak infected annotation
      geom_vline(xintercept=df$time[which.max(df$I)],
                 linetype="dashed", color="#F44336", alpha=0.4, linewidth=0.7) +
      labs(title    = paste0(input$model, " Epidemic Curves — ", mode_label()),
           subtitle = paste0("N=", comma(input$N), "  β=", input$beta,
                             "  γ=", input$gamma,
                             "  R₀=", round(R0(), 2),
                             "  |  Dashed = peak infections"),
           x="Time (days)", y="Number of individuals", color=NULL) +
      base_theme()
  })
  
  # ---- Stacked area (compartment dynamics) ------------------
  output$area_plot <- renderPlot({
    df <- sim(); req(df)
    dl <- to_long(df)
    dl$compartment <- factor(dl$compartment, levels=rev(c("S","E","I","R")))
    cols_used <- intersect(c("S","E","I","R"), unique(dl$compartment))
    
    ggplot(dl, aes(x=time, y=value, fill=compartment)) +
      geom_area(alpha=0.85, color="white", linewidth=0.25) +
      scale_fill_manual(values=comp_colors[cols_used],
                        labels=comp_labels[cols_used]) +
      scale_y_continuous(labels=comma) +
      labs(title    = paste0(input$model, " Compartment Dynamics — ", mode_label()),
           subtitle = "Stacked area — total height = N",
           x="Time (days)", y="Number of individuals", fill=NULL) +
      base_theme()
  })
  
  # ---- Log scale tab ----------------------------------------
  output$log_plot <- renderPlot({
    df  <- sim(); req(df)
    dl  <- to_long(df)
    dl  <- dl[dl$value > 0, ]
    cols_used <- intersect(c("S","E","I","R"), unique(dl$compartment))
    dl$compartment <- factor(dl$compartment, levels=c("S","E","I","R"))
    g_type <- if (input$time_mode == "discrete") geom_step else geom_line
    
    ggplot(dl, aes(x=time, y=value, color=compartment)) +
      g_type(linewidth=1.1) +
      scale_color_manual(values=comp_colors[cols_used],
                         labels=comp_labels[cols_used]) +
      scale_y_log10(labels=comma,
                    breaks=10^seq(0, ceiling(log10(input$N)))) +
      labs(title    = paste0(input$model, " — Log Scale"),
           subtitle = "Useful for visualising early exponential growth",
           x="Time (days)", y="Population (log scale)", color=NULL) +
      base_theme()
  })
  
  # ---- Side-by-side comparison ------------------------------
  output$comparison_plot <- renderPlot({
    both <- sim_both()
    req(both$cont, both$disc)
    cols_present <- intersect(c("S","E","I","R"), names(both$cont))
    
    make_plot <- function(df, title_str, is_discrete=FALSE) {
      dl <- to_long(df)
      cols_used <- intersect(c("S","E","I","R"), unique(dl$compartment))
      dl$compartment <- factor(dl$compartment, levels=c("S","E","I","R"))
      g_type <- if (is_discrete) geom_step else geom_line
      ggplot(dl, aes(x=time, y=value, color=compartment)) +
        g_type(linewidth=1.0) +
        scale_color_manual(values=comp_colors[cols_used],
                           labels=comp_labels[cols_used]) +
        scale_y_continuous(labels=comma) +
        labs(title=title_str, x="Time (days)", y="Individuals", color=NULL) +
        base_theme(12) + theme(legend.position="top")
    }
    
    p1 <- make_plot(both$cont, "📈 Continuous (ODE)", FALSE)
    p2 <- make_plot(both$disc, "📈 Discrete (Difference Eq.)", TRUE)
    
    p1 + p2 +
      plot_annotation(
        title    = paste(input$model, "— Continuous vs Discrete"),
        subtitle = paste0("β=", input$beta, "  γ=", input$gamma,
                          "  R₀=", round(R0(), 2), "  N=", comma(input$N)),
        theme = theme(plot.title    = element_text(face="bold", size=15, color="#2c3e50"),
                      plot.subtitle = element_text(color="#888", size=11))
      )
  })
  
  # ---- Theory panel -----------------------------------------
  output$theory_panel <- renderUI({
    HTML("
    <div style='max-width:740px; font-size:14px; line-height:1.85; color:#333; padding:10px'>

    <h4 style='color:#2c3e50'>🦠 Compartmental Models of Infectious Disease</h4>

    <h5>SIR Model</h5>
    <pre style='background:#f4f6f8; padding:12px; border-radius:6px; font-size:13px;'>
  dS/dt = −β·S·I/N
  dI/dt =  β·S·I/N − γ·I
  dR/dt =  γ·I
    </pre>

    <h5>SEIR Model (with latent period)</h5>
    <pre style='background:#fff8e1; padding:12px; border-radius:6px; font-size:13px;'>
  dS/dt = −β·S·I/N
  dE/dt =  β·S·I/N − σ·E
  dI/dt =  σ·E − γ·I
  dR/dt =  γ·I
    </pre>

    <h5>SIRS Model (waning immunity)</h5>
    <pre style='background:#fce4ec; padding:12px; border-radius:6px; font-size:13px;'>
  dS/dt = −β·S·I/N + ω·R
  dI/dt =  β·S·I/N − γ·I
  dR/dt =  γ·I − ω·R
    </pre>

    <h5>Key Parameters</h5>
    <table style='border-collapse:collapse; width:100%; font-size:13px'>
      <tr style='background:#eee'>
        <th style='padding:6px 10px; text-align:left'>Symbol</th>
        <th style='padding:6px 10px; text-align:left'>Meaning</th>
        <th style='padding:6px 10px; text-align:left'>Interpretation</th>
      </tr>
      <tr><td style='padding:5px 10px'>β</td>
          <td>Transmission rate</td>
          <td>Contacts per day × prob. of transmission</td></tr>
      <tr style='background:#f9f9f9'><td style='padding:5px 10px'>γ</td>
          <td>Recovery rate</td>
          <td>1/γ = mean infectious period (days)</td></tr>
      <tr><td style='padding:5px 10px'>σ</td>
          <td>Incubation rate (SEIR)</td>
          <td>1/σ = mean incubation period (days)</td></tr>
      <tr style='background:#f9f9f9'><td style='padding:5px 10px'>ω</td>
          <td>Waning rate (SIRS)</td>
          <td>1/ω = mean duration of immunity (days)</td></tr>
      <tr><td style='padding:5px 10px'>ν</td>
          <td>Vaccination coverage</td>
          <td>Fraction of N immune at t=0</td></tr>
    </table>

    <h5 style='margin-top:18px'>Basic Reproduction Number R₀</h5>
    <p><b>R₀ = β / γ</b> — average number of secondary cases caused by one infectious
    individual in a fully susceptible population.</p>
    <table style='border-collapse:collapse; width:100%; font-size:13px'>
      <tr style='background:#eee'>
        <th style='padding:6px 10px'>R₀</th><th style='padding:6px 10px'>Interpretation</th>
      </tr>
      <tr><td style='padding:5px 10px'>&lt; 1</td>
          <td>Epidemic dies out — each case infects &lt;1 person on average</td></tr>
      <tr style='background:#f9f9f9'><td style='padding:5px 10px'>= 1</td>
          <td>Endemic equilibrium — disease persists at constant level</td></tr>
      <tr><td style='padding:5px 10px'>&gt; 1</td>
          <td>Epidemic grows — higher R₀ = faster spread &amp; higher final size</td></tr>
    </table>

    <h5 style='margin-top:18px'>Herd Immunity Threshold</h5>
    <p>The fraction of the population that must be immune to prevent epidemic growth:<br>
    <b>H = 1 − 1/R₀</b><br>
    Set vaccination coverage ν ≥ H to see the epidemic suppressed.</p>

    <h5>Real-World R₀ Examples</h5>
    <table style='border-collapse:collapse; width:100%; font-size:13px'>
      <tr style='background:#eee'>
        <th style='padding:6px 10px'>Disease</th>
        <th style='padding:6px 10px'>R₀</th>
        <th style='padding:6px 10px'>Herd Immunity Threshold</th>
      </tr>
      <tr><td style='padding:5px 10px'>Seasonal flu</td>  <td>1.2 – 1.4</td><td>~25%</td></tr>
      <tr style='background:#f9f9f9'><td>COVID-19 (original)</td><td>2 – 3</td><td>50–67%</td></tr>
      <tr><td style='padding:5px 10px'>Measles</td>       <td>12 – 18</td>  <td>92–95%</td></tr>
      <tr style='background:#f9f9f9'><td>Smallpox</td>    <td>5 – 7</td>    <td>80–83%</td></tr>
    </table>

    <p style='margin-top:14px'><b>Tip:</b> Use the <b>Side-by-Side</b> tab to compare continuous vs
    discrete dynamics. At high β, discrete models can produce oscillations
    and overshoot impossible in the ODE version.</p>
    </div>
    ")
  })
}

# ================================================================
shinyApp(ui = ui, server = server)
