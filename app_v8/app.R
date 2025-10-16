# ============================================================================
#                                  SHINY APP
# ============================================================================

# Source all modules (loaded once at startup)
source("global.R", local = TRUE)      # Libraries, data, config
source("helpers.R", local = TRUE)     # Data processing functions  
source("modeling.R", local = TRUE)    # Statistical modeling functions

# ============================================================================
#                                USER INTERFACE
# ============================================================================

ui <- fluidPage(
  useShinyjs(),
  
  # Enable Enter key to trigger Run buttons based on active tab
  tags$script(HTML("
    $(document).on('keypress', function(e) {
      if(e.which == 13) {  // Enter key pressed
        var activeTab = $('.nav-tabs .active').text().trim();
        if(activeTab === 'Regression Analysis') {
          $('#run_regression').click();
        } else if(activeTab === 'Mixed Effects Analysis') {
          $('#run_lmer').click();
        } else if(activeTab === 'DiD Analysis') {
          $('#run_did').click();
        }
      }
    });
  ")),
  
  titlePanel("The P-Value Prophet"),
  actionButton("toggleSidebar", "Toggle Sidebar"),
  
  fluidRow(
    column(
      width = 3, id = "sidebar",
      h4("Data Mode:"),
      radioButtons("active_dataset", NULL, 
                   choices = c("Assay Data" = "assay", "Physical Endpoint Data" = "physical"), 
                   selected = "assay"),
      
      # Assay controls
      conditionalPanel(
        condition = "input.active_dataset == 'assay'",
        h4("Processing mode:"),
        radioButtons("mode", NULL,
                     choices = c("Baseline (filter only)" = "baseline",
                                "Well-level CV (wrangle_v3_avg)" = "wrangle", 
                                "Inter-assay CV (inter_assay_CV)" = "inter"),
                     selected = "baseline"),
        br(),
        selectizeInput("assay", "Assay Type:", choices = endpoint_choices_assay, multiple = TRUE),
        selectizeInput("fiber", "Fiber Type:", choices = c("cotton", "pet"), multiple = TRUE),
        selectizeInput("sample_type", "Sample Type:", 
                       choices = if(exists("config") && !is.null(config$samples)) {
                         tolower(config$samples)
                       } else if(exists("final_data") && nrow(final_data) > 0) {
                         sort(unique(tolower(final_data$sample_type)))
                       } else {
                         c("tissue", "gland", "gills", "hemolymph")
                       }, 
                       multiple = TRUE),
        selectizeInput("week", "Week:", choices = week_choices_assay, multiple = TRUE),
        selectizeInput("treatment", "Treatment:", choices = c("Treated", "Untreated", "Control"), multiple = TRUE),
        selectizeInput("fiber_concentration", "Fiber Concentration:", choices = c("0", "100", "1000", "10000"), multiple = TRUE),
        conditionalPanel(condition = "input.mode != 'inter'",
                        selectizeInput("plate_replicate", "Plate Replicate:", choices = NULL, multiple = TRUE),
                        sliderInput("cv_threshold", "CV Threshold (%)", min = 0, max = 100, value = 15)),
        selectizeInput("x_axis_assay", "X-axis:", 
                      choices = c("Fiber Group" = "fiber_group", "Week" = "week", "Sample Type" = "sample_type", 
                                 "Treatment" = "treatment", "Fiber Concentration" = "fiber_concentration"), 
                      selected = "fiber_group"),
        actionButton("run_analysis", "Run Analysis"),
        br(), downloadButton("download_data", "Download Results")
      ),
      
      # Physical controls  
      conditionalPanel(
        condition = "input.active_dataset == 'physical'",
        radioButtons("mode_phys", "Mode:", choices = c("Baseline" = "baseline", "Inter CV" = "inter"), selected = "baseline"),
        selectizeInput("fiber_type_phys", "Fiber Type", choices = c("cotton", "pet"), multiple = TRUE),
        selectizeInput("week_phys", "Week", choices = sort(unique(physical_master$week)), multiple = TRUE),
        selectizeInput("endpoint", "Endpoint", choices = endpoint_choices_physical, multiple = TRUE),
        uiOutput("tissue_type_ui"),
        selectizeInput("fiber_concentration_phys", "Fiber Concentration", choices = c("0","100","1000","10000"), multiple = TRUE),
        selectizeInput("treatment_phys", "Treatment", choices = c("Control","Untreated","Treated"), multiple = TRUE),
        selectizeInput("x_axis_phys", "X Axis", 
                      choices = c("Fiber Group"="fiber_group", "Week"="week", "Fiber Type"="fiber_type", 
                                 "Endpoint"="endpoint", "Tissue Type"="tissue_type"), selected = "week"),
        actionButton("run_phys_analysis", "Run Physical Analysis"),
        downloadButton("download_phys_data", "Download Physical Data")
      )
    ),
    
    column(9,
      tabsetPanel(
        id = "analysis_tabs",
        
        # Visualization Tab
        tabPanel("Visualization",
                h4("Results"), DT::dataTableOutput("table"), br(),
                tags$div(style = "width: 100%; max-width: 1200px; aspect-ratio: 4/3;",
                        plotOutput("dist_plot", width = "100%", height = "100%")), br(), 
                tags$div(style = "width: 100%; max-width: 1200px; aspect-ratio: 4/3;",
                        plotOutput("dot_plot", width = "100%", height = "100%")), br(),
                verbatimTextOutput("info")),
        
        # Regression Analysis Tab 
        tabPanel("Regression Analysis",
                sidebarLayout(
                  sidebarPanel(width = 4,
                              h4("Multiple Linear Regression"),
                              wellPanel(
                                h5("Model Settings"),
                                selectInput("regression_endpoint", "Select Endpoint:", choices = NULL),
                                selectInput("regression_dataset", "Dataset:", 
                                           choices = c("Assay Data" = "assay", "Physical Data" = "physical"), selected = "assay"),
                                checkboxGroupInput("weeks_include", "Include Weeks:",
                                                   choices = c("1" = "1", "3" = "3", "5" = "5"),
                                                   selected = c("1","3","5"), inline = TRUE),
                                checkboxInput("include_three_way", "Include 3-way: fiber × treatment × week", TRUE),
                                checkboxInput("dose_by_fiber", "Include dose × fiber_type", FALSE),
                                checkboxInput("dose_by_treat", "Include dose × chem_treatment", FALSE),
                                actionButton("run_regression", "Run Regression Model", icon = icon("play"))
                              ),
                              wellPanel(
                                h5("Post-hoc Comparisons"),
                                selectInput("emmeans_by", "Calculate EM Means by:",
                                           choices = list("Treatment Groups" = "treatment", "Fiber Type" = "fiber_type",
                                                         "Week" = "week", "Treatment × Week" = "treatment_week",
                                                         "Treatment × Fiber Type" = "treatment_fiber", 
                                                         "Fiber Type × Week" = "fiber_week",
                                                         "Treatment × Fiber × Week" = "all_interactions"),
                                           selected = "treatment"),
                                selectInput("emmeans_dose", "Evaluate at dose:",
                                           choices = c("Controls (0 mf/L)" = "0", "Low dose (100 mf/L)" = "2",
                                                      "Medium dose (1000 mf/L)" = "3", "High dose (10000 mf/L)" = "4"),
                                           selected = "0"),
                                selectInput("comparison_type", "Comparison Type:",
                                           choices = c("All Pairwise" = "pairwise", 
                                                       "Treatment vs Control" = "control"), 
                                           selected = "pairwise"),
                                selectInput("adjustment_method", "P-value Adjustment:",
                                           choices = c("Tukey HSD" = "tukey", "Holm" = "holm", 
                                                      "BH (FDR)" = "fdr", "Bonferroni" = "bonferroni"), selected = "tukey"),
                                actionButton("run_emmeans", "Calculate EM Means", icon = icon("chart-line"))
                              )),
                  mainPanel(width = 8,
                            h5("Model Summary"), 
                            verbatimTextOutput("regression_model_summary"),
                            
                            h5("ANOVA Table"), 
                            DT::dataTableOutput("regression_anova"),
                            
                            h5("Estimated Marginal Means"), 
                            plotOutput("emmeans_plot", height = "420px"),
                            DT::dataTableOutput("emmeans_table"),
                            
                            h5("Pairwise Comparisons"), 
                            DT::dataTableOutput("pairwise_table"),
                            
                            # ============================================================================
                            # TIME TRENDS ANALYSIS - In app.R mainPanel
                            # ============================================================================
                            hr(),
                            h4("Optional: Time Trends Analysis"),
                            
                            wellPanel(
                              style = "background-color: #f9f9f9;",
                              
                              checkboxInput("show_trends", 
                                            "Analyze linear trends over time", 
                                            value = FALSE),
                              
                              helpText(
                                strong("What this does:"), 
                                "Tests whether outcomes change linearly across weeks for each fiber × treatment combination.",
                                br(),
                                "This analysis estimates the slope of change per week and compares slopes between groups.",
                                br(), br(),
                                strong("When to use:"), 
                                "Only useful when 'week' is included in your model and you want to know if trends differ between groups.",
                                style = "color: #555;"
                              ),
                              
                              conditionalPanel(
                                condition = "input.show_trends",
                                
                                hr(),
                                
                                # Filter options for trends
                                h5("Trend Analysis Options"),
                                fluidRow(
                                  column(6,
                                         selectInput("trend_filter_type", "Filter by:",
                                                     choices = c("All Groups" = "none",
                                                                 "Fiber Type" = "fiber",
                                                                 "Treatment" = "treatment",
                                                                 "Both" = "both"),
                                                     selected = "none")
                                  ),
                                  column(6,
                                         conditionalPanel(
                                           condition = "input.trend_filter_type == 'fiber' || input.trend_filter_type == 'both'",
                                           selectInput("trend_fiber_select", "Select Fiber:",
                                                       choices = c("cotton", "pet"),
                                                       selected = "cotton")
                                         ),
                                         conditionalPanel(
                                           condition = "input.trend_filter_type == 'treatment' || input.trend_filter_type == 'both'",
                                           selectInput("trend_treatment_select", "Select Treatment:",
                                                       choices = c("treated", "untreated"),
                                                       selected = "untreated")
                                         )
                                  )
                                ),
                                
                                hr(),
                                
                                h5("Linear Trends Output"),
                                verbatimTextOutput("regression_trends"),
                                
                                br(),
                                
                                h5("Trends Visualization"),
                                plotOutput("regression_trend_plot", height = "500px")
                              )
                            )
                  ))
                
                
              ),
        
        # ---------------------------------------------------------------------------
        # Combined Treatment Analysis tab (UI)
        # ---------------------------------------------------------------------------
        tabPanel(
          "Combined Treatment Analysis",
          icon = icon("layer-group"),
          br(),
          h3("14-Level Combined Treatment Model", style = "color:#2c3e50; font-weight: bold;"),
          
          # Data selection and setup
          wellPanel(
            h5(icon("filter"), "Data Selection", style = "margin-top:0;"),
            fluidRow(
              column(
                4,
                # Endpoint selector
                selectInput(
                  "combined_endpoint",
                  "Select Endpoint",
                  choices = c("Loading..." = "loading")  # updated in server
                )
              ),
              column(
                4,
                # Weeks to include
                checkboxGroupInput(
                  "combined_weeks",
                  "Select Weeks",
                  choices = c("1" = "1", "3" = "3", "5" = "5"),
                  selected = c("1", "3", "5"),
                  inline = TRUE
                )
              ),
              column(
                4,
                # Quick data status
                textOutput("combined_data_info"),
                tags$style("#combined_data_info { font-weight: bold; margin-top: 25px; font-size: 14px; }")
              )
            ),
            hr(),
            h4("Model Setup"),
            fluidRow(
              column(
                4,
                # Reference level for combined factor
                selectInput(
                  "combined_reference",
                  "Reference Level",
                  choices = c(
                    "Control Cotton" = "Control Cotton",
                    "Control PET" = "Control PET",
                    "Untreated Cotton 100 mfL" = "Untreated Cotton 100 mfL",
                    "Treated Cotton 100 mfL"   = "Treated Cotton 100 mfL"
                  ),
                  selected = "Control Cotton"
                )
              ),
              column(
                4,
                checkboxInput("combined_include_week", "Include Week (Time) Variable", value = FALSE)
              ),
              column(
                4,
                # Only show the weeks control once; the primary one above is used
                # Left intentionally blank to keep layout aligned
                div()
              )
            ),
            fluidRow(
              column(
                6,
                actionButton(
                  "run_combined_model",
                  "Run Combined Treatment Model",
                  icon = icon("play"),
                  class = "btn-primary btn-lg",
                  style = "width:100%; margin-top:10px;"
                )
              ),
              column(
                6,
                # Enable download only after a model has been run
                conditionalPanel(
                  condition = "input.run_combined_model > 0",
                  downloadButton(
                    "download_combined_results",
                    "Download Results",
                    class = "btn-success btn-lg",
                    style = "width:100%; margin-top:10px;"
                  )
                )
              )
            )
          ),
          
          # Model Results gated by button click to avoid circular dependency on outputs
          h4("Model Results"),
          conditionalPanel(
            condition = "input.run_combined_model > 0",
            tabsetPanel(
              id = "combined_results_tabs",
              
              tabPanel(
                "Model Summary",
                br(),
                h5("Combined Treatment Variable Levels"),
                verbatimTextOutput("combined_treatment_levels"),
                br(),
                h5("Model Coefficients"),
                verbatimTextOutput("combined_model_summary"),
                br(),
                h5("ANOVA Table"),
                verbatimTextOutput("combined_anova")
              ),
              
              tabPanel(
                "Estimated Means",
                br(),
                h5("Estimated Marginal Means by Treatment Level"),
                DT::dataTableOutput("combined_emm"),
                br(),
                plotOutput("combined_emm_plot", height = "600px")
              ),
              
              tabPanel(
                "Pairwise Comparisons",
                br(),
                fluidRow(
                  column(4,
                         selectInput(
                           "combined_comparison_filter",
                           "Filter Comparisons",
                           choices = c(
                             "All Comparisons" = "all",
                             "vs. Reference Only" = "ref",
                             "Within Fiber vs Control" = "fiber_control",
                             "Within Fiber Type" = "fiber",
                             "Within Concentration" = "conc"
                           ),
                           selected = "fiber_control"
                         )
                  ),
                  column(4,
                         numericInput(
                           "combined_alpha",
                           "Significance Level",
                           value = 0.05,
                           min = 0.001,
                           max = 0.1,
                           step = 0.01
                         )
                  )
                ),
                DT::dataTableOutput("combined_pairwise"),
                br(),
                plotOutput("combined_pairwise_plot", height = "700px")
              ),
              
              tabPanel(
                "Diagnostics",
                br(),
                plotOutput("combined_diagnostics", height = "800px")
              )
            )
          )
        ),
        
        # Mixed Effects Analysis Tab
        tabPanel("Mixed Effects Analysis",
                sidebarLayout(
                  sidebarPanel(width = 4,
                              h4("Linear Mixed Effects Regression"),
                              wellPanel(
                                h5("Model Settings"),
                                checkboxInput("lmer_include_three_way", "Include 3-way: fiber × treatment × week", TRUE),
                                checkboxInput("lmer_dose_by_fiber", "Include dose × fiber_type", FALSE),
                                checkboxInput("lmer_dose_by_treat", "Include dose × chem_treatment", FALSE),
                                checkboxGroupInput("lmer_weeks_include", "Include Weeks:",
                                                   choices = c("1" = "1", "3" = "3", "5" = "5"),
                                                   selected = c("1","3","5"), inline = TRUE),
                                checkboxGroupInput("lmer_random_terms", "Random intercepts:",
                                                  choices = c("Tank" = "tank", "Experiment batch" = "experiment_batch"),
                                                  selected = c("tank","experiment_batch")),
                                actionButton("run_lmer", "Run Mixed Model", icon = icon("play"))
                              ),
                              wellPanel(
                                h5("Post-hoc Comparisons"),
                                selectInput("lmer_emmeans_by", "Calculate EM Means by:",
                                           choices = list("Treatment Groups" = "treatment", "Fiber Type" = "fiber_type", 
                                                         "Week" = "week", "Treatment × Week" = "treatment_week",
                                                         "Treatment × Fiber Type" = "treatment_fiber",
                                                         "Fiber Type × Week" = "fiber_week", 
                                                         "Treatment × Fiber × Week" = "all_interactions"), selected = "treatment"),
                                selectInput("lmer_emmeans_dose", "Evaluate at dose:",
                                           choices = c("Controls (0 mf/L)" = "0", "Low dose (100 mf/L)" = "2",
                                                      "Medium dose (1000 mf/L)" = "3", "High dose (10000 mf/L)" = "4"), selected = "0"),
                                selectInput("lmer_comparison_type", "Comparison Type:",
                                           choices = c("All Pairwise" = "pairwise", 
                                                       "Treatment vs Control" = "control"),
                                           selected = "pairwise"),
                                selectInput("lmer_adjustment_method", "P-value Adjustment:",
                                           choices = c("Tukey HSD" = "tukey", "Holm" = "holm",
                                                      "BH (FDR)" = "fdr", "Bonferroni" = "bonferroni"), selected = "tukey"),
                                actionButton("lmer_run_emmeans", "Calculate EM Means", icon = icon("chart-line"))
                              )),
                  mainPanel(width = 8,
                           h5("Mixed Model Summary"), verbatimTextOutput("lmer_model_summary"),
                           h5("ANOVA Table"), DT::dataTableOutput("lmer_anova"), 
                           h5("Estimated Marginal Means"), plotOutput("lmer_emmeans_plot", height = "420px"),
                           DT::dataTableOutput("lmer_emmeans_table"),
                           h5("Pairwise / Custom Contrasts"), DT::dataTableOutput("lmer_pairwise_table"),
                           h5("Model Comparison (lm vs lmer)"), DT::dataTableOutput("model_compare")))
                ),
        
        # ============================================================================
        # RECOVERY & TRANSLOCATION TAB (standalone, no helpers)
        # ============================================================================
        tabPanel(
          "Recovery & Translocation Analysis",
          fluidRow(
            column(
              4,
              wellPanel(
                h4("Recovery Analysis (Week vs Week Comparison)"),
                
                # Dataset and endpoint
                fluidRow(
                  column(
                    6,
                    selectInput(
                      "recovery_dataset", "Dataset:",
                      choices = c("Assay Data" = "assay", "Physical Data" = "physical"),
                      selected = "assay"
                    )
                  ),
                  column(
                    6,
                    selectInput("recovery_endpoint", "Endpoint:", choices = NULL)
                  )
                ),
                
                # Fixed weeks
                shinyjs::disabled(selectInput("recovery_baseline_week", "Baseline Week:", choices = 5, selected = 5)),
                shinyjs::disabled(selectInput("recovery_recovery_week", "Recovery Week:", choices = 6, selected = 6)),
                
                h5("Recovery-Specific Filters:"),
                
                # Fiber and Treatment
                fluidRow(
                  column(
                    6,
                    checkboxGroupInput(
                      "recovery_fiber_types", "Fiber Types:",
                      choices = c("Cotton" = "cotton", "PET" = "pet"),
                      selected = c("cotton", "pet")
                    )
                  ),
                  column(
                    6,
                    checkboxGroupInput(
                      "recovery_treatments", "Treatments:",
                      choices = c("Untreated" = "untreated", "Treated" = "treated"),
                      selected = c("untreated", "treated")
                    )
                  )
                ),
                
                # Assay: samples + discrete dose
                conditionalPanel(
                  condition = "input.recovery_dataset == 'assay'",
                  checkboxGroupInput(
                    "recovery_samples", "Sample Types:",
                    choices = c("Hemolymph", "Gills", "Gland"),
                    selected = c("Hemolymph", "Gills", "Gland")
                  ),
                  checkboxGroupInput(
                    "recovery_concentration", "Fiber concentration (mf/L):",
                    choices  = c("0" = "0", "100" = "100", "1000" = "1000", "10000" = "10000"),
                    selected = c("0","100","1000","10000"),
                    inline   = TRUE
                  )
                ),
                
                # Physical: tissues + dose for mf_counts
                conditionalPanel(
                  condition = "input.recovery_dataset == 'physical'",
                  checkboxGroupInput("recovery_tissues", "Tissues:",
                                     choices = c("Gills", "Gland", "Tissue"),
                                     selected = c("Gills", "Gland", "Tissue")),
                  checkboxGroupInput(
                    "recovery_concentration_physical", "Fiber concentration (mf/L):",
                    choices  = c("0" = "0", "100" = "100", "1000" = "1000", "10000" = "10000"),
                    selected = c("0","100","1000","10000"),
                    inline   = TRUE
                  )
                ),
                
                # Run controls
                fluidRow(
                  column(6, actionButton("run_recovery", "Run Recovery Analysis", class = "btn-primary")),
                  column(6, checkboxInput("recovery_include_emmeans", "Include EM Means Analysis", value = TRUE))
                ),
                
                helpText("Recovery compares Week 6 vs Week 5 using filters independent of the main sidebar.")
              ),
              
              # Translocation stub
              wellPanel(
                h4("Translocation Analysis"),
                selectInput("translocation_endpoint", "Endpoint:", choices = NULL),
                actionButton("run_translocation", "Run Translocation Analysis")
              )
            ),
            
            # Results
            column(
              8,
              tabsetPanel(
                tabPanel(
                  "Recovery Model",
                  h4("Model Summary"),
                  verbatimTextOutput("recovery_model_summary"),
                  br(),
                  h4("Pairwise Comparisons"),
                  DT::dataTableOutput("recovery_comparisons")
                ),
                tabPanel("Recovery Plot", plotOutput("recovery_plot", height = "500px")),
                tabPanel(
                  "Recovery EM Means",
                  conditionalPanel(
                    condition = "input.recovery_include_emmeans",
                    h4("Estimated Marginal Means"),
                    DT::dataTableOutput("recovery_emmeans_table"),
                    br(),
                    h4("EM Means Plot"),
                    plotOutput("recovery_emmeans_plot", height = "500px")
                  ),
                  conditionalPanel(
                    condition = "!input.recovery_include_emmeans",
                    div(class = "alert alert-info", "Enable 'Include EM Means Analysis' to view emmeans results.")
                  )
                ),
                tabPanel("Translocation Summary", DT::dataTableOutput("translocation_summary")),
                tabPanel("Translocation Plots",
                         plotOutput("translocation_plot", height = "400px"),
                         br(), plotOutput("tissue_comparison_plot", height = "400px"))
              )
            )
          )
        ),
        
        # ============================================================================
        #                                  DID TAB 
        # ============================================================================
        
        tabPanel(
          "Difference-in-Differences",
          fluidRow(
            column(
              4,
              wellPanel(
                h4("DiD Setup"),
                
                # Dataset + endpoint
                fluidRow(
                  column(
                    6,
                    selectInput(
                      "did_dataset", "Dataset:",
                      choices = c("Assay Data" = "assay", "Physical Data" = "physical"),
                      selected = "assay"
                    )
                  ),
                  column(
                    6,
                    selectInput("did_endpoint", "Endpoint:", choices = NULL)
                  )
                ),
                
                conditionalPanel(
                  condition = "input.active_dataset == 'assay'",
                  
                  # DiD Mode Selection
                  radioButtons("did_mode", "DiD Mode:",
                               choices = c(
                                 "Across Fiber (at week)" = "across_fiber",
                                 "Across Week (within fiber)" = "across_week",
                                 "Across Treatment (within week)" = "across_treatment"
                               ),
                               selected = "across_fiber"),
                  
                  # Sample/Tissue selection - only for across_fiber and across_week
                  conditionalPanel(
                    condition = "input.did_mode == 'across_fiber' || input.did_mode == 'across_week'",
                    conditionalPanel(
                      condition = "input.active_dataset == 'assay'",
                      checkboxGroupInput("did_samples", "Sample Types:",
                                         choices = c("Hemolymph", "Gills", "Gland"), 
                                         selected = c("Hemolymph", "Gills", "Gland"),
                                         inline = TRUE)
                    ),
                    conditionalPanel(
                      condition = "input.active_dataset == 'physical'",
                      checkboxGroupInput("did_tissues", "Tissues:",
                                         choices = c("Gills", "Gland", "Tissue"),
                                         selected = c("Gills", "Gland", "Tissue"),
                                         inline = TRUE)
                    )
                  ),
                  
                  # Week selection - only for across_fiber
                  conditionalPanel(
                    condition = "input.did_mode == 'across_fiber'",
                    selectInput("did_week_select", "Select Week:",
                                choices = c("1", "3", "5"), selected = "5")
                  ),
                  
                  # Fiber selection - only for across_week and across_treatment
                  conditionalPanel(
                    condition = "input.did_mode == 'across_week' || input.did_mode == 'across_treatment'",
                    selectInput("did_fiber_select", "Select Fiber:",
                                choices = c("Cotton" = "cotton", "PET" = "pet"), 
                                selected = "cotton")
                  ),
                  
                  # Treatment selection - only for across_treatment
                  conditionalPanel(
                    condition = "input.did_mode == 'across_treatment'",
                    selectInput("did_treatment_select", "Select Treatment:",
                                choices = c("Treated" = "treated", "Untreated" = "untreated"), 
                                selected = "treated")
                  ),
                  
                  # Contrast type - only for across_fiber
                  conditionalPanel(
                    condition = "input.did_mode == 'across_fiber'",
                    radioButtons("did_across_fiber_contrast", "Contrast Type:",
                                 choices = c("DiD (fiber effect on treatment)" = "did",
                                             "Single-treatment fiber contrast" = "single_trt"),
                                 selected = "did")
                  ),
                  
                  # Single treatment selector - only when single_trt is selected
                  conditionalPanel(
                    condition = "input.did_mode == 'across_fiber' && input.did_across_fiber_contrast == 'single_trt'",
                    selectInput("did_single_treatment", "Select Treatment:",
                                choices = c("Treated" = "treated", "Untreated" = "untreated"),
                                selected = "treated")
                  ),
                  
                  # Fiber A and B - only for across_fiber DiD mode
                  conditionalPanel(
                    condition = "input.did_mode == 'across_fiber' && input.did_across_fiber_contrast == 'did'",
                    selectInput("did_fiber_a", "Fiber A:", 
                                choices = c("Cotton" = "cotton", "PET" = "pet"), 
                                selected = "pet"),
                    selectInput("did_fiber_b", "Fiber B:",
                                choices = c("Cotton" = "cotton", "PET" = "pet"), 
                                selected = "cotton")
                  ),
                  
                  # Week A and B - only for across_week
                  conditionalPanel(
                    condition = "input.did_mode == 'across_week'",
                    selectInput("did_week_a", "Week A:", 
                                choices = c("1", "3", "5"), selected = "3"),
                    selectInput("did_week_b", "Week B:", 
                                choices = c("1", "3", "5"), selected = "5")
                  ),
                  
                  # Treatment A and B - only for across_treatment
                  conditionalPanel(
                    condition = "input.did_mode == 'across_treatment'",
                    selectInput("did_trt_a", "Treatment A:",
                                choices = c("Treated" = "treated", "Untreated" = "untreated"),
                                selected = "treated"),
                    selectInput("did_trt_b", "Treatment B:",
                                choices = c("Treated" = "treated", "Untreated" = "untreated"),
                                selected = "untreated")
                  ),
                  
                  # Run button
                  hr(),
                  actionButton("run_did", "Run DiD Analysis", class = "btn-primary")
                )
              )
            ),
            
            column(
              8,
              tabsetPanel(
                tabPanel(
                  "DiD Result",
                  h4("Custom DiD Contrast"),
                  verbatimTextOutput("did_summary"),
                  br(),
                  DT::dataTableOutput("did_table")
                ),
                tabPanel(
                  "EMMeans Grid",
                  DT::dataTableOutput("did_emm_table")
                ),
                tabPanel(
                  "Visualization",
                  plotOutput("did_plot", height = "480px")
                )
              )
            )
          )
        )
      )
    )
  )
)

# ============================================================================
#                             SERVER FUNCTION 
# ============================================================================

server <- function(input, output, session) {
  
  # Toggle sidebar
  observeEvent(input$toggleSidebar, { toggle("sidebar") })
  
  # Update plate replicate choices from pre-loaded data
  observe({
    updateSelectizeInput(session, "plate_replicate", 
                        choices = sort(unique(final_data$plate_replicate)))
  })
  
  # Dynamic tissue filter UI for mf_counts endpoint
  output$tissue_type_ui <- renderUI({
    req(input$endpoint)
    if ("mf_counts" %in% input$endpoint) {
      selectizeInput("tissue_type_phys", "Tissue Type", 
                    choices = unique(na.omit(physical_master$tissue_type)), multiple = TRUE)
    }
  })
  
  # ADD THIS NEW BLOCK HERE - STEP 2: Dynamic sample type choices
  observe({
    # Get available sample types from actual data
    if (exists("final_data") && nrow(final_data) > 0) {
      actual_samples <- sort(unique(tolower(final_data$sample_type)))
      
      # Get config samples (convert to lowercase for consistency)  
      config_samples <- if(exists("config") && !is.null(config$samples)) {
        tolower(config$samples)
      } else {
        c("tissue", "gland", "gills", "hemolymph")
      }
      
      # Use intersection of config and actual data, or just actual data
      available_samples <- if(any(config_samples %in% actual_samples)) {
        intersect(config_samples, actual_samples)
      } else {
        actual_samples
      }
      
      # Update the choices
      updateSelectizeInput(session, "sample_type", 
                           choices = available_samples)
    }
  })
  
  # ---- Helpers for normalization ----
  normalize_fiber <- function(x) {
    x_low <- tolower(x)
    dplyr::case_when(
      x_low %in% c("cotton", "cot") ~ "cotton",
      x_low %in% c("pet", "polyester") ~ "pet",
      TRUE ~ x_low
    )
  }
  
  to_chem_treatment <- function(x) {
    x_low <- tolower(x)
    dplyr::case_when(
      x_low %in% c("untreated", "control") ~ "untreated",
      x_low %in% c("treated") ~ "treated",
      TRUE ~ NA_character_
    )
  }
  
  
  # ============================================================================
  # DATA FILTERING - FIXED: Use reactive() with bindCache(), then bindEvent()
  # ============================================================================
  
  assay_filtered_cached <- reactive({
    df <- final_data  # Pre-processed data
    
    if (length(input$assay)) {
      df <- df %>% filter(assay_type %in% input$assay)
      message("After assay filter: ", nrow(df), " rows")
    }
    
    # FIXED: Case-insensitive fiber type filtering
    if (length(input$fiber)) {
      df <- df %>% filter(tolower(fiber_type) %in% tolower(input$fiber))
      message("After fiber filter: ", nrow(df), " rows")
    }
    
    # FIXED: Case-insensitive sample type filtering
    if (length(input$sample_type)) {
      df <- df %>% filter(tolower(sample_type) %in% tolower(input$sample_type))
      message("After sample_type filter: ", nrow(df), " rows")
    }
    
    if (length(input$week)) {
      df <- df %>% filter(week %in% as.integer(input$week))
      message("After week filter: ", nrow(df), " rows")
    }
    
    if (length(input$fiber_concentration)) {
      df <- df %>% filter(fiber_concentration %in% input$fiber_concentration)
      message("After fiber_concentration filter: ", nrow(df), " rows")
    }
    
    if (length(input$treatment)) {
      if ("Control" %in% input$treatment) {
        df <- df %>% filter(fiber_concentration == "0" | treatment %in% setdiff(input$treatment, "Control"))
      } else {
        df <- df %>% filter(treatment %in% input$treatment)
      }
      message("After treatment filter: ", nrow(df), " rows")
    }
    
    if (input$mode != "inter" && length(input$plate_replicate)) {
      df <- df %>% filter(plate_replicate %in% input$plate_replicate)
      message("After plate_replicate filter: ", nrow(df), " rows")
    }
    
    message("Final filtered data rows: ", nrow(df))
    df
  }) %>% bindCache(input$assay, input$fiber, input$sample_type, input$week, 
                   input$fiber_concentration, input$treatment, input$plate_replicate, input$mode)
  
  assay_filtered <- assay_filtered_cached %>% bindEvent(input$run_analysis)
  
  #Physical data filter - Use reactive with proper caching  
  phys_filtered_cached <- reactive({
    df <- physical_master  # Pre-processed data
    
    if (length(input$fiber_type_phys)) df <- df %>% filter(fiber_type %in% input$fiber_type_phys)
    if (length(input$week_phys)) df <- df %>% filter(week %in% input$week_phys)
    if (length(input$endpoint)) df <- df %>% filter(endpoint %in% input$endpoint)
    if ("mf_counts" %in% input$endpoint && length(input$tissue_type_phys)) {
      df <- df %>% filter(tissue_type %in% input$tissue_type_phys)
    }
    if (length(input$fiber_concentration_phys)) df <- df %>% filter(fiber_concentration %in% input$fiber_concentration_phys)
    if (length(input$treatment_phys)) df <- df %>% filter(treatment %in% input$treatment_phys)
    df
  }) %>% bindCache(input$fiber_type_phys, input$week_phys, input$endpoint, 
                   input$tissue_type_phys, input$fiber_concentration_phys, input$treatment_phys)
  
  # Now use bindEvent
  phys_filtered <- phys_filtered_cached %>% bindEvent(input$run_phys_analysis)
  
  # ============================================================================  
  # DATA SUMMARIZATION - Using helper functions from helpers.R
  # ============================================================================
  
  summarized <- reactive({
    req(input$active_dataset == "assay")
    df <- assay_filtered()
    if (nrow(df) == 0) return(NULL)
    
    threshold <- input$cv_threshold
    if (input$mode == "baseline") {
      summarize_baseline(df, threshold)
    } else if (input$mode == "wrangle") {
      summarize_wrangle(df, threshold)  
    } else {
      summarize_inter(df, threshold)
    }
  })
  
  phys_summarized <- reactive({
    df <- phys_filtered()
    req(nrow(df) > 0)
    
    # Apply blank correction for mf_counts
    df <- df %>%
      group_by(fiber_type, week, tissue_type, endpoint) %>%
      mutate(
        avg_control = if_else(endpoint == "mf_counts", mean(value[treatment == "Control"], na.rm = TRUE), 0),
        value_corr = if_else(endpoint == "mf_counts" & treatment != "Control", value - avg_control, value),
        value_corr = if_else(endpoint == "mf_counts", round(value_corr), value_corr),
        value = if_else(endpoint == "mf_counts", round(value), value)
      ) %>% ungroup()
    
    if (input$mode_phys == "baseline") {
      df
    } else {
      df %>%
        mutate(use_val = if_else(endpoint == "mf_counts", value_corr, value)) %>%
        group_by(fiber_group, week, endpoint, tissue_type, fiber_concentration) %>%
        summarise(
          mean_value = mean(use_val, na.rm = TRUE), sd_value = sd(use_val, na.rm = TRUE),
          n = n(), cv = ifelse(mean_value != 0, 100 * sd_value / mean_value, NA_real_),
          .groups = "drop"
        ) %>%
        mutate(
          mean_value = if_else(abs(mean_value) < 1, round(mean_value, 4), round(mean_value, 2)),
          sd_value = signif(sd_value, 5), cv = round(cv, 2)
        )
    }
  })
  
  # ============================================================================
  # OUTPUTS - Optimized rendering  
  # ============================================================================
  
  active_table <- reactive({
    if (input$active_dataset == "assay") summarized() else phys_summarized()
  })
  
  output$table <- DT::renderDataTable({
    df <- active_table()
    if (is.null(df) || nrow(df) == 0) {
      return(datatable(data.frame(Message = "No data for selected filters / mode")))
    }
    
    if("endpoint" %in% names(df) && "value" %in% names(df)) {
      df <- df %>% mutate(value = if_else(endpoint == "mf_counts", round(value), value))
    }
    
    datatable(df, filter = "top", options = list(pageLength = 15))
  })
  
  # Distribution plot - optimized
  output$dist_plot <- renderPlot({
    df <- active_table()
    req(df, nrow(df) > 0)
    
    if (input$active_dataset == "assay") {
      x_axis_var <- input$x_axis_assay %||% "fiber_group"
      if (x_axis_var == "week") df$week <- factor(df$week)
      
      ggplot(df, aes_string(x = x_axis_var, y = "mean_activity", fill = "fiber_concentration")) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
        facet_wrap(~ assay_type) + theme_minimal() +
        labs(title = glue("Assay Activity by {str_to_title(gsub('_', ' ', x_axis_var))}"),
             x = str_to_title(gsub("_", " ", x_axis_var)), y = "Mean Activity", fill = "Fiber Concentration") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      x_axis_var <- input$x_axis_phys %||% "week"
      x_axis_factor <- paste0("factor(", x_axis_var, ")")
      
      if (input$mode_phys == "baseline") {
        ggplot(df, aes_string(x = x_axis_factor, y = "value", fill = "fiber_concentration")) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.18, alpha = 0.6, size = 2) +
          facet_wrap(~ endpoint, scales = "free_y") + theme_minimal() +
          labs(title = glue("Physical Baseline by {str_to_title(gsub('_', ' ', x_axis_var))}"),
               x = str_to_title(gsub("_", " ", x_axis_var)), y = "Value", fill = "Fiber Concentration") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        ggplot(df, aes_string(x = x_axis_factor, y = "mean_value", fill = "fiber_concentration")) +
          geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
          geom_point(aes_string(group = "fiber_concentration"), position = position_dodge(width = 0.8), 
                    size = 2, color = "black", show.legend = FALSE) +
          facet_wrap(~ endpoint, scales = "free_y") + theme_minimal() +
          labs(title = glue("Physical Inter CV by {str_to_title(gsub('_', ' ', x_axis_var))}"),
               x = str_to_title(gsub("_", " ", x_axis_var)), y = "Mean Value", fill = "Fiber Concentration") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
  })
  
  # Dot plot - assay only
  output$dot_plot <- renderPlot({
    req(input$active_dataset == "assay")
    df <- assay_filtered()
    req(!is.null(df), nrow(df) > 0)
    
    required_cols <- c("fiber_group", "calculated_concentration", "fiber_concentration", "assay_type")
    if (!all(required_cols %in% colnames(df))) return(NULL)
    
    ggplot(df, aes(x = fiber_group, y = calculated_concentration, color = fiber_concentration)) +
      geom_jitter(width = 0.3, alpha = 0.7, size = 2) + facet_wrap(~ assay_type) + theme_minimal() +
      labs(title = "Sample-Level Distribution by Group and Fiber Concentration",
           x = "Fiber Group", y = "Calculated Concentration", color = "Fiber Concentration") +
      theme(plot.title = element_text(size = 14, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Info text
  output$info <- renderText({
    df <- active_table()
    if (is.null(df)) return("No data for current filters / mode.")
    
    if (input$active_dataset == "assay" && "cv_flag" %in% colnames(df)) {
      flagged <- sum(coalesce(df$cv_flag, FALSE), na.rm = TRUE)
      total <- nrow(df)
      glue("Assay mode: {input$mode}. Showing {total} summary rows. {flagged} rows exceed CV threshold ({input$cv_threshold}%).")
    } else {
      glue("Physical mode: {input$mode_phys}. Showing {nrow(df)} rows. Filtered endpoint: {paste(unique(df$endpoint), collapse = ', ')}")
    }
  })
  
  # Download handlers
  output$download_data <- downloadHandler(
    filename = function() { paste0("ECOTOX_", input$mode, "_", Sys.Date(), ".xlsx") },
    content = function(file) {
      df <- summarized()
      writexl::write_xlsx(if (is.null(df)) data.frame(Message = "No data") else df, file)
    }
  )
  
  output$download_phys_data <- downloadHandler(
    filename = function() { paste0("ECOTOX_physical_", Sys.Date(), ".xlsx") },
    content = function(file) {
      df <- phys_summarized()
      writexl::write_xlsx(if (is.null(df)) data.frame(Message = "No data") else df, file)
    }
  )
  
  # ============================================================================
  # REGRESSION ANALYSIS MODULE - Using modeling.R functions
  # ============================================================================
  
  # Dynamic endpoint selection for regression
  observe({
    req(input$regression_dataset)
    if (identical(input$regression_dataset, "assay")) {
      if (exists("final_data") && nrow(final_data) > 0) {
        choices <- sort(unique(final_data$assay_type))
      } else {
        choices <- ""
      }
    } else {
      if (exists("physical_master") && nrow(physical_master) > 0) {
        choices <- sort(unique(physical_master$endpoint))
      } else {
        choices <- ""
      }
    }
    
    if (length(choices) == 0) choices <- ""
    updateSelectInput(session, "regression_endpoint", choices = choices, selected = choices[1])
  })
  
  # Dynamic endpoint selection for recovery
  observe({
    req(input$recovery_dataset)
    if (identical(input$recovery_dataset, "assay")) {
      if (exists("final_data") && nrow(final_data) > 0) {
        choices <- sort(unique(final_data$assay_type))
      } else {
        choices <- ""
      }
    } else {
      if (exists("physical_master") && nrow(physical_master) > 0) {
        choices <- sort(unique(physical_master$endpoint))
      } else {
        choices <- ""
      }
    }
    
    if (length(choices) == 0) choices <- ""
    updateSelectInput(session, "recovery_endpoint", choices = choices, selected = choices[1])
  })
  
  # FIXED: Model-ready data for regression with proper caching
  regression_data_cached <- reactive({
    req(input$regression_dataset, input$regression_endpoint)
    
    if (input$regression_dataset == "assay") {
      df <- final_data %>%
        create_enhanced_treatment_categories() %>%
        mutate(outcome = calculated_concentration) %>%
        filter(assay_type == input$regression_endpoint) %>%
        average_assay_replicates(outcome)
    } else {
      df <- physical_master %>%
        create_enhanced_treatment_categories() %>%
        mutate(outcome = value) %>%
        filter(endpoint == input$regression_endpoint)
    }
    
    wks <- input$weeks_include %||% c("1","3","5")
    df <- df %>% filter(as.character(week) %in% wks)
    
    df %>% normalize_controls_and_dose() %>% droplevels()
  }) %>% bindCache(input$regression_dataset, input$regression_endpoint, input$weeks_include)
  
  regression_data <- regression_data_cached
  
  # FIXED: Mixed effects model data with proper caching
  lmer_data_cached <- reactive({
    req(input$regression_dataset, input$regression_endpoint)
    
    if (input$regression_dataset == "assay") {
      df <- final_data %>%
        create_experiment_batch() %>%
        create_enhanced_treatment_categories() %>%
        mutate(
          outcome = calculated_concentration,
          tank = as.factor(tank),
          week = as.character(week),
          sample_type = as.factor(sample_type)
        ) %>%
        filter(assay_type == input$regression_endpoint) %>%
        average_assay_replicates(outcome) %>%
        mutate(grouping_var = tank)
    } else {
      df <- physical_master %>%
        create_experiment_batch() %>%
        create_enhanced_treatment_categories() %>%
        filter(endpoint == input$regression_endpoint) %>%
        mutate(
          outcome = value,
          week = as.character(week),
          tissue_type = as.factor(tissue_type)
        )

      if (identical(input$regression_endpoint, "mf_counts")) {
        df <- df %>% mutate(
          individual_id = stringr::str_extract(sample, "^[0-9]+"),
          tissue_replicate = sample
        )
      } else {
        df <- df %>% mutate(
          individual_id = as.character(sample)
        )
      }
      df <- df %>% mutate(grouping_var = as.factor(individual_id))
    }

    wks <- input$lmer_weeks_include %||% c("1","3","5")
    df <- df %>%
      filter(as.character(week) %in% wks) %>%
      mutate(
        week = factor(as.character(week), levels = sort(unique(as.character(week))))
      )

    normalize_controls_and_dose(df) %>% droplevels()
  }) %>% bindCache(input$regression_dataset, input$regression_endpoint, input$lmer_weeks_include)
  
  lmer_data <- lmer_data_cached
  
  # Regression model fitting
  regression_model <- eventReactive(input$run_regression, {
    tryCatch({
      df <- regression_data()
      req(nrow(df) > 0)
      
      # Build formula
      form <- outcome ~ fiber_type * chem_treatment + week + 
        fiber_type:week + chem_treatment:week + is_control * fiber_type + dose_log10
      
      if (isTRUE(input$include_three_way)) form <- update(form, . ~ . + fiber_type:chem_treatment:week)
      if (isTRUE(input$dose_by_fiber)) form <- update(form, . ~ . + dose_log10:fiber_type)
      if (isTRUE(input$dose_by_treat)) form <- update(form, . ~ . + dose_log10:chem_treatment)
      
      wks_selected <- input$weeks_include %||% c("1","3","5")
      if (isTRUE(input$include_recovery_interaction) && "6" %in% wks_selected) {
        form <- update(form, . ~ . + is_recovery + is_recovery:chem_treatment + is_recovery:dose_log10)
      }
      
      lm(form, data = df)
    }, error = function(e) { 
      list(error = e$message) 
    })
  })
  
  # Mixed effects model fitting  
  lmer_model <- eventReactive(input$run_lmer, {
    tryCatch({
      df <- lmer_data()
      req(nrow(df) > 0)
      
      message("Building lmer model with ", nrow(df), " observations")
      
      # ENHANCED: More stable model building approach
      
      # Start with basic formula
      fixed_formula <- "outcome ~ fiber_type + chem_treatment + week + dose_log10 + is_control"
      
      # Add interactions carefully to avoid over-parameterization
      if (input$lmer_include_three_way && length(unique(df$week)) > 2) {
        fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment + fiber_type:week + chem_treatment:week")
        
        # Only add 3-way if we have enough data
        if (nrow(df) > 100) {
          fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment:week")
        }
      } else {
        # Add basic 2-way interactions
        fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment")
      }
      
      if (input$lmer_dose_by_fiber) {
        fixed_formula <- paste(fixed_formula, "+ dose_log10:fiber_type")
      }
      
      if (input$lmer_dose_by_treat) {
        fixed_formula <- paste(fixed_formula, "+ dose_log10:chem_treatment")
      }
      
      # SMART: Build random effects based on data availability
      random_terms <- input$lmer_random_terms
      random_formula <- ""
      
      if ("tank" %in% random_terms && "tank" %in% names(df)) {
        # Check if we have enough tanks
        n_tanks <- length(unique(df$tank))
        if (n_tanks > 3) {
          random_formula <- paste(random_formula, "+ (1|tank)")
        } else {
          message("Too few tanks (", n_tanks, ") for random effect, skipping tank")
        }
      }
      
      if ("experiment_batch" %in% random_terms && "experiment_batch" %in% names(df)) {
        # Check if we have enough batches  
        n_batches <- length(unique(df$experiment_batch))
        if (n_batches > 1) {
          random_formula <- paste(random_formula, "+ (1|experiment_batch)")
        } else {
          message("Only one experiment batch, skipping batch random effect")
        }
      }
      
      # Fallback random effect
      if (random_formula == "" && "grouping_var" %in% names(df)) {
        n_groups <- length(unique(df$grouping_var))
        if (n_groups > 3) {
          random_formula <- "+ (1|grouping_var)"
          message("Using grouping_var as fallback random effect")
        }
      }
      
      # Construct final formula
      if (random_formula != "") {
        full_formula <- paste(fixed_formula, random_formula)
      } else {
        # If no random effects possible, fall back to simple lm
        message("No suitable random effects found, this may cause issues")
        full_formula <- fixed_formula
      }
      
      formula_obj <- as.formula(full_formula)
      message("Final formula: ", deparse(formula_obj))
      
      # Fit model with enhanced control parameters for stability
      model <- lmer(formula_obj, data = df, 
                    control = lmerControl(
                      optimizer = "bobyqa",      # More stable optimizer
                      optCtrl = list(maxfun = 20000),  # More iterations
                      check.conv.singular = "warning",  # Don't fail on singular fits
                      calc.derivs = FALSE        # Skip expensive derivative calculations
                    ))
      
      # Check for convergence warnings
      if (length(model@optinfo$conv$lme4$messages) > 0) {
        message("Model convergence warnings: ", paste(model@optinfo$conv$lme4$messages, collapse = "; "))
      }
      
      model
      
    }, error = function(e) {
      message("LMER model fitting error: ", e$message)
      list(error = e$message)
    })
  })
  
  
  # ============================================================================
  # ENHANCED EMMEANS WITH COLORED SIGNIFICANCE & DID ANALYSIS
  # ============================================================================
  
  # EMMEANS calculation for regular regression 
  emmeans_results <- eventReactive(input$run_emmeans, {
    tryCatch({
      # Require that the regression model exists first
      model <- regression_model()
      validate(need(!is.null(model) && !inherits(model, "list"), "Run the regression model first"))
      
      # Calculate dose value based on selection
      dose_val <- switch(input$emmeans_dose,
                         "0" = 0,           # Control level (log10 scale)
                         "2" = log10(100),  # Low dose
                         "3" = log10(1000), # Medium dose  
                         "4" = log10(10000)) # High dose
      
      is_ctrl <- ifelse(dose_val == 0, 1, 0)
      at_list <- list(dose_log10 = dose_val, is_control = is_ctrl)
      
      # Get the grouping specification
      emmeans_spec <- switch(input$emmeans_by,
                             "treatment" = ~ chem_treatment | fiber_type,
                             "fiber_type" = ~ fiber_type,
                             "week" = ~ week,
                             "treatment_week" = ~ chem_treatment | week + fiber_type,
                             "treatment_fiber" = ~ chem_treatment | fiber_type,
                             "fiber_week" = ~ fiber_type | week,
                             "all_interactions" = ~ chem_treatment | fiber_type + week)
      
      # Get the emmeans
      emm <- emmeans(model, specs = emmeans_spec, at = at_list)
      
      list(
        emmeans = emm,
        emmeans_df = as.data.frame(emm)
      )
    }, error = function(e) {
      message("Emmeans error: ", e$message)
      list(error = e$message)
    })
  })
  
  # PAIRWISE COMPARISONS 
  output$pairwise_table <- DT::renderDataTable({
    # Require a fitted regression model
    m <- regression_model()
    validate(need(!is.null(m), "Run the regression model first."))
    
    # Multiple-comparison adjustment (fallback to Tukey)
    adj_method <- if (!is.null(input$adjustment_method) && nzchar(input$adjustment_method)) {
      input$adjustment_method
    } else {
      "tukey"
    }
    
    # Pull emmeans object prepared elsewhere in the regression workflow
    emm_res <- emmeans_results()
    validate(need(!is.null(emm_res$emmeans), "Calculate EM Means first"))
    emm <- emm_res$emmeans
    
    # Contrasts: control (trt.vs.ctrl with untreated as reference) or all pairwise
    contrasts <- if (identical(input$comparison_type, "control")) {
      emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "untreated", adjust = adj_method)
    } else {
      emmeans::contrast(emm, method = "pairwise", adjust = adj_method)
    }
    
    # Coerce to data.frame and ensure p.value exists for consistent formatting
    dfp <- as.data.frame(contrasts)
    if (!"p.value" %in% names(dfp)) dfp$p.value <- NA_real_
    
    # Format numbers and add a significance label
    numeric_cols <- intersect(
      c("estimate", "SE", "df", "t.ratio", "p.value", "lower.CL", "upper.CL"),
      names(dfp)
    )
    
    dfp <- dfp %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(numeric_cols), ~ signif(.x, 5)),
        Significant = dplyr::case_when(
          is.na(p.value)        ~ "Unable to compute",
          p.value < 0.001       ~ "Yes (***)",
          p.value < 0.01        ~ "Yes (**) ",
          p.value < 0.05        ~ "Yes (*)",
          TRUE                  ~ "No"
        ),
        Adjustment = adj_method
      )
    
    # Render styled table
    DT::datatable(
      dfp,
      options = list(pageLength = 20, scrollX = TRUE, dom = "Blfrtip"),
      filter = "top",
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        "Significant",
        backgroundColor = DT::styleEqual(
          c("Yes (***)", "Yes (**) ", "Yes (*)", "No", "Unable to compute"),
          c("#28a745",  "#28a745",    "#28a745", "#ffffff", "#ffe6cc")
        ),
        color = DT::styleEqual(
          c("Yes (***)", "Yes (**) ", "Yes (*)", "No", "Unable to compute"),
          c("#ffffff",  "#ffffff",    "#ffffff", "#000000", "#000000")
        )
      ) %>%
      { if (all(c("estimate","SE","t.ratio") %in% names(dfp)))
        DT::formatRound(., columns = c("estimate","SE","t.ratio"), digits = 4) else . } %>%
      { if ("p.value" %in% names(dfp))
        DT::formatRound(., columns = "p.value", digits = 6) else . }
  })
  
  # ENHANCED EMMEANS TABLE WITH SIGNIFICANCE COLORS  
  output$emmeans_table <- DT::renderDataTable({
    results <- emmeans_results()
    validate(need(!is.null(results$emmeans_df), "Click 'Calculate EM Means' to generate results"))
    
    df_emm <- results$emmeans_df %>%
      mutate(across(where(is.numeric), ~ signif(.x, 4)))
    
    DT::datatable(
      df_emm, 
      options = list(pageLength = 15, scrollX = TRUE, dom = 'Blfrtip'), 
      filter = "top", 
      rownames = FALSE
    ) %>%
      DT::formatRound(columns = c("emmean", "SE", "lower.CL", "upper.CL"), digits = 4)
  })
  
  # ENHANCED EMMEANS PLOT
  output$emmeans_plot <- renderPlot({
    results <- emmeans_results()
    validate(need(!is.null(results$emmeans), "Calculate EM Means first"))
    
    tryCatch({
      # USE THE SAME APPROACH AS YOUR OLD APP - built-in emmeans plot
      plot(results$emmeans) + 
        theme_minimal() + 
        labs(
          title = "Estimated Marginal Means with 95% Confidence Intervals",
          subtitle = paste("Dose level:", input$emmeans_dose, "| Adjustment:", input$adjustment_method),
          x = "Estimated Marginal Mean"
        ) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          panel.grid.minor = element_blank()
        )
    }, error = function(e) {
      # Fallback if emmeans plot fails - use FIXED geom_errorbar (not geom_errorbarh)
      df <- results$emmeans_df
      ggplot(df, aes(x = emmean, y = rownames(df))) +
        geom_point(size = 3, color = "#2c3e50") +
        # FIXED: Use geom_errorbar with orientation instead of deprecated geom_errorbarh
        geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), 
                      width = 0.2, orientation = "y") +
        theme_minimal() +
        labs(
          title = "Estimated Marginal Means", 
          x = "Estimated Mean", 
          y = "Groups"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
  })
  
  # ============================================================================
  # LINEAR TRENDS ANALYSIS - SERVER CODE (OPTION 2: REFIT WITH NUMERIC WEEK)
  # ============================================================================
  
  # Helper: Get filtered trend specification
  get_filtered_trend_spec <- reactive({
    req(regression_model())
    
    filter_type <- input$trend_filter_type
    
    if (filter_type == "none") {
      return(~ fiber_type + chem_treatment)
    } else if (filter_type == "fiber") {
      return(list(~ chem_treatment, fiber_type = input$trend_fiber_select))
    } else if (filter_type == "treatment") {
      return(list(~ fiber_type, chem_treatment = input$trend_treatment_select))
    } else if (filter_type == "both") {
      return(list(fiber_type = input$trend_fiber_select, 
                  chem_treatment = input$trend_treatment_select))
    }
  })
  
  # ===========================
  # LINEAR TRENDS TEXT OUTPUT
  # ===========================
  output$regression_trends <- renderPrint({
    req(input$show_trends == TRUE)
    req(regression_model())
    
    tryCatch({
      mdl <- regression_model()
      
      # Check if week is in the model
      model_terms <- attr(terms(mdl), "term.labels")
      has_week <- any(grepl("week", model_terms, ignore.case = TRUE))
      
      if (!has_week) {
        cat("Week is not included in the current model.\n")
        cat("To use time trends analysis, select week checkboxes in 'Include Weeks'.\n")
        return(invisible())
      }
      
      # Get the original data and formula from the model
      model_data <- mdl$model
      original_formula <- formula(mdl)
      
      # Check if week variable exists
      if (!"week" %in% names(model_data)) {
        cat("Week variable not found in model data.\n")
        return(invisible())
      }
      
      # Create a CLEAN dataset with numeric week
      trend_data <- model_data
      
      if (is.factor(trend_data$week)) {
        trend_data$week <- as.numeric(as.character(trend_data$week))
        cat("Note: Week converted from factor to numeric for trend analysis\n")
        cat("Week values: ", paste(sort(unique(trend_data$week)), collapse = ", "), "\n\n")
      }
      
      cat("=== LINEAR TRENDS OVER TIME ===\n\n")
      
      # Display what's being analyzed
      if (input$trend_filter_type == "none") {
        cat("Analyzing trends for: ALL fiber types and treatments\n\n")
      } else if (input$trend_filter_type == "fiber") {
        cat("Analyzing trends for: ", toupper(input$trend_fiber_select), 
            " fiber (both treatments)\n\n")
      } else if (input$trend_filter_type == "treatment") {
        cat("Analyzing trends for: ", toupper(input$trend_treatment_select), 
            " (both fiber types)\n\n")
      } else {
        cat("Analyzing trends for: ", toupper(input$trend_fiber_select), 
            " fiber, ", toupper(input$trend_treatment_select), " treatment\n\n")
      }
      
      # Fit a NEW model from scratch with numeric week
      cat("Refitting model with numeric week...\n")
      mdl_numeric <- lm(original_formula, data = trend_data)
      cat("Model refit successful.\n\n")
      
      # Compute linear trends with the numeric week model
      if (input$trend_filter_type == "none") {
        trends <- emmeans::emtrends(mdl_numeric, ~ fiber_type + chem_treatment, 
                                    var = "week")
      } else if (input$trend_filter_type == "fiber") {
        trends <- emmeans::emtrends(mdl_numeric, ~ chem_treatment, 
                                    var = "week",
                                    at = list(fiber_type = input$trend_fiber_select))
      } else if (input$trend_filter_type == "treatment") {
        trends <- emmeans::emtrends(mdl_numeric, ~ fiber_type, 
                                    var = "week",
                                    at = list(chem_treatment = input$trend_treatment_select))
      } else {
        trends <- emmeans::emtrends(mdl_numeric, ~ 1, 
                                    var = "week",
                                    at = list(fiber_type = input$trend_fiber_select,
                                              chem_treatment = input$trend_treatment_select))
      }
      
      cat("--- Estimated Linear Trends (Slopes) ---\n")
      cat("How much does the outcome change per 1-week increase?\n\n")
      print(summary(trends, infer = TRUE))
      
      if (input$trend_filter_type == "none" || 
          (input$trend_filter_type != "both")) {
        cat("\n\n--- Pairwise Comparisons of Trends ---\n")
        cat("Are slopes significantly different between groups?\n\n")
        trend_pairs <- pairs(trends)
        print(summary(trend_pairs))
      }
      
      cat("\n\n=== INTERPRETATION GUIDE ===\n")
      cat("- week.trend: Change in outcome per 1-week increase\n")
      cat("- SE: Standard error of the trend estimate\n")
      cat("- t.ratio & p.value: Test if trend differs from zero\n")
      cat("- lower.CL & upper.CL: 95% confidence interval\n")
      cat("\nPositive trend = outcome increases over time\n")
      cat("Negative trend = outcome decreases over time\n")
      cat("\nNote: Trends are averaged across all dose levels and concentrations.\n")
      
    }, error = function(e) {
      cat("Error computing trends:\n")
      cat(conditionMessage(e), "\n\n")
      cat("This may occur if:\n")
      cat("- Model cannot be refit with numeric week\n")
      cat("- Insufficient data for selected filter\n")
      cat("- Model lacks necessary interaction terms\n")
    })
  })
  
  # ===========================
  # TREND VISUALIZATION PLOT
  # ===========================
  output$regression_trend_plot <- renderPlot({
    req(input$show_trends == TRUE)
    req(regression_model())
    
    tryCatch({
      mdl <- regression_model()
      
      # Check if week is in the model
      model_terms <- attr(terms(mdl), "term.labels")
      has_week <- any(grepl("week", model_terms, ignore.case = TRUE))
      
      if (!has_week) {
        plot.new()
        text(0.5, 0.5, "Week not included in model\nSelect week checkboxes in 'Include Weeks'", 
             cex = 1.3, col = "#d9534f", font = 2)
        return(invisible())
      }
      
      # Get estimated marginal means for plotting
      # (Uses original model - plotting doesn't need numeric week)
      if (input$trend_filter_type == "none") {
        emm_week <- emmeans::emmeans(mdl, ~ week | fiber_type + chem_treatment)
      } else if (input$trend_filter_type == "fiber") {
        emm_week <- emmeans::emmeans(mdl, ~ week | chem_treatment,
                                     at = list(fiber_type = input$trend_fiber_select))
      } else if (input$trend_filter_type == "treatment") {
        emm_week <- emmeans::emmeans(mdl, ~ week | fiber_type,
                                     at = list(chem_treatment = input$trend_treatment_select))
      } else {
        emm_week <- emmeans::emmeans(mdl, ~ week,
                                     at = list(fiber_type = input$trend_fiber_select,
                                               chem_treatment = input$trend_treatment_select))
      }
      
      emm_df <- as.data.frame(emm_week)
      
      # Convert week to numeric for plotting
      if (is.factor(emm_df$week)) {
        emm_df$week_num <- as.numeric(as.character(emm_df$week))
      } else {
        emm_df$week_num <- emm_df$week
      }
      
      # Create plot based on filter selection
      library(ggplot2)
      
      if (input$trend_filter_type == "none") {
        # Full plot with facets
        p <- ggplot(emm_df, aes(x = week_num, y = emmean, 
                                color = chem_treatment, 
                                linetype = chem_treatment,
                                shape = chem_treatment)) +
          geom_line(size = 1.2) +
          geom_point(size = 3.5) +
          geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                        width = 0.2, size = 0.8, alpha = 0.7) +
          facet_wrap(~ fiber_type, nrow = 1, scales = "free_y") +
          theme_bw(base_size = 14) +
          theme(
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 12),
            legend.text = element_text(size = 11),
            strip.background = element_rect(fill = "#5bc0de"),
            strip.text = element_text(face = "bold", size = 13, color = "white"),
            panel.grid.minor = element_blank()
          ) +
          labs(
            title = "Estimated Marginal Means: Trends Over Time",
            subtitle = "Linear trends show rate of change per week. Error bars = 95% CI",
            x = "Week", 
            y = "Estimated Mean Outcome",
            color = "Treatment", 
            linetype = "Treatment",
            shape = "Treatment"
          )
        
      } else if (input$trend_filter_type == "fiber") {
        # Filter by fiber only
        p <- ggplot(emm_df, aes(x = week_num, y = emmean, 
                                color = chem_treatment, 
                                linetype = chem_treatment,
                                shape = chem_treatment)) +
          geom_line(size = 1.3) +
          geom_point(size = 4) +
          geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                        width = 0.2, size = 0.9, alpha = 0.7) +
          theme_bw(base_size = 15) +
          theme(
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 13),
            legend.text = element_text(size = 12),
            panel.grid.minor = element_blank()
          ) +
          labs(
            title = paste0("Trends Over Time: ", toupper(input$trend_fiber_select), " Fiber"),
            subtitle = "Linear trends show rate of change per week. Error bars = 95% CI",
            x = "Week", 
            y = "Estimated Mean Outcome",
            color = "Treatment", 
            linetype = "Treatment",
            shape = "Treatment"
          )
        
      } else if (input$trend_filter_type == "treatment") {
        # Filter by treatment only
        p <- ggplot(emm_df, aes(x = week_num, y = emmean, 
                                color = fiber_type, 
                                linetype = fiber_type,
                                shape = fiber_type)) +
          geom_line(size = 1.3) +
          geom_point(size = 4) +
          geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                        width = 0.2, size = 0.9, alpha = 0.7) +
          theme_bw(base_size = 15) +
          theme(
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = 13),
            legend.text = element_text(size = 12),
            panel.grid.minor = element_blank()
          ) +
          labs(
            title = paste0("Trends Over Time: ", toupper(input$trend_treatment_select), " Treatment"),
            subtitle = "Linear trends show rate of change per week. Error bars = 95% CI",
            x = "Week", 
            y = "Estimated Mean Outcome",
            color = "Fiber Type", 
            linetype = "Fiber Type",
            shape = "Fiber Type"
          )
        
      } else {
        # Both filters applied - single line
        p <- ggplot(emm_df, aes(x = week_num, y = emmean)) +
          geom_line(size = 1.5, color = "#0073C2") +
          geom_point(size = 5, color = "#0073C2") +
          geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                        width = 0.2, size = 1, color = "#0073C2", alpha = 0.7) +
          theme_bw(base_size = 16) +
          theme(
            panel.grid.minor = element_blank()
          ) +
          labs(
            title = paste0("Trend: ", toupper(input$trend_fiber_select), " + ", 
                           toupper(input$trend_treatment_select)),
            subtitle = "Linear trend shows rate of change per week. Error bars = 95% CI",
            x = "Week", 
            y = "Estimated Mean Outcome"
          )
      }
      
      # Add custom x-axis breaks at actual week values
      p <- p + scale_x_continuous(breaks = unique(emm_df$week_num))
      print(p)
      
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error creating plot:\n", conditionMessage(e)), 
           cex = 1.1, col = "#d9534f")
    })
  })
  
  # =============================================================================
  # Combined Treatment Analysis (server)
  # - Uses snake_case IDs
  # - Harmonizes legacy/snake_case columns
  # - Avoids circular gating by rendering even when hidden and gating by button
  # =============================================================================
  
  # Small helper to format a quick data status line in the UI
  output$combined_data_info <- renderText({
    if (is.null(input$combined_endpoint) || identical(input$combined_endpoint, "loading")) {
      return("No data loaded")  # initial state
    }
    # Try to count rows safely
    n <- tryCatch({
      df <- combined_base_data()
      if (is.null(df)) 0L else nrow(df)
    }, error = function(e) 0L)
    if (n > 0L) paste0(format(n, big.mark = ","), " observations loaded") else "No data for selected endpoint"
  })
  
  # -----------------------------------------------------------------------------
  # REACTIVE: Base data for combined analysis
  # -----------------------------------------------------------------------------
  combined_base_data <- reactive({
    req(input$combined_endpoint)
    req(input$active_dataset)  # ensure your sidebar id uses snake_case
    
    cat("\n=== combined_base_data reactive triggered ===\n")
    cat("Dataset:", input$active_dataset, "\n")
    cat("Endpoint:", input$combined_endpoint, "\n")
    
    if (identical(input$active_dataset, "assay")) {
      req(final_data)
      df <- final_data %>%
        create_enhanced_treatment_categories() %>%        # adds fiber/treatment normalization
        dplyr::mutate(outcome = calculated_concentration)
      
      # CRITICAL: Convert outcome to numeric and drop invalid rows BEFORE averaging
      df <- df %>%
        dplyr::mutate(
          outcome = suppressWarnings(as.numeric(outcome))
        ) %>%
        dplyr::filter(!is.na(outcome))
      
      # Now filter and average replicates
      df <- df %>%
        dplyr::filter(assay_type == input$combined_endpoint) %>%
        average_assay_replicates(outcome)               # NSE, no quotes
      
      cat("Assay data filtered. Rows:", nrow(df), "\n")
      
    } else if (identical(input$active_dataset, "physical")) {
      req(physical_master)
      df <- physical_master %>%
        create_enhanced_treatment_categories() %>%
        dplyr::mutate(outcome = value) %>%
        dplyr::filter(endpoint == input$combined_endpoint)
      cat("Physical data filtered. Rows:", nrow(df), "\n")
      
    } else {
      stop("Invalid dataset selection")
    }
    
    # Weeks filter and normalization
    wks <- input$combined_weeks
    if (is.null(wks) || length(wks) == 0) {
      wks <- c("1", "3", "5")
      cat("No weeks selected, using default: 1, 3, 5\n")
    }
    
    df <- df %>%
      dplyr::filter(as.character(week) %in% wks) %>%
      normalize_controls_and_dose() %>%
      droplevels()
    
    cat("After week filter and normalization. Final rows:", nrow(df), "\n")
    cat("Columns:", paste(names(df), collapse = ", "), "\n\n")
    df
  })
  
  # =============================================================================
  # Combined Treatment Analysis: Harmonize and build 14-level factor
  # =============================================================================
  combined_data <- reactive({
    df <- combined_base_data()
    req(df, nrow(df) > 0)
    
    # 1) Harmonize column names
    if ("fibertype" %in% names(df) && !"fiber_type" %in% names(df)) {
      df <- dplyr::rename(df, fiber_type = fibertype)
    }
    if ("chemtreatment" %in% names(df) && !"chem_treatment" %in% names(df)) {
      df <- dplyr::rename(df, chem_treatment = chemtreatment)
    }
    if ("iscontrol" %in% names(df) && !"is_control" %in% names(df)) {
      df <- dplyr::rename(df, is_control = iscontrol)
    }
    
    # 2) Build fiber_concentration_numeric ONLY from columns that exist
    if (!"fiber_concentration_numeric" %in% names(df)) {
      # Start with a fallback column of zeros
      fcn <- rep(0, nrow(df))
      
      # Overwrite with first available numeric source
      if ("fiber_concentration" %in% names(df)) {
        tmp <- suppressWarnings(as.numeric(df[["fiber_concentration"]]))
        if (length(tmp) == nrow(df)) fcn <- dplyr::coalesce(tmp, fcn)
      }
      # Do NOT reference dose_log10/doselog10 unless they exist
      # (combined_base_data likely doesn't have them)
      
      df[["fiber_concentration_numeric"]] <- fcn
    }
    
    # 3) Call the wrapper to create the combined factor
    df_combined <- create_combined_treatment(
      data = df,
      reference_level = input$combined_reference
    )
    
    # 4) Normalize final names
    if ("treatmentcombined" %in% names(df_combined) && !"treatment_combined" %in% names(df_combined)) {
      df_combined <- dplyr::rename(df_combined, treatment_combined = treatmentcombined)
    }
    
    # 5) *** FIXED: TRIM WHITESPACE FROM treatment_combined FACTOR ***
    # This fixes the "Control Cotton " (with trailing space) issue
    if ("treatment_combined" %in% names(df_combined)) {
      df_combined <- df_combined %>%
        dplyr::mutate(
          # Convert to character, trim whitespace, then back to factor
          treatment_combined = stringr::str_trim(as.character(treatment_combined)),
          # Re-factorize with cleaned levels
          treatment_combined = factor(
            treatment_combined,
            levels = unique(stringr::str_trim(as.character(treatment_combined)))
          )
        )
      
      # Remove "Other" or NA levels if present
      df_combined <- df_combined %>%
        dplyr::filter(!is.na(treatment_combined), treatment_combined != "Other")
      
      # Drop unused factor levels after filtering
      df_combined$treatment_combined <- droplevels(df_combined$treatment_combined)
    }
    
    # 6) Final validation
    req("outcome" %in% names(df_combined), "treatment_combined" %in% names(df_combined))
    req(nrow(df_combined) > 0)  # Also check we have data left after filtering
    
    df_combined
  })
  
  # =============================================================================
  #         OBSERVER TO DYNAMICALLY UPDATE THE REFERENCE LEVEL CHOICES
  # =============================================================================
  
  observe({
    # Trigger when combined_base_data changes (endpoint/week selection)
    req(input$combined_endpoint, length(input$combined_weeks) > 0)
    
    tryCatch({
      # Get the actual data
      data <- combined_data()
      
      if (!is.null(data) && "treatment_combined" %in% names(data)) {
        # Get actual levels from the data
        all_levels <- levels(data$treatment_combined)
        
        # Filter to only control groups (levels that start with "Control")
        control_levels <- all_levels[stringr::str_detect(all_levels, "^Control")]
        
        # If no control levels found, fall back to first level
        if (length(control_levels) == 0) {
          control_levels <- all_levels[1]
        }
        
        # Update the selectInput with only control levels
        updateSelectInput(
          session,
          "combined_reference",
          choices = control_levels,
          selected = control_levels[1]  # Default to first control
        )
        
        cat("✓ Reference level options updated:", paste(control_levels, collapse = ", "), "\n")
      }
    }, error = function(e) {
      # Silently fail if data not ready
      NULL
    })
  }) %>%
    bindEvent(input$combined_endpoint, input$combined_weeks)
  
  # --- Fit the model -----------------------------------------------------------
  combined_model <- eventReactive(input$run_combined_model, {
    data <- combined_data()
    req(data, nrow(data) > 0)
    
    # Get reference level from input
    ref_level <- input$combined_reference
    
    # Ensure reference level exists in data
    if (!(ref_level %in% unique(data$treatment_combined))) {
      showNotification(
        paste("Reference level", ref_level, "not found in data."),
        type = "warning"
      )
      ref_level <- levels(data$treatment_combined)[1]
    }
    
    # Relevel the factor with the selected reference
    data$treatment_combined <- stats::relevel(
      data$treatment_combined,
      ref = ref_level
    )
    
    # Check that we have enough levels
    n_lvls <- nlevels(droplevels(data$treatment_combined))
    validate(
      need(n_lvls >= 2, 
           sprintf("Not enough treatment levels after filtering (have %d). Try different weeks or endpoint.", n_lvls))
    )
    
    # Build formula
    f <- if (isTRUE(input$combined_include_week)) {
      stats::as.formula("outcome ~ treatment_combined + week")
    } else {
      stats::as.formula("outcome ~ treatment_combined")
    }
    
    # Fit model
    stats::lm(f, data = data)
  })
  
  # --- Outputs: keep ids consistent with UI and use validate() -----------------
  output$combined_treatment_levels <- renderPrint({
    data <- combined_data()
    validate(need(nrow(data) > 0, "No rows available for the selected filters."))
    level_summary <- data %>%
      dplyr::group_by(treatment_combined) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean_value = mean(outcome, na.rm = TRUE),
        sd_value   = stats::sd(outcome, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(treatment_combined)
    print(level_summary, n = Inf)
  })
  
  output$combined_model_summary <- renderPrint({
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Click Run to fit the model, or check data/endpoint selection."))
    print(summary(mdl))
  })
  
  output$combined_anova <- renderPrint({
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Model not available."))
    print(stats::anova(mdl))
  })
  
  # ============================================================================
  # COMBINED MODEL: EMMEANS CALCULATION (Reactive)
  # ============================================================================
  
  combined_emmeans_data <- eventReactive(input$run_combined_model, {
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Model not available"))
    
    tryCatch({
      # Calculate emmeans for treatment_combined factor
      emm <- emmeans::emmeans(mdl, specs = ~ treatment_combined)
      
      # Convert to data frame
      emm_df <- as.data.frame(emm) %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~signif(.x, 4)))
      
      list(
        emmeans = emm,
        emmeans_df = emm_df
      )
    }, error = function(e) {
      message("Emmeans calculation error: ", e$message)
      list(error = e$message)
    })
  })
  
  # ============================================================================
  # COMBINED MODEL: PAIRWISE COMPARISONS (Reactive)
  # ============================================================================
  
  combined_pairwise_data <- eventReactive(input$run_combined_model, {
    emm_res <- combined_emmeans_data()
    validate(need(!is.null(emm_res$emmeans), "Calculate emmeans first"))
    
    tryCatch({
      emm <- emm_res$emmeans
      ref_level <- input$combined_reference
      filter_type <- input$combined_comparison_filter
      
      # Calculate all pairwise comparisons
      pw <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
      pw_df <- as.data.frame(pw)
      
      # Add significance indicator
      pw_df <- pw_df %>%
        dplyr::mutate(
          significant = dplyr::case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ ""
          )
        )
      
      # Apply filtering based on user selection
      if (filter_type == "ref") {
        # Only comparisons vs reference
        pw_df <- pw_df %>%
          dplyr::filter(stringr::str_detect(contrast, stringr::fixed(ref_level)))
      } else if (filter_type == "fiber") {
        # Within fiber type: Cotton vs Cotton, PET vs PET
        pw_df <- pw_df %>%
          dplyr::filter(
            (stringr::str_detect(contrast, "Cotton") & 
               stringr::str_count(contrast, "Cotton") == 2) |
              (stringr::str_detect(contrast, "PET|Pet") & 
                 stringr::str_count(contrast, "PET|Pet") >= 2)
          )
      } else if (filter_type == "conc") {
        # Within concentration only
        pw_df <- pw_df %>%
          dplyr::mutate(
            # Extract concentration from each side of contrast
            conc1 = stringr::str_extract(contrast, "\\d+(?= mfL)"),
            conc2 = stringr::str_extract(
              stringr::str_remove(contrast, "^[^-]+ - "), 
              "\\d+(?= mfL)"
            )
          ) %>%
          dplyr::filter(conc1 == conc2 | (is.na(conc1) & is.na(conc2))) %>%
          dplyr::select(-conc1, -conc2)
      }
      # else filter_type == "all", keep everything
      
      # Format numeric columns
      pw_df <- pw_df %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~signif(.x, 4)))
      
      list(
        pairwise = pw,
        pairwise_df = pw_df
      )
    }, error = function(e) {
      message("Pairwise comparison error: ", e$message)
      list(error = e$message)
    })
  })
  
  # ============================================================================
  #         CUSTOM CONTRAST OPTION FOR WITHIN-FIBER COMPARISONS
  # ============================================================================
  
  combined_fiber_contrasts <- eventReactive(input$run_combined_model, {
    emm_res <- combined_emmeans_data()
    validate(need(!is.null(emm_res$emmeans), "Calculate emmeans first"))
    
    tryCatch({
      emm <- emm_res$emmeans
      emm_df <- as.data.frame(emm)
      
      # Separate Cotton and PET groups
      cotton_levels <- emm_df %>%
        dplyr::filter(stringr::str_detect(treatment_combined, "Cotton")) %>%
        dplyr::pull(treatment_combined)
      
      pet_levels <- emm_df %>%
        dplyr::filter(stringr::str_detect(treatment_combined, "PET|Pet")) %>%
        dplyr::pull(treatment_combined)
      
      # Build contrast lists
      contrast_list <- list()
      
      # Cotton vs Cotton Control
      cotton_control <- cotton_levels[stringr::str_detect(cotton_levels, "Control")]
      if (length(cotton_control) > 0 && length(cotton_levels) > 1) {
        other_cotton <- setdiff(cotton_levels, cotton_control)
        for (level in other_cotton) {
          contrast_name <- paste(level, "-", cotton_control)
          contrast_list[[contrast_name]] <- emmeans::contrast(
            emm,
            method = list(setNames(
              c(-1, 1)[match(c(cotton_control, level), levels(emm@levels$treatment_combined))],
              c(cotton_control, level)
            ))
          )
        }
      }
      
      # PET vs PET Control
      pet_control <- pet_levels[stringr::str_detect(pet_levels, "Control")]
      if (length(pet_control) > 0 && length(pet_levels) > 1) {
        other_pet <- setdiff(pet_levels, pet_control)
        for (level in other_pet) {
          contrast_name <- paste(level, "-", pet_control)
          contrast_list[[contrast_name]] <- emmeans::contrast(
            emm,
            method = list(setNames(
              c(-1, 1)[match(c(pet_control, level), levels(emm@levels$treatment_combined))],
              c(pet_control, level)
            ))
          )
        }
      }
      
      # Combine into single data frame if we have contrasts
      if (length(contrast_list) > 0) {
        # Use Tukey-adjusted pairwise within each fiber type
        cotton_contrast <- if (length(other_cotton) > 0) {
          emmeans::contrast(
            emm,
            method = "trt.vs.ctrl",
            ref = which(levels(emm@levels$treatment_combined) == cotton_control)
          ) %>% as.data.frame()
        } else NULL
        
        pet_contrast <- if (length(other_pet) > 0) {
          emmeans::contrast(
            emm,
            method = "trt.vs.ctrl",
            ref = which(levels(emm@levels$treatment_combined) == pet_control)
          ) %>% as.data.frame()
        } else NULL
        
        contrast_df <- dplyr::bind_rows(
          cotton_contrast,
          pet_contrast
        ) %>%
          dplyr::mutate(
            significant = dplyr::case_when(
              p.value < 0.001 ~ "***",
              p.value < 0.01 ~ "**",
              p.value < 0.05 ~ "*",
              TRUE ~ ""
            )
          ) %>%
          dplyr::mutate(dplyr::across(where(is.numeric), ~signif(.x, 4)))
        
        return(list(success = TRUE, contrast_df = contrast_df))
      } else {
        return(list(success = FALSE, message = "No within-fiber contrasts available"))
      }
      
    }, error = function(e) {
      message("Fiber contrast error: ", e$message)
      list(success = FALSE, error = e$message)
    })
  })
  
  # ============================================================================
  # OUTPUT: COMBINED EMMEANS TABLE
  # ============================================================================
  
  output$combined_emm <- DT::renderDataTable({
    emm_res <- combined_emmeans_data()
    validate(need(!is.null(emm_res$emmeans_df), "Click Run to calculate emmeans"))
    
    df <- emm_res$emmeans_df
    
    DT::datatable(
      df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = "Blfrtip",
        columnDefs = list(
          list(className = 'dt-center', targets = "_all")
        )
      ),
      rownames = FALSE,
      caption = "Estimated Marginal Means for Combined Treatment Levels"
    ) %>%
      DT::formatRound(
        columns = c("emmean", "SE", "lower.CL", "upper.CL"),
        digits = 4
      )
  })
  
  # ============================================================================
  # OUTPUT: COMBINED EMMEANS PLOT
  # ============================================================================
  
  output$combined_emm_plot <- renderPlot({
    emm_res <- combined_emmeans_data()
    validate(need(!is.null(emm_res$emmeans_df), "Model not available"))
    
    df <- emm_res$emmeans_df
    
    # Add color coding by fiber type
    df <- df %>%
      dplyr::mutate(
        fiber_type = dplyr::case_when(
          stringr::str_detect(treatment_combined, "Cotton") ~ "Cotton",
          stringr::str_detect(treatment_combined, "Pet") ~ "PET",  # Changed from "PET"
          TRUE ~ "Control"
        ),
        treatment_type = dplyr::case_when(
          stringr::str_detect(treatment_combined, "Treated") ~ "Treated",
          stringr::str_detect(treatment_combined, "Untreated") ~ "Untreated",
          TRUE ~ "Control"
        )
      )
    
    # Check what fiber types actually exist in the data
    actual_fiber_types <- unique(df$fiber_type)
    
    # Create color palette only for existing types
    color_palette <- c()
    if ("Cotton" %in% actual_fiber_types) color_palette["Cotton"] <- "#2E7D32"
    if ("PET" %in% actual_fiber_types) color_palette["PET"] <- "#1565C0"
    if ("Control" %in% actual_fiber_types) color_palette["Control"] <- "#757575"
    
    # Reorder for better visualization
    df$treatment_combined <- forcats::fct_reorder(df$treatment_combined, df$emmean)
    
    ggplot2::ggplot(df, ggplot2::aes(x = emmean, y = treatment_combined, color = fiber_type)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = lower.CL, xmax = upper.CL),
        height = 0.3,
        linewidth = 1
      ) +
      ggplot2::scale_color_manual(values = color_palette) +  # Use dynamic palette
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(
        title = "Estimated Marginal Means with 95% Confidence Intervals",
        subtitle = paste("Reference:", input$combined_reference),
        x = "Estimated Mean (Outcome)",
        y = "Treatment Group",
        color = "Fiber Type"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 11),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "right"
      )
  })
  
  # ============================================================================
  # OUTPUT: COMBINED PAIRWISE TABLE
  # ============================================================================
  
  output$combined_pairwise <- DT::renderDataTable({
    filter_type <- input$combined_comparison_filter
    
    # Use custom fiber contrasts if selected
    if (filter_type == "fiber_control") {
      fiber_res <- combined_fiber_contrasts()
      validate(
        need(fiber_res$success, "Unable to create within-fiber contrasts"),
        need(!is.null(fiber_res$contrast_df), "No contrast data available")
      )
      df <- fiber_res$contrast_df
    } else {
      # Use regular pairwise comparisons
      pw_res <- combined_pairwise_data()
      validate(need(!is.null(pw_res$pairwise_df), "Run model first"))
      df <- pw_res$pairwise_df
    }
    
    alpha <- input$combined_alpha
    
    DT::datatable(
      df,
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        scrollY = "500px",
        dom = "Blfrtip",
        columnDefs = list(
          list(className = 'dt-center', targets = "_all")
        )
      ),
      rownames = FALSE,
      caption = paste0(
        "Pairwise Comparisons (α = ", alpha, ") - ",
        "Filter: ", filter_type
      )
    ) %>%
      DT::formatRound(
        columns = c("estimate", "SE", "t.ratio", "p.value"),
        digits = 4
      ) %>%
      DT::formatStyle(
        'p.value',
        backgroundColor = DT::styleInterval(
          c(0.001, 0.01, 0.05),
          c('#ffebee', '#fff3e0', '#fffde7', 'white')
        )
      ) %>%
      DT::formatStyle(
        'significant',
        fontWeight = 'bold',
        color = DT::styleEqual(
          c('***', '**', '*', ''),
          c('red', 'orange', 'blue', 'black')
        )
      )
  })
  
  # ============================================================================
  # OUTPUT: COMBINED PAIRWISE PLOT
  # ============================================================================
  
  output$combined_pairwise_plot <- renderPlot({
    pw_res <- combined_pairwise_data()
    validate(need(!is.null(pw_res$pairwise_df), "Run model first"))
    
    df <- pw_res$pairwise_df %>%
      dplyr::arrange(estimate) %>%
      dplyr::mutate(
        contrast = forcats::fct_reorder(contrast, estimate),
        sig_level = dplyr::case_when(
          p.value < 0.001 ~ "p < 0.001",
          p.value < 0.01 ~ "p < 0.01",
          p.value < 0.05 ~ "p < 0.05",
          TRUE ~ "Not significant"
        )
      )
    
    # Limit to top 30 comparisons if too many
    if (nrow(df) > 30) {
      df <- df %>%
        dplyr::slice_head(n = 30)
    }
    
    ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = contrast, color = sig_level)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = estimate - SE * 1.96, xmax = estimate + SE * 1.96),
        height = 0.2
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "p < 0.001" = "#D32F2F",
          "p < 0.01" = "#F57C00",
          "p < 0.05" = "#1976D2",
          "Not significant" = "#757575"
        )
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = "Pairwise Comparison Estimates",
        subtitle = paste("Filter:", input$combined_comparison_filter, "| Adjustment: Tukey HSD"),
        x = "Estimated Difference",
        y = "Contrast",
        color = "Significance"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 9),
        legend.position = "bottom"
      )
  })
  
  output$combined_diagnostics <- renderPlot({
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Run model first to see diagnostics"))
    
    # Create diagnostic plots
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    
    # 1. Residuals vs Fitted
    plot(
      fitted(mdl),
      residuals(mdl),
      xlab = "Fitted Values",
      ylab = "Residuals",
      main = "Residuals vs Fitted",
      pch = 16,
      col = alpha("steelblue", 0.6)
    )
    abline(h = 0, col = "red", lty = 2, lwd = 2)
    
    # 2. Normal Q-Q Plot
    qqnorm(
      residuals(mdl),
      main = "Normal Q-Q Plot",
      pch = 16,
      col = alpha("steelblue", 0.6)
    )
    qqline(residuals(mdl), col = "red", lwd = 2, lty = 2)
    
    # 3. Scale-Location (sqrt of standardized residuals vs fitted)
    sqrt_std_resid <- sqrt(abs(scale(residuals(mdl))))
    plot(
      fitted(mdl),
      sqrt_std_resid,
      xlab = "Fitted Values",
      ylab = expression(sqrt("|Standardized Residuals|")),
      main = "Scale-Location",
      pch = 16,
      col = alpha("steelblue", 0.6)
    )
    lines(
      lowess(fitted(mdl), sqrt_std_resid),
      col = "red",
      lwd = 2
    )
    
    # 4. Cook's Distance
    cooks_d <- cooks.distance(mdl)
    n <- length(cooks_d)
    plot(
      1:n,
      cooks_d,
      type = "h",
      xlab = "Observation Index",
      ylab = "Cook's Distance",
      main = "Cook's Distance (Influence)",
      col = ifelse(cooks_d > 4 / n, "red", "steelblue"),
      lwd = 2
    )
    abline(h = 4 / n, col = "red", lty = 2, lwd = 1.5)
    text(
      x = n * 0.7,
      y = max(cooks_d) * 0.9,
      labels = paste("Threshold:", round(4 / n, 4)),
      col = "red",
      cex = 0.8
    )
  })
  
  # -----------------------------------------------------------------------------
  # Download handler
  # -----------------------------------------------------------------------------
  output$download_combined_results <- downloadHandler(
    filename = function() paste0("combined_treatment_analysis_", Sys.Date(), ".txt"),
    content = function(file) {
      mdl <- combined_model()
      sink(file)
      cat("COMBINED TREATMENT ANALYSIS RESULTS\n")
      cat("Date:", as.character(Sys.Date()), "\n")
      cat("Dependent variable: outcome\n")
      cat("Reference level:", input$combined_reference, "\n\n")
      cat("MODEL SUMMARY\n")
      print(summary(mdl))
      cat("\nANOVA\n")
      print(stats::anova(mdl))
      sink()
    }
  )
  
  # ============================================================================
  # DYNAMIC ENDPOINT SELECTOR FOR COMBINED TREATMENT ANALYSIS
  # ============================================================================
  
  observe({
    req(input$active_dataset)  # CORRECT: active_dataset with underscore!
    
    cat("\n=== Updating combined_endpoint choices ===\n")
    cat("Active dataset:", input$active_dataset, "\n")
    
    choices <- character(0)
    
    if (identical(input$active_dataset, "assay")) {
      # Try to get assay endpoints from final_data
      choices <- tryCatch({
        req(final_data)
        sort(unique(final_data$assay_type))
      }, error = function(e) {
        cat("Error accessing final_data:", e$message, "\n")
        c("ACP", "AChE", "ALP", "Bradford", "CAT", "SOD")
      })
      cat("Assay endpoints:", paste(choices, collapse = ", "), "\n")
      
    } else if (identical(input$active_dataset, "physical")) {
      # Try to get physical endpoints from physical_master
      choices <- tryCatch({
        req(physical_master)
        sort(unique(physical_master$endpoint))
      }, error = function(e) {
        cat("Error accessing physical_master:", e$message, "\n")
        c("bci", "clearance_rate", "mf_counts", "respiration_rate")
      })
      cat("Physical endpoints:", paste(choices, collapse = ", "), "\n")
    }
    
    # Update the dropdown
    if (length(choices) > 0) {
      updateSelectInput(
        session, 
        "combined_endpoint", 
        choices = choices, 
        selected = choices[1]
      )
      cat("✓ Dropdown updated successfully with", length(choices), "choices\n")
    } else {
      updateSelectInput(
        session, 
        "combined_endpoint", 
        choices = c("No endpoints available" = "none"), 
        selected = "none"
      )
      cat("⚠ WARNING: No endpoints found!\n")
    }
  })
  
  # Status indicator
  output$combined_data_info <- renderText({
    if (is.null(input$combined_endpoint) || 
        identical(input$combined_endpoint, "none") ||
        identical(input$combined_endpoint, "loading")) {
      return("⚠ No data loaded")
    }
    
    data_count <- tryCatch({
      req(combined_base_data())
      nrow(combined_base_data())
    }, error = function(e) 0)
    
    if (data_count > 0) {
      paste0("✓ ", data_count, " observations loaded")
    } else {
      "⚠ No data for selected endpoint"
    }
  })
  
  # EMMEANS calculation for mixed effects models
  lmer_emmeans_results <- eventReactive(input$lmer_run_emmeans, {
    tryCatch({
      # Require that the lmer model exists first
      model <- lmer_model()
      validate(need(!is.null(model) && !inherits(model, "list"), "Run the mixed effects model first"))
      
      # Calculate dose value based on selection
      dose_val <- switch(input$lmer_emmeans_dose,
                         "0" = 0,           # Control level
                         "2" = log10(100),  # Low dose
                         "3" = log10(1000), # Medium dose
                         "4" = log10(10000)) # High dose
      
      is_ctrl <- ifelse(dose_val == 0, 1, 0)
      at_list <- list(dose_log10 = dose_val, is_control = is_ctrl)
      
      # FIXED: Better emmeans specifications to avoid convergence issues
      emmeans_spec <- switch(input$lmer_emmeans_by,
                             "treatment" = ~ chem_treatment,           # Simple: just treatment
                             "fiber_type" = ~ fiber_type,             # Simple: just fiber
                             "week" = ~ week,                         # Simple: just week
                             "treatment_week" = ~ chem_treatment + week,      # Additive, not interaction
                             "treatment_fiber" = ~ chem_treatment + fiber_type, # Additive, not interaction
                             "fiber_week" = ~ fiber_type + week,       # Additive, not interaction
                             "all_interactions" = ~ chem_treatment + fiber_type + week) # Additive only
      
      # ENHANCED: Try different approaches if emmeans fails
      emm <- tryCatch({
        # Primary approach: Use full specification
        emmeans(model, specs = emmeans_spec, at = at_list)
      }, error = function(e1) {
        message("Primary emmeans failed, trying simplified approach: ", e1$message)
        tryCatch({
          # Fallback 1: Remove the at= specification (use model's reference levels)
          emmeans(model, specs = emmeans_spec)
        }, error = function(e2) {
          message("Simplified emmeans failed, trying minimal approach: ", e2$message)
          # Fallback 2: Use simplest specification possible
          emmeans(model, specs = ~ chem_treatment)
        })
      })
      
      list(
        emmeans = emm,
        emmeans_df = as.data.frame(emm),
        specification_used = deparse(emmeans_spec)
      )
    }, error = function(e) {
      message("All emmeans approaches failed: ", e$message)
      list(error = e$message)
    })
  })
  
  # Mixed effects emmeans outputs  
  output$lmer_emmeans_table <- DT::renderDataTable({
    results <- lmer_emmeans_results()
    validate(need(!is.null(results$emmeans_df), "Click 'Calculate EM Means' for mixed effects"))
    
    df_emm <- results$emmeans_df %>%
      mutate(across(where(is.numeric), ~ signif(.x, 4)))
    
    DT::datatable(
      df_emm,
      options = list(pageLength = 10, scrollX = TRUE),
      filter = "top", rownames = FALSE
    ) %>%
      formatRound(columns = c("emmean", "SE", "lower.CL", "upper.CL"), digits = 4)
  })
  
  output$lmer_pairwise_table <- DT::renderDataTable({
    m <- lmer_model()
    validate(need(!is.null(m), "Run the mixed effects model first."))
    
    adj_method <- input$lmer_adjustment_method %||% "tukey"
    
    if (input$lmer_comparison_type == "did") {
      # DID analysis for mixed effects - simplified version
      tryCatch({
        dose_val <- as.numeric(input$lmer_emmeans_dose %||% "0")
        is_ctrl <- ifelse(dose_val == 0, 1, 0)
        
        emm_tw <- emmeans(m, 
                          specs = ~ chem_treatment | fiber_type + week,
                          at = list(dose_log10 = dose_val, is_control = is_ctrl))
        
        contrasts <- contrast(emm_tw, method = "pairwise", adjust = adj_method)
        dfp <- as.data.frame(contrasts)
        
      }, error = function(e) {
        dfp <- data.frame(
          Message = paste0("DiD Error for Mixed Effects: ", e$message),
          stringsAsFactors = FALSE
        )
      })
      
    } else {
      # Regular pairwise or treatment vs control for mixed effects
      emm_res <- lmer_emmeans_results()
      validate(need(!is.null(emm_res$emmeans), "Calculate EM Means first"))
      
      emm <- emm_res$emmeans
      
      if (input$lmer_comparison_type == "control") {
        contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "untreated", adjust = adj_method)
      } else {
        contrasts <- contrast(emm, method = "pairwise", adjust = adj_method)
      }
      
      dfp <- as.data.frame(contrasts)
    }
    
    # Add significance column and formatting
    if (!"p.value" %in% names(dfp)) {
      dfp$p.value <- NA_real_
    }
    
    dfp <- dfp %>%
      mutate(
        across(where(is.numeric), ~ signif(.x, 5)),
        Significant = case_when(
          is.na(p.value) ~ "Unable to compute",
          p.value < 0.001 ~ "Yes (***)",
          p.value < 0.01 ~ "Yes (**)", 
          p.value < 0.05 ~ "Yes (*)",
          TRUE ~ "No"
        ),
        Adjustment = adj_method
      )
    
    # Create colored DataTable with FIXED column names
    dt <- DT::datatable(
      dfp,
      options = list(pageLength = 15, scrollX = TRUE), 
      filter = "top", rownames = FALSE
    ) %>%
      DT::formatStyle(
        "Significant",
        backgroundColor = DT::styleEqual(
          c("Yes (***)", "Yes (**)", "Yes (*)", "No", "Unable to compute"),
          c("#28a745", "#28a745", "#28a745", "#ffffff", "#ffe6cc")  # Green for significant
        ),
        color = DT::styleEqual(
          c("Yes (***)", "Yes (**)", "Yes (*)", "No", "Unable to compute"),
          c("#ffffff", "#ffffff", "#ffffff", "#000000", "#000000")   # White text on green
        )
      ) %>%
      DT::formatRound(columns = "p.value", digits = 6)
    
    # FIXED: Check which columns actually exist before formatting
    col_names <- names(dfp)
    
    if ("estimate" %in% col_names) {
      dt <- dt %>% DT::formatRound(columns = "estimate", digits = 4)
    }
    if ("SE" %in% col_names) {
      dt <- dt %>% DT::formatRound(columns = "SE", digits = 4)
    }
    if ("t.ratio" %in% col_names) {
      dt <- dt %>% DT::formatRound(columns = "t.ratio", digits = 4)
    }
    if ("z.ratio" %in% col_names) {
      dt <- dt %>% DT::formatRound(columns = "z.ratio", digits = 4)
    }
    
    dt
  })
  
  # ============================================================================
  # 3. ENHANCED MIXED EFFECTS PLOT WITH BETTER ERROR HANDLING
  # ============================================================================
  
  output$lmer_emmeans_plot <- renderPlot({
    results <- lmer_emmeans_results()
    validate(need(!is.null(results$emmeans), "Calculate EM Means for mixed effects first"))
    validate(need(is.null(results$error), paste("Emmeans error:", results$error)))
    
    tryCatch({
      # SAME AS OLD APP - use built-in emmeans plot first
      plot(results$emmeans) + 
        theme_minimal() + 
        labs(title = "Mixed Effects Estimated Marginal Means with 95% CI",
             subtitle = paste("Dose level:", input$lmer_emmeans_dose, 
                              "| Adjustment:", input$lmer_adjustment_method),
             x = "Estimated Marginal Mean") +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          panel.grid.minor = element_blank()
        )
    }, error = function(e) {
      message("Built-in emmeans plot failed, using manual approach: ", e$message)
      
      # Enhanced fallback with better error handling
      df_plot <- results$emmeans_df
      
      # Remove rows with missing or invalid confidence intervals
      valid_rows <- !is.na(df_plot$emmean) & 
        !is.na(df_plot$lower.CL) & 
        !is.na(df_plot$upper.CL) &
        is.finite(df_plot$lower.CL) & 
        is.finite(df_plot$upper.CL)
      
      if (sum(valid_rows) == 0) {
        # Fallback: create simple point plot without confidence intervals
        ggplot(df_plot, aes(x = emmean, y = rownames(df_plot))) +
          geom_point(size = 4, color = "#e74c3c") +
          theme_minimal() +
          labs(title = "Mixed Effects Estimated Marginal Means (No CI due to model issues)",
               subtitle = paste("Dose level:", input$lmer_emmeans_dose, "| Some confidence intervals invalid"),
               x = "Estimated Marginal Mean", y = "Groups") +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11, color = "red"),
            axis.text.y = element_text(size = 10)
          )
      } else {
        # Use valid rows only with FIXED geom_errorbar
        df_clean <- df_plot[valid_rows, ]
        
        ggplot(df_clean, aes(x = emmean, y = factor(rownames(df_clean)))) +
          geom_point(size = 3, color = "#2c3e50") +
          # FIXED: Use geom_errorbar with orientation instead of deprecated geom_errorbarh
          geom_errorbar(aes(xmin = pmax(lower.CL, emmean - 3*SE), 
                            xmax = pmin(upper.CL, emmean + 3*SE)), 
                        width = 0.2, orientation = "y", color = "#3498db") +
          theme_minimal() +
          labs(title = "Mixed Effects Estimated Marginal Means with 95% CI",
               subtitle = paste("Dose level:", input$lmer_emmeans_dose,
                                "| Adjustment:", input$lmer_adjustment_method),
               x = "Estimated Marginal Mean", y = "Groups") +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            axis.text.y = element_text(size = 10),
            panel.grid.minor = element_blank()
          )
      }
    })
  })
  
  # ============================================================================
  #                           REGRESSION OUTPUTS
  # ============================================================================
  
  output$regression_model_summary <- renderPrint({
    m <- regression_model()
    validate(need(!is.null(m) && !inherits(m, "list"), "Click 'Run Regression Model' to fit the model."))
    summary(m)
  })
  
  output$regression_anova <- DT::renderDataTable({
    m <- regression_model()
    validate(need(!is.null(m) && !inherits(m, "list"), ""))
    a <- broom::tidy(car::Anova(m))
    DT::datatable(a, options = list(pageLength = 10), filter = "top", rownames = FALSE)
  })
  
  # ============================================================================
  #                          MIXED EFFECTS OUTPUTS
  # ============================================================================
  
  output$lmer_model_summary <- renderPrint({
    m <- lmer_model()
    validate(need(!is.null(m) && !inherits(m, "list"), "Click 'Run Mixed Model' to fit the model."))
    summary(m)
  })
  
  output$lmer_anova <- DT::renderDataTable({
    m <- lmer_model()
    validate(need(!is.null(m) && !inherits(m, "list"), ""))
    a <- broom.mixed::tidy(car::Anova(m))
    DT::datatable(a, options = list(pageLength = 10), filter = "top", rownames = FALSE)
  })
  
  # Model comparison output 
  output$model_compare <- DT::renderDataTable({
    # Simple model comparison between lm and lmer
    lm_model <- regression_model()
    lmer_model_obj <- lmer_model()
    
    validate(need(!is.null(lm_model) && !inherits(lm_model, "list") && 
                    !is.null(lmer_model_obj) && !inherits(lmer_model_obj, "list"), 
                  "Run both models to compare"))
    
    # Extract model comparison metrics
    lm_aic <- AIC(lm_model)
    lmer_aic <- AIC(lmer_model_obj)
    
    comparison_df <- data.frame(
      Model = c("Linear Model (lm)", "Mixed Effects (lmer)"),
      AIC = c(lm_aic, lmer_aic),
      BIC = c(BIC(lm_model), BIC(lmer_model_obj)),
      LogLik = c(as.numeric(logLik(lm_model)), as.numeric(logLik(lmer_model_obj))),
      df = c(length(coef(lm_model)), length(fixef(lmer_model_obj))),
      Better_AIC = c(lm_aic < lmer_aic, lmer_aic < lm_aic)
    )
    
    DT::datatable(comparison_df, options = list(pageLength = 5, dom = 't'), rownames = FALSE) %>%
      formatRound(columns = c("AIC", "BIC", "LogLik"), digits = 2)
  })
  
  
  # ============================================================================
  # RECOVERY SERVER (standalone, no helpers)
  # ============================================================================
  
  # Reactive state
  recovery_state <- reactiveVal(NULL)
  
  # Mini utilities (local to Recovery, no external helpers)
  canon_fiber <- function(x) {
    xl <- tolower(trimws(as.character(x)))
    dplyr::case_when(
      xl %in% c("cot", "cotton") ~ "cotton",
      xl %in% c("pet", "polyester", "poly") ~ "pet",
      TRUE ~ xl
    )
  }
  
  canon_treatment <- function(x) {
    xl <- tolower(trimws(as.character(x)))
    dplyr::case_when(
      xl %in% c("untreated", "control", "ctrl") ~ "untreated",
      xl %in% c("treated") ~ "treated",
      TRUE ~ NA_character_
    )
  }
  
  canon_conc_labels <- function(x) {
    v <- suppressWarnings(as.numeric(as.character(x)))
    allowed <- c(0, 100, 1000, 10000)
    snapped <- vapply(v, function(y) {
      if (!is.finite(y)) return(NA_real_)
      idx <- which.min(abs(allowed - y))
      tol <- max(1, 0.01 * allowed[idx])
      if (abs(allowed[idx] - y) <= tol) allowed[idx] else NA_real_
    }, numeric(1))
    as.character(as.integer(snapped))
  }
  
  nlev2 <- function(x) {
    length(unique(stats::na.omit(as.character(x))))
  }
  
  # Populate endpoints on dataset change
  observe({
    if (identical(input$recovery_dataset, "assay")) {
      req(exists("final_data"))
      eps <- sort(unique(as.character(final_data$assay_type)))
    } else {
      req(exists("physical_master"))
      eps <- sort(unique(as.character(physical_master$endpoint)))
    }
    updateSelectInput(session, "recovery_endpoint", choices = eps, selected = eps[1])
  })
  
  # Main recovery analysis
  observeEvent(input$run_recovery, {
    message("=== RECOVERY START (scratch) ===")
    recovery_state(NULL)
    
    wk_base <- 5L; wk_reco <- 6L
    ds <- input$recovery_dataset
    ep <- input$recovery_endpoint
    
    # Bring in base data
    if (identical(ds, "assay")) {
      req(exists("final_data"))
      df <- final_data %>%
        dplyr::filter(.data$assay_type == ep) %>%
        dplyr::mutate(
          outcome = calculated_concentration,
          week = as.integer(week),
          fiber_type = canon_fiber(fiber_type),
          chem_treatment = canon_treatment(treatment),
          sample_type = as.character(sample_type),
          fiber_concentration_label = canon_conc_labels(fiber_concentration)
        )
    } else {
      req(exists("physical_master"))
      df <- physical_master %>%
        dplyr::filter(.data$endpoint == ep) %>%
        dplyr::mutate(
          outcome = value,
          week = as.integer(week),
          fiber_type = canon_fiber(fiber_type),
          chem_treatment = canon_treatment(treatment),
          sample_type = dplyr::coalesce(as.character(sample_type), as.character(tissue_type)),
          fiber_concentration_label = canon_conc_labels(fiber_concentration)
        )
    }
    message("[RECOVERY] Base rows: ", nrow(df))
    
    # Filters
    if (!is.null(input$recovery_fiber_types) && length(input$recovery_fiber_types))
      df <- dplyr::filter(df, fiber_type %in% input$recovery_fiber_types)
    
    if (!is.null(input$recovery_treatments) && length(input$recovery_treatments))
      df <- dplyr::filter(df, chem_treatment %in% input$recovery_treatments)
    
    if (identical(ds, "assay")) {
      if (!is.null(input$recovery_samples) && length(input$recovery_samples))
        df <- dplyr::filter(df, tolower(sample_type) %in% tolower(input$recovery_samples))
    } else {
      if (!is.null(input$recovery_tissues) && length(input$recovery_tissues))
        df <- dplyr::filter(df, tolower(sample_type) %in% tolower(input$recovery_tissues))
    }
    
    # Concentration (discrete labels)
    if ("fiber_concentration_label" %in% names(df)) {
      if (identical(ds, "assay")) {
        if (!is.null(input$recovery_concentration) && length(input$recovery_concentration))
          df <- dplyr::filter(df, fiber_concentration_label %in% input$recovery_concentration)
      } else {
        if (!is.null(input$recovery_concentration_physical) && length(input$recovery_concentration_physical))
          df <- dplyr::filter(df, fiber_concentration_label %in% input$recovery_concentration_physical)
      }
      df <- df %>%
        dplyr::mutate(fiber_concentration = suppressWarnings(as.numeric(fiber_concentration_label)))
    }
    
    message("[RECOVERY] After filters: ", nrow(df))
    validate(need(nrow(df) > 0, "No data after filters"))
    
    # Week subset
    df <- df %>%
      dplyr::filter(week %in% c(wk_base, wk_reco)) %>%
      dplyr::mutate(
        is_recovery = factor(ifelse(week == wk_reco, "Week6", "Week5"), levels = c("Week5","Week6")),
        fiber_type = factor(fiber_type),
        chem_treatment = factor(chem_treatment),
        sample_type = factor(sample_type)
      ) %>%
      droplevels()
    
    validate(need(dplyr::n_distinct(df$is_recovery) == 2, "Both Week 5 and Week 6 required"))
    
    # Hierarchical week coverage
    present <- function(v) v %in% names(df) && nlev2(df[[v]]) > 1
    strict_vars <- c(if (present("fiber_type")) "fiber_type",
                     if (present("chem_treatment")) "chem_treatment",
                     if (present("sample_type") && identical(ds,"assay")) "sample_type")
    candidates <- list(
      strict_vars,
      setdiff(strict_vars, "sample_type"),
      setdiff(strict_vars, c("sample_type","fiber_type")),
      character(0)
    )
    
    picked <- NULL
    for (gv in candidates) {
      df_try <- if (length(gv)) {
        df %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(gv))) %>%
          dplyr::filter(dplyr::n_distinct(is_recovery) == 2) %>%
          dplyr::ungroup()
      } else {
        if (dplyr::n_distinct(df$is_recovery) == 2) df else df[0,]
      }
      if (nrow(df_try) > 0) { picked <- gv; df <- df_try; break }
    }
    msg_pick <- if (is.null(picked) || length(picked) == 0) "pooled" else paste(picked, collapse="×")
    message("[RECOVERY] Week coverage by: ", msg_pick, " | rows: ", nrow(df))
    validate(need(nrow(df) > 0, "No strata retain both weeks; widen filters or include another dose"))
    
    # Estimable formula
    include_trt  <- nlev2(df$chem_treatment) >= 2
    include_fib  <- nlev2(df$fiber_type)    >= 2
    include_samp <- nlev2(df$sample_type)   >= 2
    
    rhs <- c("is_recovery")
    if (include_trt)  rhs <- c(rhs, "chem_treatment", "is_recovery:chem_treatment")
    if (include_fib)  rhs <- c(rhs, "fiber_type", "is_recovery:fiber_type")
    if (include_samp) rhs <- c(rhs, "sample_type", "is_recovery:sample_type")
    
    form <- stats::as.formula(paste("outcome ~", paste(rhs, collapse = " + ")))
    message("[RECOVERY] Formula: ", deparse(form))
    
    mdl <- tryCatch(stats::lm(form, data = df), error = function(e) e)
    if (inherits(mdl, "error")) {
      recovery_state(list(
        error = paste("Model fitting failed:", mdl$message),
        data = df, formula = form, endpoint = ep,
        dataset = ds, baseline_week = wk_base, recovery_week = wk_reco
      ))
      return()
    }
    
    # EMMeans
    emm_spec <- if (include_trt && include_fib && include_samp && identical(ds, "assay")) {
      ~ is_recovery | chem_treatment + fiber_type + sample_type
    } else if (include_trt && include_fib) {
      ~ is_recovery | chem_treatment + fiber_type
    } else if (include_trt) {
      ~ is_recovery | chem_treatment
    } else if (include_fib) {
      ~ is_recovery | fiber_type
    } else if (include_samp) {
      ~ is_recovery | sample_type
    } else {
      ~ is_recovery
    }
    
    emm <- tryCatch(emmeans::emmeans(mdl, emm_spec),
                    error = function(e) { message("[RECOVERY] emmeans failed: ", e$message); NULL })
    delta <- if (!is.null(emm)) {
      tryCatch(emmeans::contrast(emm, "revpairwise"), error = function(e) NULL)
    } else NULL
    emm_df <- if (!is.null(emm)) tryCatch(as.data.frame(emm), error = function(e) NULL) else NULL
    
    recovery_state(list(
      model = mdl,
      data = df,
      endpoint = ep,
      dataset = ds,
      baseline_week = wk_base,
      recovery_week = wk_reco,
      emmeans = list(emmeans = emm, emmeans_df = emm_df, delta = delta),
      formula = form
    ))
    message("=== RECOVERY DONE ===")
  })
  
  # ---- Outputs ----
  output$recovery_model_summary <- renderPrint({
    res <- recovery_state()
    validate(need(!is.null(res), "Click 'Run Recovery Analysis'"))
    validate(need(is.null(res$error), res$error))
    cat("Dataset:", res$dataset, "\n")
    cat("Endpoint:", res$endpoint, "\n")
    cat("Weeks: ", res$baseline_week, " vs ", res$recovery_week, "\n", sep = "")
    cat("Formula:", deparse(res$formula), "\n\n")
    print(summary(res$model))
  })
  
  output$recovery_comparisons <- DT::renderDataTable({
    res <- recovery_state()
    validate(need(!is.null(res), ""))
    validate(need(is.null(res$error), res$error))
    validate(need(!is.null(res$emmeans$delta), "No pairwise Week6-Week5 contrasts available"))
    DT::datatable(as.data.frame(res$emmeans$delta), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$recovery_plot <- renderPlot({
    res <- recovery_state()
    validate(need(!is.null(res), ""))
    validate(need(is.null(res$error), res$error))
    df <- res$data
    ggplot2::ggplot(df, ggplot2::aes(x = is_recovery, y = outcome, color = fiber_type, linetype = chem_treatment)) +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 3) +
      ggplot2::stat_summary(fun.data = ~ c(y = mean(.), ymin = mean(.) - sd(.)/sqrt(length(.)),
                                           ymax = mean(.) + sd(.)/sqrt(length(.))), geom = "errorbar", width = 0.2) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Recovery (Week6 vs Week5)", x = "Time point", y = "Outcome")
  })
  
  output$recovery_emmeans_table <- DT::renderDataTable({
    res <- recovery_state()
    validate(need(!is.null(res), ""))
    validate(need(is.null(res$error), res$error))
    validate(need(!is.null(res$emmeans$emmeans_df), "No EM Means grid available"))
    DT::datatable(res$emmeans$emmeans_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$recovery_emmeans_plot <- renderPlot({
    res <- recovery_state()
    validate(need(!is.null(res), ""))
    validate(need(is.null(res$error), res$error))
    if (!is.null(res$emmeans$emmeans)) {
      return(plot(res$emmeans$emmeans) + ggplot2::theme_minimal() +
               ggplot2::labs(title = "Estimated Marginal Means", x = "EMMean"))
    }
    validate(need(FALSE, "No emmeans object to plot"))
  })
  
  # ---- Translocation analysis (mf_counts only) ----
  translocation_analysis <- eventReactive(input$run_translocation, {
    df <- physical_master %>%
      create_enhanced_treatment_categories() %>%       # helpers.R
      dplyr::filter(endpoint == "mf_counts")
    
    # summary by fiber, week, tissue, treatment
    summary_df <- df %>%
      dplyr::group_by(fiber_type, week, tissue_type, treatment_category) %>%
      dplyr::summarise(
        mean_count = mean(value, na.rm = TRUE),
        sd_count   = sd(value, na.rm = TRUE),
        n          = dplyr::n(),
        prop_nonzero = mean(value > 0, na.rm = TRUE),
        .groups = "drop"
      )
    summary_df
  })
  
  # translocation table
  output$translocation_summary <- DT::renderDataTable({
    tbl <- translocation_analysis()
    validate(need(nrow(tbl) > 0, "Click 'Run Translocation Analysis'"))
    DT::datatable(tbl, options = list(pageLength = 15, scrollX = TRUE), filter = "top", rownames = FALSE)
  })
  
  # translocation bar plot
  output$translocation_plot <- renderPlot({
    tbl <- translocation_analysis()
    validate(need(nrow(tbl) > 0, "Run Translocation Analysis"))
    ggplot(tbl, aes(x = factor(week), y = mean_count, fill = tissue_type)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
      facet_grid(fiber_type ~ treatment_category) +
      theme_minimal() +
      labs(title = "Microfiber Translocation (mf_counts)",
           x = "Week", y = "Mean Count", fill = "Tissue")
  })
  
  # tissue comparison plot (proportion non-zero)
  output$tissue_comparison_plot <- renderPlot({
    tbl <- translocation_analysis()
    validate(need(nrow(tbl) > 0, "Run Translocation Analysis"))
    ggplot(tbl, aes(x = tissue_type, y = prop_nonzero, fill = treatment_category)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.85) +
      facet_wrap(~ week) +
      theme_minimal() +
      labs(title = "Proportion of Non-zero Counts by Tissue",
           x = "Tissue", y = "Proportion > 0", fill = "Treatment") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # ============================================================================
  #                  TRANSLOCATION OUTPUTS
  # ============================================================================
  
  
  output$translocation_summary <- DT::renderDataTable({
    analysis <- translocation_analysis()
    validate(need(!is.null(analysis), "Click 'Run Translocation Analysis' to generate summary."))
    
    DT::datatable(analysis, options = list(pageLength = 15), filter = "top", rownames = FALSE)
  })
  
  output$translocation_plot <- renderPlot({
    analysis <- translocation_analysis()
    validate(need(!is.null(analysis), ""))
    
    ggplot(analysis, aes(x = week, y = mean_count, fill = tissue_type)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.7) +
      facet_grid(fiber_type ~ treatment_category) +
      theme_minimal() +
      labs(title = "Microfiber Translocation Across Tissues",
           x = "Week", y = "Mean Count", fill = "Tissue Type") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$tissue_comparison_plot <- renderPlot({
    analysis <- translocation_analysis()
    validate(need(!is.null(analysis), ""))
    
    ggplot(analysis, aes(x = tissue_type, y = prop_nonzero, fill = treatment_category)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.7) +
      facet_wrap(~ week) +
      theme_minimal() +
      labs(title = "Proportion of Non-zero Counts by Tissue",
           x = "Tissue Type", y = "Proportion Non-zero", fill = "Treatment") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # ============================================================================
  #                             DID SERVER
  # ============================================================================
  
  # Small local utilities (reuse pattern from Recovery)
  did_canon_fiber <- function(x) {
    xl <- tolower(trimws(as.character(x)))
    dplyr::case_when(
      xl %in% c("cot", "cotton") ~ "cotton",
      xl %in% c("pet", "polyester", "poly") ~ "pet",
      TRUE ~ xl
    )
  }
  
  did_canon_treatment <- function(x) {
    xl <- tolower(trimws(as.character(x)))
    dplyr::case_when(
      xl %in% c("untreated", "control", "ctrl") ~ "untreated",
      xl %in% c("treated") ~ "treated",
      TRUE ~ NA_character_
    )
  }
  
  did_canon_conc_labels <- function(x) {
    v <- suppressWarnings(as.numeric(as.character(x)))
    allowed <- c(0, 100, 1000, 10000)
    snapped <- vapply(v, function(y) {
      if (!is.finite(y)) return(NA_real_)
      idx <- which.min(abs(allowed - y))
      tol <- max(1, 0.01 * allowed[idx])
      if (abs(allowed[idx] - y) <= tol) allowed[idx] else NA_real_
    }, numeric(1))
    as.character(as.integer(snapped))
  }
  
  nlev2 <- function(x) length(unique(stats::na.omit(as.character(x))))
  
  # Populate endpoint and week choices
  observe({
    if (identical(input$did_dataset, "assay")) {
      req(exists("final_data"))
      eps <- sort(unique(as.character(final_data$assay_type)))
      wks <- sort(unique(as.integer(final_data$week)))
    } else {
      req(exists("physical_master"))
      eps <- sort(unique(as.character(physical_master$endpoint)))
      wks <- sort(unique(as.integer(physical_master$week)))
    }
    updateSelectInput(session, "did_endpoint", choices = eps, selected = eps[1])
    updateSelectInput(session, "did_week_single", choices = wks, selected = wks[1])
    updateSelectInput(session, "did_week_single2", choices = wks, selected = wks[1])
    updateSelectInput(session, "did_week_a", choices = wks, selected = wks[1])
    updateSelectInput(session, "did_week_b", choices = rev(wks)[1])
  })
  
  did_state <- reactiveVal(NULL)
  
  observeEvent(input$run_did, {
    message("=== DID START (patched) ===")
    did_state(NULL)
    
    ds   <- input$did_dataset
    ep   <- input$did_endpoint
    mode <- input$did_mode
    
    # ---- Base data (normalize sample_type) ----
    if (identical(ds, "assay")) {
      req(exists("final_data"))
      df0 <- final_data %>%
        dplyr::filter(.data$assay_type == ep) %>%
        dplyr::mutate(
          outcome = calculated_concentration,
          week = as.integer(week),
          fiber_type = did_canon_fiber(fiber_type),
          chem_treatment = did_canon_treatment(treatment),
          # Normalize sample_type once at the base
          sample_type = tolower(trimws(as.character(sample_type))),
          fiber_concentration_label = did_canon_conc_labels(fiber_concentration)
        )
    } else {
      req(exists("physical_master"))
      df0 <- physical_master %>%
        dplyr::filter(.data$endpoint == ep) %>%
        dplyr::mutate(
          outcome = value,
          week = as.integer(week),
          fiber_type = did_canon_fiber(fiber_type),
          chem_treatment = did_canon_treatment(treatment),
          # Coalesce and normalize for physical data too
          sample_type = tolower(trimws(dplyr::coalesce(as.character(sample_type), as.character(tissue_type)))),
          fiber_concentration_label = did_canon_conc_labels(fiber_concentration)
        )
    }
    
    # Keep an unfiltered copy for smart pooling
    df <- df0
    
    # ---- Sample/tissue filter (apply, but allow smart fallback) ----
    applied_sample_filter <- FALSE
    if (identical(ds, "assay")) {
      if (!is.null(input$did_samples) && length(input$did_samples)) {
        sel_samples <- tolower(trimws(input$did_samples))
        df <- dplyr::filter(df, sample_type %in% sel_samples)
        applied_sample_filter <- TRUE
      }
    } else {
      if (!is.null(input$did_tissues) && length(input$did_tissues)) {
        sel_tissues <- tolower(trimws(input$did_tissues))
        df <- dplyr::filter(df, sample_type %in% sel_tissues)
        applied_sample_filter <- TRUE
      }
    }
    
    # ---- Dose filter (always safe)
    if ("fiber_concentration_label" %in% names(df) &&
        !is.null(input$did_concentration) && length(input$did_concentration)) {
      df <- dplyr::filter(df, fiber_concentration_label %in% input$did_concentration) %>%
        dplyr::mutate(fiber_concentration = suppressWarnings(as.numeric(fiber_concentration_label)))
    }
    
    # ---- Fiber/treatment filters with mode-aware logic
    # Across-fiber requirements
    needs_both_fibers_for_mode <-
      identical(mode, "across_fiber") || identical(mode, "across_treatment")
    # Across-fiber (DiD) and across-week need both treatments
    needs_both_trt_for_mode <-
      (identical(mode, "across_fiber") && identical(input$did_across_fiber_contrast, "did")) ||
      identical(mode, "across_week")
    
    # Fiber filter: only apply when not required to have both fibers
    if (!needs_both_fibers_for_mode) {
      if (!is.null(input$did_fibers_filter) && length(input$did_fibers_filter))
        df <- dplyr::filter(df, fiber_type %in% input$did_fibers_filter)
    }
    
    # Treatment filter:
    if (identical(mode, "across_fiber") && identical(input$did_across_fiber_contrast, "single_trt")) {
      # Single-treatment fiber contrast: enforce the chosen treatment
      df <- dplyr::filter(df, chem_treatment == input$did_single_treatment)
    } else if (!needs_both_trt_for_mode) {
      # Modes where a subset of treatments is OK
      if (!is.null(input$did_treatments_filter) && length(input$did_treatments_filter))
        df <- dplyr::filter(df, chem_treatment %in% input$did_treatments_filter)
    }
    validate(need(nrow(df) > 0, "No data after DiD filters."))
    
    # ---- Mode-specific subsetting
    if (identical(mode, "across_fiber")) {
      wk <- as.integer(input$did_week_single)
      df <- dplyr::filter(df, week == wk, fiber_type %in% c(input$did_fiber_a, input$did_fiber_b))
    } else if (identical(mode, "across_week")) {
      wk_a <- as.integer(input$did_week_a)
      wk_b <- as.integer(input$did_week_b)
      df <- dplyr::filter(df, week %in% c(wk_a, wk_b), fiber_type == input$did_fiber_single)
    } else { # across_treatment
      wk2 <- as.integer(input$did_week_single2)
      df <- dplyr::filter(df, week == wk2, fiber_type %in% c("cotton","pet"),
                          chem_treatment %in% c(input$did_trt_a, input$did_trt_b))
    }
    validate(need(nrow(df) > 0, "No rows after mode-specific subsetting."))
    
    # ---- Smart sample pooling: if coverage is insufficient, revert the sample filter
    has_both <- function(v) v %in% names(df) && nlev2(df[[v]]) >= 2
    need_fibers <- needs_both_fibers_for_mode
    need_trt    <- needs_both_trt_for_mode
    
    coverage_ok <- (!need_fibers || has_both("fiber_type")) &&
      (!need_trt    || has_both("chem_treatment"))
    
    if (!coverage_ok && applied_sample_filter) {
      # Rebuild df WITHOUT the sample/tissue filter, then reapply the same
      # dose + mode-specific subsetting to rescue contrast estimability.
      df <- df0
      if ("fiber_concentration_label" %in% names(df) &&
          !is.null(input$did_concentration) && length(input$did_concentration)) {
        df <- dplyr::filter(df, fiber_concentration_label %in% input$did_concentration) %>%
          dplyr::mutate(fiber_concentration = suppressWarnings(as.numeric(fiber_concentration_label)))
      }
      if (identical(mode, "across_fiber")) {
        wk <- as.integer(input$did_week_single)
        df <- dplyr::filter(df, week == wk, fiber_type %in% c(input$did_fiber_a, input$did_fiber_b))
        if (identical(input$did_across_fiber_contrast, "single_trt")) {
          df <- dplyr::filter(df, chem_treatment == input$did_single_treatment)
        }
      } else if (identical(mode, "across_week")) {
        wk_a <- as.integer(input$did_week_a); wk_b <- as.integer(input$did_week_b)
        df <- dplyr::filter(df, week %in% c(wk_a, wk_b), fiber_type == input$did_fiber_single)
      } else {
        wk2 <- as.integer(input$did_week_single2)
        df <- dplyr::filter(df, week == wk2, fiber_type %in% c("cotton","pet"),
                            chem_treatment %in% c(input$did_trt_a, input$did_trt_b))
      }
      # Re-evaluate coverage
      coverage_ok <- (!need_fibers || has_both("fiber_type")) &&
        (!need_trt    || has_both("chem_treatment"))
      validate(need(coverage_ok, "Not enough data to estimate the requested contrast; widen filters."))
    }
    
    # ---- Factors
    df <- df %>%
      dplyr::mutate(
        week = factor(week),
        fiber_type = factor(fiber_type, levels = c("cotton","pet")),
        chem_treatment = factor(chem_treatment, levels = c("untreated","treated")),
        sample_type = factor(sample_type)
      ) %>% droplevels()
    
    # ---- Model
    form <- switch(
      mode,
      "across_fiber" = outcome ~ chem_treatment * fiber_type + sample_type,
      "across_week"  = outcome ~ chem_treatment * week + sample_type,
      "across_treatment" = outcome ~ fiber_type * chem_treatment + sample_type
    )
    mdl <- tryCatch(stats::lm(form, data = df), error = function(e) e)
    validate(need(!inherits(mdl, "error"), paste("Model failed:", mdl$message)))
    
    # ---- EMMeans and contrasts (fixed branching) ----
    if (identical(mode, "across_fiber")) {
      if (identical(input$did_across_fiber_contrast, "did")) {
        # True DiD: (treated − untreated) within each fiber, then Fiber A − Fiber B
        emm <- emmeans::emmeans(mdl, ~ chem_treatment | fiber_type)
        within_fiber <- emmeans::contrast(
          emm,
          method = list(trt_minus_ctrl = c(treated = 1, untreated = -1))
        )
        a <- input$did_fiber_a
        b <- input$did_fiber_b
        did <- emmeans::contrast(
          within_fiber,
          method = list(DiD = stats::setNames(c(1, -1), c(a, b))),
          by = NULL
        )
        did_tbl <- as.data.frame(did)
        emm_tbl <- as.data.frame(emm)
      } else {
        # Single-treatment fiber contrast at selected week: PET − Cotton within selected chem_treatment
        emm <- emmeans::emmeans(mdl, ~ fiber_type | chem_treatment)
        sel <- input$did_single_treatment
        fiber_con <- emmeans::contrast(
          emm,
          method = list(pet_minus_cotton = c(pet = 1, cotton = -1)),
          by = "chem_treatment"
        )
        fiber_con_df <- as.data.frame(fiber_con)
        if (!any(fiber_con_df$chem_treatment == sel)) {
          did_state(list(
            model = mdl, data = df, mode = mode,
            emmeans_tbl = as.data.frame(emm),
            did_tbl = data.frame(
              Message = sprintf("No %s rows available for fiber contrast at this week; widen filters.", sel),
              stringsAsFactors = FALSE
            )
          ))
          message("[DiD] Single-treatment contrast: no rows for selected treatment after filtering")
          return()
        }
        did_tbl <- fiber_con_df[fiber_con_df$chem_treatment == sel, , drop = FALSE]
        emm_tbl <- as.data.frame(emm)
      }
      
    } else if (identical(mode, "across_week")) {
      # (treated − untreated) within fiber at Week B minus Week A
      emm <- emmeans::emmeans(mdl, ~ chem_treatment | week)
      within_week <- emmeans::contrast(
        emm,
        method = list(trt_minus_ctrl = c(treated = 1, untreated = -1))
      )
      wa <- as.character(input$did_week_a)
      wb <- as.character(input$did_week_b)
      did <- emmeans::contrast(
        within_week,
        method = list(DiD = stats::setNames(c(1, -1), c(wb, wa))),
        by = NULL
      )
      did_tbl <- as.data.frame(did)
      emm_tbl <- as.data.frame(emm)
      
    } else { # across_treatment
      # (PET − Cotton) within Treatment A minus within Treatment B
      emm <- emmeans::emmeans(mdl, ~ fiber_type | chem_treatment)
      within_trt <- emmeans::contrast(
        emm,
        method = list(pet_minus_cotton = stats::setNames(c(1, -1), c("pet","cotton")))
      )
      ta <- input$did_trt_a
      tb <- input$did_trt_b
      did <- emmeans::contrast(
        within_trt,
        method = list(DiD = stats::setNames(c(1, -1), c(ta, tb))),
        by = NULL
      )
      did_tbl <- as.data.frame(did)
      emm_tbl <- as.data.frame(emm)
    }
    
    did_state(list(model = mdl, data = df, emmeans_tbl = emm_tbl, did_tbl = did_tbl, mode = mode))
    message("=== DID DONE (patched) ===")
  
  output$did_summary <- renderPrint({
    res <- did_state()
    validate(need(!is.null(res), "Click 'Run DiD Analysis'"))
    if (!is.null(res$error)) {
      cat("Error:", res$error, "\n")
      return(invisible(NULL))
    }
    cat("DiD Mode:", res$mode, "\n")
    if (!is.null(res$formula)) {
      cat("Formula:", deparse(res$formula), "\n")
    }
    cat("Rows used:", if (!is.null(res$data)) nrow(res$data) else 0, "\n\n")
    if (!is.null(res$model)) {
      print(summary(res$model))
    } else {
      cat("No model fit (insufficient coverage or filters too narrow).\n")
    }
  })
  
  output$did_table <- DT::renderDataTable({
    res <- did_state()
    validate(need(!is.null(res), ""))
    tbl <- res$did_tbl
    validate(need(!is.null(tbl) && nrow(tbl) > 0, "No DiD rows for current settings."))
    DT::datatable(tbl, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$did_emm_table <- DT::renderDataTable({
    res <- did_state()
    validate(need(!is.null(res), ""))
    tbl <- res$emmeans_tbl
    validate(need(!is.null(tbl) && nrow(tbl) > 0, "No EMMeans grid for current settings."))
    DT::datatable(tbl, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$did_plot <- renderPlot({
    res <- did_state()
    validate(need(!is.null(res), ""))
    dfp <- res$data
    validate(need(!is.null(dfp) && nrow(dfp) > 0, "No data to plot."))
    ggplot2::ggplot(dfp, ggplot2::aes(x = chem_treatment, y = outcome, color = fiber_type)) +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 3,
                            position = ggplot2::position_dodge(width = 0.35)) +
      ggplot2::stat_summary(
        fun.data = ~ c(y = mean(.), ymin = mean(.) - sd(.)/sqrt(length(.)),
                       ymax = mean(.) + sd(.)/sqrt(length(.))),
        geom = "errorbar", width = 0.2,
        position = ggplot2::position_dodge(width = 0.35)
      ) +
      ggplot2::facet_wrap(~ week, nrow = 1, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "DiD visualization", x = "Treatment", y = "Outcome")
  })
  
  # Memory cleanup
  onStop(function() {
    rm(list = ls())
    gc()
  })
})
}
# ============================================================================
# RUN APPLICATION  
# ============================================================================

shinyApp(ui = ui, server = server)