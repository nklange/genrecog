source("R/helper.R")

ui <-
  navbarPage(title = "RECOGNITION MEMORY DATA",

             theme = "bootstrap.css",
             shinyjs::useShinyjs(),

             tabPanel(title = "Individual Datasets",
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,

                          radioButtons("radio", h3("Datasets..."),
                                       choices = list("included in meta-analysis" = 1, "not included" = 3),
                                       selected = 1),
                          conditionalPanel(
                            condition = "input.radio == 1",
                            title = "Select dataset",
                            selectInput("exp1", "Select experiment:",
                                        choices = list("AP2007_e1",   "AP2007_e2",   "AP2007_e3",   "APH2016_e1",  "BKSR2013_e1", "BSG2014_e1",
                                                       "BTL2013_e1",  "CSR2015_e1",  "D2007_e1a",   "D2007_e1b",   "DR2012_e1b",  "DW2020_e1",
                                                       "FBH2013_e1",  "FGR2019_e1",  "FO2016_e2a",  "FO2016_e2b",  "FO2016_e3",   "GKH1999_e1",
                                                       "GKH1999_e2",  "GKH1999_e3",  "GKH1999_e4",  "HDM2006_e1",  "HDM2006_e2",  "HUW2015_e1",
                                                       "JCD2012_e1a", "JCD2012_e1b", "JCM2019_e1",  "JW2019_e1" ,  "JWH2009_e1",  "KAWY2013_e2",
                                                       "KAWY2013_e3", "KAWY2013_e4", "KFH2013_e1",  "KK2015_e1",   "KL2012_e1",   "KL2012_e2",
                                                       "KL2012_e3",   "KL2012_e4",   "KUO2017_e2",  "KY2010_e1",   "KY2011_e1",   "KY2016_e1",
                                                       "LBA2019_e1",  "LBA2019_e2",  "LM2020_e1",   "LP2006_e2",   "LP2013_e1",
                                                       "LP2013_e2","MKG2018_e2",
                                                       "NFR2013_e1",  "NFR2013_e2",  "NFR2013_e3",  "OBD2017_e1",  "OZH2010_e1",  "PCM2006_e1",
                                                       "PRM2010_e1",  "QGM2021_e1",  "RS2009_e1",   "RS2009_e2",   "RSM2009_e1",  "RSP2012_e1",
                                                       "SB2020_e1",   "SB2020_e2",   "SB2020_e3",   "SBT2018_e1",  "SCR2019_e1",  "SD2004_e2",
                                                       "SD2014_e1",   "SD2014_e2",   "SHJ2005_e1",  "TMP2014_e1",  "TR2017_e1",   "USS2015_e1",
                                                       "WHD2018_e1",  "WKH2020_e1",  "WKS2020_e1",  "WMS2012_e1",  "ZMD2011_e1",  "ZOL2021_e1" ),selected="AP2007_e1")
                          ),
                          # conditionalPanel(
                          #   condition = "input.radio == 1",
                          #   title = "Select experiment",
                          #   selectInput("exp1", "Select exp:",
                          #               choices = list("APH2016_e1",  "BKSR2013_e1", "BSG2014_e1" , "CSR2015_e1" , "D2007_e1a",   "D2007_e1b",   "DW2020_e1" ,
                          #                              "FBH2013_e1",  "FGR2019_e1" , "FO2016_e2a",  "FO2016_e2b" , "FO2016_e3" ,  "JW2019_e1" ,  "KAWY2013_e2",
                          #                              "KAWY2013_e3", "KAWY2013_e4", "KFH2013_e1",  "KK2015_e1",   "KUO2017_e2" , "KY2010_e1" ,  "KY2011_e1" ,
                          #                              "LBA2019_e1",  "LBA2019_e2" , "LM2020_e1" ,  "LP2013_e1",   "MKG2018_e2" , "OBD2017_e1" , "PCM2006_e1" ,
                          #                              "PRM2010_e1" , "RSP2012_e1" , "SB2020_e1" ,  "SB2020_e2",   "SB2020_e3",   "SBT2018_e1" , "SCR2019_e1",
                          #                              "SD2014_e1"  , "SD2014_e2"  , "SHJ2005_e1" , "TMP2014_e1",  "USS2015_e1",  "WHD2018_e1" , "WKS2020_e1" ,
                          #                              "WMS2012_e1"  ,"ZMD2011_e1" , "ZOL2021_e1"), selected = "APH2016_e1")
                          # ),
                          conditionalPanel(
                            condition = "input.radio == 3",
                            title = "Select experiment",
                            selectInput("exp3", "Select exp:",
                                        choices = list("JCD2012_e2","JKD2013_e1",

                                                       "SSW2015_e3","SWR2012_e1",
                                                       "SWR2012_e2","SWS2011_e1",
                                                       "WLH2011_e1","WLM2009_e1",
                                        "ZSB2021_e1"), selected = "ZSB2021_e1")
                          ),
                          fluidRow(
                            column(3,
                                   div(style = "padding:20px",

                                       downloadButton("downloadData2", "Download.zip"),

                                   )

                            )
                          )


                        ),
                        mainPanel(width=9,
                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",
                                               h2("Information"),
                                               htmlOutput("information")
                                           )
                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",
                                               h2("ROC"),
                                               div(style = "text-align:right",
                                                   a(id = "toggleROC", "Show explanation", href = "#")
                                               ),
                                               shinyjs::hidden(
                                                 div(id = "ROCinfo",
                                                     htmlOutput("ROCinfo")
                                                 )

                                               ),

                                               uiOutput("uilim")


                                           )



                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",
                                               h2("Download full data"),

                                               htmlOutput("fullinformation")
                                               )

                                    )
                                  ),

                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",

                                               downloadButton("downloadData", "Download",disabled = "disabled"),

                                           )

                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           div(style = "padding:20px",

                                               dataTableOutput("datpreview")
                                           )

                                    )
                                  ),
                        )
                      )
             )


  )
server <- function(input, output, session) {


  output$information <- renderUI({


     if (input$radio == 1){

      selectee <- input$exp1
    } else {
      selectee <- input$exp3
    }


    ids <- expinfo %>% filter(Experiment == selectee)

    str1 <- paste0("<b>Exp ID:</b> ",ids$Experiment)
    str2 <- paste0("<b>Paper</b>: ", ids$Citation)
    str3 <- paste0("<b>Item information:</b> ", ids$Iteminformation)
    str8 <- paste0("<b>Modification for fitting:</b> ", ids$MetaModification)
    str4 <- paste0("<b>Fitted experimental design:</b> ",ifelse(ids$MetaManipulation == "None","No manipulation",ids$MetaManipulation))
    str5 <- paste0("<b>Key manipulation:</b> ",ids$MetaKeymanipulation)
    str6 <- paste0("<b>Number of participants fitted:</b> ", ids$MetaParticipants, " / ", ids$TotParticipants)
    str7 <- paste0("<b>Stimuli:</b> ", ids$MetaStimuli)


    HTML(paste(str1, str2, str3, str8, str4, str5, str6, str7, sep = '<br/>'))

  })

  output$fullinformation <- renderUI({

    if (input$radio == 1){

      selectee <- input$exp1
    } else {

      selectee <- input$exp3
    }


    ids <- expinfo %>% filter(Experiment == selectee)

    str1 <- paste0("<b>Exp ID:</b> ",ids$Experiment)
    str2 <- paste0("<b>Paper</b>: ", ids$Citation)
    str3 <- paste0("<b>Number of participants:</b> ", ids$TotParticipants)
    str4 <- paste0("<b>Experimental design:</b> ",ifelse(ids$TotManipulation == "None","No manipulation",ids$TotManipulation))
    str5 <- paste0("<b>Factorial design:</b> ",ids$TotDesign)
    str6 <- paste0("<b>Trials per participant:</b> ", ids$Trialsperperson)
    str7 <- paste0("<b>Delay to test:</b> ", ids$Delaytotest)
    str8 <- paste0("<b>Stimuli:</b> ", ids$TotStimuli)
    str9 <- paste0("<b>Length confidence scale:</b> ", ids$RatingScaleSize)
    str10 <- paste0("<b>Original test format:</b> ", ids$Confidenceprocedure)
    str11 <- paste0("<b>Item information:</b> ",ids$Iteminformation)
    str12 <- paste0("<b>contains RT:</b> ", ids$ContainsRT)
    str13 <- paste0("<b>Trial order preserved:</b> ", ids$Trialorderpreserved)



    HTML(paste(str1, str2, str3, str4, str5, str6, str7,
               str8, str9, str10, str11, str12, str13,sep = '<br/>'))

  })


  shinyjs::onclick("toggleROC",
                   shinyjs::toggle(id = "ROCinfo", anim = TRUE))

  output$ROCinfo <- renderUI({


    str_e1 <- "ROC curves were calculated separately for all conditions (panels) and participants (lines). The coding of the conditions corresponds to the (dummy) coding of conditions during the model fit. Where condition varied between blocks and/or participants this is relatively self-explanatory. For within-subjects manipulations where two types of new items were manipulated within the same test block (i.e., condition != block), we used dummy coding and designated one set of new items as 'new' and the other set of new items as 'old'. Where this is the case, the panel heading carries the suffix '_0' (where '_1' stands for 'true' old items)."
    str_e2 <- "In case of multiple conditions, we collapsed conditions into a single condition with multiple levels. The panel heading shows this where two terms are connected by an underscore, i.e., 'HighF_Speed_0' indicating high frequency novel items that were responded to under speed-instructions. Each of these re-jigged conditions is shown in its own panel."
    str_e3 <- "Where there is only a single panel with a panel heading other than 'None', we fitted only one condition of a given experiment (with the experiment information suggestion that there was neither a between or within-subjects manipulation)."

    HTML(paste(str_e1, str_e2, sep = '<p/><p/>'))

  })


 output$uilim <- renderUI({

  output$rocplot1 <- renderPlot({


    if(input$radio == 1){
      exps <- exps1
      dexp <- d1
      selectee <- input$exp1
      inp<-which(explist==selectee)
      makeROC( dexp[[inp]])
    } else if (input$radio == 3) {

      NULL
    }



  })

  if(input$radio == 1) {

    exps <- exps1
    dexp <- d1
    selectee <- input$exp1

  inp<-which(exps==selectee)

  ln <-  length(unique(dexp[[inp]] %>% .$condition))

  nn <- ifelse((ln > 3 & ln <= 6) , 2, ifelse(ln > 6,3,1))

  plotOutput("rocplot1", height = 400 * nn)
  }

 })

 selectee <- reactive({
   if (input$radio == 1){

   sel <-  input$exp1

   } else {
     sel <- input$exp3
   }

   sel
 })

 output$datpreview <- renderDataTable(

  makeHead(selectee()),
   rownames= FALSE,
   options = list(searching = FALSE,
                  paging = TRUE,
                  pageLength = 10,
                  initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                    "}")
   )


 )


 observeEvent(selectee(), {
   if(!(selectee() %in% c("ZMD2011_e1","RSM2009_e1"))){
     enable("downloadData")
     runjs("$('#dwnbutton').removeAttr('title');")
   }else{
     disable("downloadData")
     runjs("$('#dwnbutton').attr('title', 'Data not available');")
   }
 })

 output$downloadData <-

   downloadHandler(
   filename = function() {
     paste(selectee(), ".csv", sep = "")
   },
   content = function(file) {
     write.csv(selectee(), file, row.names = FALSE)
   }
 )

 output$downloadData2 <-

   downloadHandler(
     filename = 'RecognitionData.zip',
     content = function(fname) {
       file.copy(from = "RecognitionData.zip",
                 to = fname)
     },
     contentType = "application/zip"

   )





  # shinyjs::onclick("toggleDeviance",
  #                  shinyjs::toggle(id = "devianceinfo", anim = TRUE))
  #
  #
  #
  # selectee <- reactive({
  #   if(input$radio == 2) {
  #
  #    sel <- input$exp2
  #
  #   } else {
  #
  #   sel <-  input$exp1
  #
  #   }
  #
  #   sel
  # })
  #
  # selecttype <- reactive({
  #   if(input$selectdev == 1) {
  #
  #     sel <- "Full"
  #
  #   } else if (input$selectdev == 2) {
  #
  #     sel <-  "KFCV"
  #
  #   } else {
  #
  #     sel <- "LOP"
  #
  #   }
  #
  #   sel
  # })
  #
  #
  #
  #
  # output$devianceinfo <- renderUI({
  #
  #
  #   str_e1 <- "Quantitative fit and out-of-sample prediction. The deviance is shown for all models, where the color gradient indicates the distribution of deviance across posterior samples. In the Full panel, this shows the in-sample deviance when the full data set was fitted. In the KFCV panel, this shows the out-of-sample prediction when ~10% of trials (randomly selected proportionally across participants and conditions) was held out from the training fit, summed across all K(= 10) folds. In the LOP panel, this shows the out-of-sample prediction when all trials of ~10% participants were held out from the training fit, summed across all K(= 10) folds. For the LOP out-of-sample prediction, participant-level parameters were sampled from the group-level distribution."
  #
  #
  #   HTML(paste(str_e1, sep = '<p/><p/>'))
  #
  # })
  #
  # output$devianceplot <- renderPlot({
  #
  #   makeDev(selectee())[[1]]
  #
  # })
  #
  # shinyjs::onclick("toggleDiv",
  #                  shinyjs::toggle(id = "divinfo", anim = TRUE))
  #
  #
  # output$divinfo <- renderUI({
  #
  #
  #   str_e1 <- "The table shows the percentage of divergent transitions during the model fit. For KFCV and LOP, this shows mean divergent transitions across training samples. We consider the fit problematic if this number grossly exceeds 1%."
  #
  #
  #   HTML(paste(str_e1, sep = '<p/><p/>'))
  #
  # })
  #
  # output$divplot <- renderDataTable(
  #
  #   makeDev(selectee())[[2]] %>% ungroup() %>% select(-experiment),
  #   #rownames= FALSE,
  #   options = list(searching = FALSE,
  #                  paging = FALSE,
  #
  #     initComplete = JS(
  #       "function(settings, json) {",
  #       "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
  #       "}")
  #   )
  #
  #
  # )
  #
  #
  # shinyjs::onclick("togglePpp",
  #                  shinyjs::toggle(id = "pppinfo", anim = TRUE))
  #
  #
  # output$pppinfo <- renderUI({
  #
  #
  #   str_e1 <- "The posterior predictive p-values (p) provide information on the adequacy of the model in capturing the observed data. Loosely following Klauer (2010), we judged models on their ability to capture the confidence ratings participants assigned to items. Given the different levels of the hierarchical model, the adequacy of the model can be determined at different levels. We calculated these p-values for the average participant, for all individuals, and for all trials for a given individual where item information was present in the dataset."
  #   str_e2 <- "Broadly, for all three cases statistics were calculated in the same way, with the aggregation of data at different levels. For each posterior sample, we determined the prediction that this sample makes of the data. We compared that prediction to the observed data (T1_obs), and to data simulated from that prediction (T1_pred), by calculating (ndat - nexp)^2/nexp. We then contrasted T1pred and T1obs and calculated for what proportion of samples for which T1pred < T1obs. When this proportion (p) is not small, the model is considered to adequately capture the data."
  #   str_e3 <- "For the average ppt p-value, we summed responses across all items for each participant to construct the frequency table for the number of items assigned a particular confidence response in a given item condition (i.e., new, old, etc.). We then averaged predicted, observed and expected frequency tables across participants and calculated p(T1pred < T1obs) for the averaged data across the posterior samples."
  #   str_e4 <- "For the individual-level p-value, we summed responses across all items for each participant. For each participant, we now calculated p(T1pred < T1obs). In the table we report the proportion of participants for whom this p-value > .05 across the posterior samples."
  #   str_e5 <- "For data sets with item information, we also calculated T1pred and T1obs on the item level. That is, we calculated the squared distance of observed and expected and predicted and expected response for each item, then summed these distances across items for each individual for a measure of T1pred and T1obs. We then calculated p(T1pred < T1obs) across the posterior samples for each participant, and report the proportion of participants for whom p > .05 in the table."
  #
  #
  #
  #   HTML(paste(str_e1,str_e2,str_e3, str_e4,str_e5, sep = '<p/><p/>'))
  #
  # })
  #
  # output$pppplot <- renderDataTable(
  #
  #  makeDev(selectee())[[3]]%>% ungroup() %>% select(-experiment),
  #   #rownames= FALSE,
  #   options = list(searching = FALSE,
  #                  paging = FALSE,
  #
  #                  initComplete = JS(
  #                    "function(settings, json) {",
  #                    "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
  #                    "}")
  #   )
  #
  #
  # )
  #
  #
  #
  #
  #   output$devsolo1 <- renderPlot({
  #
  #     summaryplots(input$selectdev,exps1)
  #
  #     })
  #
  #   output$devsolo2 <- renderPlot({
  #
  #     summaryplots(input$selectdev,exps2)
  #
  #   })
  #





}


shinyApp(ui, server)




