# install.shiny <- !("shiny" %in% installed.packages()[,"Package"])
# if(install.shiny) install.packages("shiny")
library(shiny, warn.conflicts = FALSE)

# install.shinydashboard <- !("shinydashboard" %in% installed.packages()[,"Package"])
# if(install.shinydashboard) install.packages("shinydashboard")
library(shinydashboard, warn.conflicts = FALSE)

# install.shinyjs <- !("shinyjs" %in% installed.packages()[,"Package"])
# if(install.shinyjs) install.packages("shinyjs")
library(shinyjs, warn.conflicts = FALSE)

# install.DT <- !("DT" %in% installed.packages()[,"Package"])
# if(install.DT) install.packages("DT")
library(DT, warn.conflicts = FALSE)

# install.ggplot2 <- !("ggplot2" %in% installed.packages()[,"Package"])
# if(install.ggplot2) install.packages("ggplot2")
library(ggplot2, warn.conflicts = FALSE)

# install.plotly <- !("plotly" %in% installed.packages()[,"Package"])
# if(install.plotly) install.packages("plotly")
library(plotly, warn.conflicts = FALSE)

# install.shinythemes <- !("shinythemes" %in% installed.packages()[,"Package"])
# if(install.shinythemes) install.packages("shinythemes")
library(shinythemes, warn.conflicts = FALSE)

# install.shinycssloaders <- !("shinycssloaders" %in% installed.packages()[,"Package"])
# if(install.shinycssloaders) install.packages("shinycssloaders")
library(shinycssloaders, warn.conflicts = FALSE)


# install.shinyalert <- !("shinyalert" %in% installed.packages()[,"Package"])
# if(install.shinyalert) install.packages("shinyalert")
library(shinyalert, warn.conflicts = FALSE)


# install.shinyBS <- !("shinyBS" %in% installed.packages()[,"Package"])
# if(install.shinyBS) install.packages("shinyBS")
library(shinyBS, warn.conflicts = FALSE)


# install.V8 <- !("V8" %in% installed.packages()[,"Package"])
# if(install.V8) install.packages("V8")
library(V8, warn.conflicts = FALSE)


# install.tidyverse <- !("tidyverse" %in% installed.packages()[,"Package"])
# if(install.tidyverse) install.packages("tidyverse")
library(tidyverse, warn.conflicts = FALSE)

# install.tidyr <- !("tidyr" %in% installed.packages()[,"Package"])
# if(install.tidyr) install.packages("tidyr")
library(tidyr, warn.conflicts = FALSE)

# install.dplyr <- !("dplyr" %in% installed.packages()[,"Package"])
# if(install.dplyr) install.packages("dplyr")
library(dplyr, warn.conflicts = FALSE)

# install.matrixStats <- !("matrixStats" %in% installed.packages()[,"Package"])
# if(install.matrixStats) install.packages("matrixStats")
library(matrixStats, warn.conflicts = FALSE)


# install.shinyWidgets <- !("shinyWidgets" %in% installed.packages()[,"Package"])
# if(install.shinyWidgets) install.packages("shinyWidgets")
library(shinyWidgets, warn.conflicts = FALSE)

# install.devtools <- !("devtools" %in% installed.packages()[,"Package"])
# if(install.devtools) install.packages("devtools")
library(devtools, warn.conflicts = FALSE)

# install.cleaver <- !("cleaver" %in% installed.packages()[,"Package"])
# if(install.cleaver) devtools::install_github("sgibb/cleaver")
library(cleaver, warn.conflicts = FALSE)


# install.Biostrings <- !("Biostrings" %in% installed.packages()[,"Package"])
# if(install.Biostrings) install.packages("Biostrings")
library(Biostrings, warn.conflicts = FALSE)


# install.stringi <- !("stringi" %in% installed.packages()[,"Package"])
# if(install.stringi) install.packages("stringi")
library(stringi, warn.conflicts = FALSE)


# install.shinyjs <- !("shinyjs" %in% installed.packages()[,"Package"])
# if(install.shinyjs) install.packages("shinyjs")
library(shinyjs, warn.conflicts = FALSE)

# install.tictoc <- !("tictoc" %in% installed.packages()[,"Package"])
# if(install.tictoc) install.packages("tictoc")
library(tictoc, warn.conflicts = FALSE)

# install.data.table <- !("data.table" %in% installed.packages()[,"Package"])
# if(install.data.table) install.packages("data.table")
library(data.table, warn.conflicts = FALSE)


# install.shinybusy <- !("shinybusy" %in% installed.packages()[,"Package"])
# if(install.shinybusy) install.packages("shinybusy")
library(shinybusy, warn.conflicts = FALSE)


# install.httr <- !("httr" %in% installed.packages()[,"Package"])
# if(install.httr) install.packages("httr")
library(httr, warn.conflicts = FALSE)


## declaration of the maximum upload size ##
options(shiny.maxRequestSize = 200*1024^2)
options(shiny.sanitize.errors = TRUE) ##FALSE TO SHOW ERRORS, ALTERNATIVELY TRUE FOR DISTRIBUTION
options(dplyr.summarise.inform = FALSE)
# options(future.globals.maxSize= 891289600)

enzyme_list <- c("Arg-C proteinase"="arg-c proteinase",
                 "Asp-N Endopeptidase"="asp-n endopeptidase",
                 "BNPS-Skatole"="bnps-skatole", 
                 "Caspase 1"="caspase1",
                 "Caspase 2"="caspase2",
                 "Caspase 3"="caspase3",
                 "Caspase 4"="caspase4",
                 "Caspase 5"="caspase5",
                 "Caspase 6"="caspase6",
                 "Caspase 7"="caspase7",
                 "Caspase 8"="caspase8",
                 "Caspase 9"="caspase9",
                 "Caspase 10"="caspase10",
                 "Chymotrypsin-high specificity"="chymotrypsin-high",
                 "Chymotrypsin-low specificity"="chymotrypsin-low",
                 "Clostripain (Clostridiopeptidase B)"="clostripain",
                 "CNBr"="cnbr",
                 "Enterokinase"="enterokinase",
                 "Factor Xa"="factor xa", 
                 "Formic acid"="formic acid",
                 "Glutamyl endopeptidase"="glutamyl endopeptidase",
                 "Granzyme B"="granzyme-b", 
                 "Hydroxylamine"="hydroxylamine",
                 "Iodosobenzoic acid"="iodosobenzoic acid",
                 "LysC"="lysc",
                 "LysN"="lysn",
                 "Neutrophil elastase"="neutrophil elastase",
                 "NTCB (2-nitro-5-thiocyanobenzoic acid)"="ntcb", 
                 "Pepsin (pH == 1.3)"="pepsin1.3",
                 "Pepsin (pH > 2)"="pepsin",
                 "Proline endopeptidase"="proline endopeptidase",
                 "Proteinase k"="proteinase k",
                 "Staphylococcal Peptidase I"="staphylococcal peptidase i",  
                 "Thermolysin"="thermolysin",
                 "Thrombin"="thrombin",
                 "Trypsin"="trypsin")

# example_dataset <- HTML(paste0("P01023\nQ99758\nO15439"))
example_dataset <- HTML(paste0("P01023\nQ99758\nO15439\nO43184\nQ13444\nP82987\nP04083\nQ7Z5R6\nP27540\nQ13315\nP36543\nQ13535\nQ96GD4\nQ13145\nQ07812\nP56945\nP10415\nP11274\nP51587\nP55211\nP22681\nP15529\nP42771\nP49454\nQ99490\nP28329\nO14757\nO00533\nP02452\nP08123\nP02461\nP12110\nP49674\nP35222\nP09668\nQ9NQC7\nO15528\nP08684\nO43602\nQ16760\nP49619\nQ14185\nQ03001\nP00533\nQ09472\nQ96L91\nP21860\nQ969H0\nP11362\nP21802\nP17948\nP35916\nP02751\nO43524\nP98177\nP42345\nP09958\nP23771\nP63096\nP63092\nP42262\nQ9NRY4\nP07900\nQ13418\nQ9UNL4\nP06213\nP46940\nA1A580\nC9JFL3\nG5E9R7\nP60896\nQ9Y6X1"))


mol_weightAA <- c('A'=71.0788, 'C'=103.1388, 'D'=115.0886, 'E'=129.1155, 'F'=147.1766, 'G'=57.0519,
                  'H'=137.1411, 'I'=113.1594, 'K'=128.1741, 'L'=113.1594, 'M'=131.1926, 'N'=114.1038,
                  'O'=237.3018, 'P'=97.1167, 'Q'=128.1307, 'R'=156.1875, 'S'=87.0782, 'T'=101.1051,
                  'U'=150.0388, 'V'=99.1326, 'W'=186.2132, 'Y'=163.1760, 'H2O'=18.01524)

rowCallback <- c(
  "function(row, data){",
  "  for(var i=0; i<data.length; i++){",
  "    if(data[i] === null){",
  "      $('td:eq('+i+')', row).html('NA')",
  "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
  "    }",
  "  }",
  "}"
)