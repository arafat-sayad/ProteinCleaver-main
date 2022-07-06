customMenuSubItem <- function (text, tabName = NULL, href = NULL, newtab = TRUE, 
                          icon = shiny::icon("angle-double-right"), 
                          selected = NULL, badgeLabel = NULL, badgeColor = "green") {
    if (!is.null(href) && !is.null(tabName)) {
        stop("Can't specify both href and tabName")
    }
    isTabItem <- FALSE
    target <- NULL
    if (!is.null(badgeLabel)) {
        badgeTag <- tags$small(class = paste0("badge pull-right bg-", 
            badgeColor), badgeLabel)
    }
    else {
        badgeTag <- NULL
    }

    if (!is.null(tabName)) {
        shinydashboard:::validateTabName(tabName)
        isTabItem <- TRUE
        href <- paste0("#shiny-tab-", tabName)
    }
    else if (is.null(href)) {
        href <- "#"
    }
    else {
        if (newtab) 
            target <- "_blank"
    }
    tags$li(a(href = href, `data-toggle` = if (isTabItem) 
        "tab", `data-value` = if (!is.null(tabName)) 
        tabName, `data-start-selected` = if (isTRUE(selected)) 
        1
    else NULL, target = target, icon, span(text), badgeTag))
}