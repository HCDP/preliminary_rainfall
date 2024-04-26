library(httr)

notify <- function(recepients, source, type, message) {
    recepients_cat = paste(recepients, collapse = "\",\"")
    notification <- paste0('{"recepients": ["', recepients_cat, '"],"source": "', source, '","type": "', type, '","message": "', message, '"}')
    res <- POST("https://api.hcdp.ikewai.org//notify", content_type_json(), body = notification, add_headers(`Authorization` = paste0("Bearer ", Sys.getenv("NOTIFY_TOKEN"))))
}