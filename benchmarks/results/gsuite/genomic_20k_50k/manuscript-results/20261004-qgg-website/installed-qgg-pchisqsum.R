function (x, df, a, lower.tail = TRUE) 
{
    sat <- satterthwaite(a, df)
    guess <- pchisq(x/sat$scale, sat$df, lower.tail = lower.tail)
    for (i in seq(length = length(x))) {
        lambda <- rep(a, df)
        sad <- sapply(x, saddle, lambda = lambda)
        if (lower.tail) 
            sad <- 1 - sad
        guess <- ifelse(is.na(sad), guess, sad)
    }
    return(guess)
}
