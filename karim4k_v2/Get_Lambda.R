Get_Lambda <-
function (K)
{
    out.s <- eigen(K, symmetric = TRUE)
    lambda1 <- out.s$values
    IDX1 <- which(lambda1 >= 0)
    IDX2 <- which(lambda1 > mean(lambda1[IDX1])/1e+05)
    if (length(IDX2) == 0) {
        stop("No Eigenvalue is bigger than 0!!")
    }
    lambda <- lambda1[IDX2]
    lambda
}

