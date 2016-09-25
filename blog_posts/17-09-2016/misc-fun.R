##======================================================================
## Auxilar function for GLM models por Cesar Taconeli
##======================================================================

envelope <- function(modelo) {
    dados <- na.omit(modelo$data)
    nsim <- 100
    n <- modelo$df.null + 1
    r1 <- sort(rstandard(modelo, type = "deviance"))
    m1 <- matrix(0, nrow = n, ncol = nsim)
    a2 <- simulate(modelo, nsim = nsim)
    for (i in 1:nsim) {
        dados$y <- a2[, i]
        aj <- update(modelo, y ~ ., data = dados)
        m1[, i] <- sort(rstandard(aj, type = "deviance"))
    }
    li <- apply(m1, 1, quantile, 0.025)
    m <- apply(m1, 1, quantile, 0.5)
    ls <- apply(m1, 1, quantile, 0.975)
    quantis <- qnorm((1:n - 0.5)/n)
    plot(rep(quantis, 2), c(li, ls), type = "n",
         xlab = "Percentis teóricos",
         ylab = "Resíduos")
    ## title('Gráfico Normal de Probabilidades')
    lines(quantis, li, type = "l")
    lines(quantis, m, type = "l", lty = 2)
    lines(quantis, ls, type = "l")
    points(quantis, r1, pch = 16, cex = 0.75)
}

CHosmer <- function(modelo, g) {
    respostas <- modelo$y
    ## Dataframe com probabilidades estimadas e respostas para cada
    ## indivíduo.
    preditos <- predict(modelo, type = "response")
    ## Ordenando as linhas do dataframe da menor para a maior
    ## probabilidade estimada.
    dpred <- data.frame(preditos, respostas)
    ## Calculando os quantis para as probabilidades estimadas, para
    ## posterior formação dos grupos.
    dpred <- dpred[order(dpred[, 1]), ]
    ## Formando g grupos, de tamanhos (aproximadamente) iguais, com
    ## probabilidades estimadas semelhantes.
    cortes <- quantile(dpred[, 1], probs = seq(0, 1, 1/g),
                       include.lowest = TRUE)
    c1 <- cut(dpred[, 1], breaks = cortes, include.lowest = T)
    ## Obs é um vetor com o número observado de respostas em cada um dos
    ## g grupos.
    Obs <- tapply(dpred[, 2], c1, sum)
    ## pi é um vetor com as médias das probabilidades estimadas em cada
    ## um dos g grupos.
    pi <- tapply(dpred[, 1], c1, mean)
    ## n é um vetor com os tamanhos de amostras em cada grupo.
    n <- tapply(dpred[, 1], c1, length)
    ## Estatística do teste de qualidade proposto por Hosmer e Lemeshow e
    Cchap <- sum(((Obs - n * pi)^2)/(n * pi * (1 - pi)))
    ## p-valor correspondente.
    pv <- 1 - pchisq(Cchap, g - 2)
    return(list(Cchap = Cchap, pvalue = pv))
}

diag.binom <- function(fit.model){
    X <- model.matrix(fit.model)
    n <- nrow(X)
    p <- ncol(X)
    w <- fit.model$weights
    W <- diag(w)
    H <- solve(t(X)%*%W%*%X)
    H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
    h <- diag(H)
    ts <- resid(fit.model,type="pearson")/sqrt(1-h)
    td <- resid(fit.model,type="deviance")/sqrt(1-h)
    di <- (h/(1-h))*(ts^2)
    a <- max(td)
    b <- min(td)
    par(mfrow=c(1,3))
    plot(td,xlab="Índice", ylab="Resíduo Componente do Desvio",
         ylim=c(b-1,a+1), pch=20)
    title(sub="(a)")
    abline(h=c(2, -2),lty=2)
    plot(fitted(fit.model), td,xlab="Valor Ajustado",
         ylab="Resíduo Componente do Desvio", pch=20)
    title(sub="(b)")
    abline(h=c(-2, 0, 2), lty=2)
    envelope(fit.model)
    title(sub="(c)")
}

diag.pois <- function(fit.model, da, ni){
    X <- model.matrix(fit.model)
    n <- nrow(X)
    p <- ncol(X)
    w <- fit.model$weights
    W <- diag(w)
    H <- solve(t(X)%*%W%*%X)
    H <- sqrt(W)%*%X%*%H%*%t(X)%*%sqrt(W)
    h <- diag(H)
    ts <- resid(fit.model,type="pearson")/sqrt(1-h)
    td <- resid(fit.model,type="deviance")/sqrt(1-h)
    di <- (h/(1-h))*(ts^2)
    par(mfrow=c(3,2))
    a <- min(td)
    b <- max(td)
    plot(fitted(fit.model), h,xlab="Valor Ajustado", ylab="Medida h",
         pch=16)
    ## identify(fitted(fit.model), h, n=1)
    ##
    envelope(m3o)
    ## identify(di, n=1)
    ##
    plot(td,xlab="Índice", ylab="Resíduo Componente do Desvio",
         ylim=c(a-1,b+1), pch=16)
    abline(2,0,lty=2)
    abline(-2,0,lty=2)
    ## identify(td, n=1)
    ##
    w <- fit.model$weights
    eta <- predict(fit.model)
    z <- eta + resid(fit.model, type="pearson")/sqrt(w)
    plot(predict(fit.model),z,xlab="Preditor Linear",
         ylab="Variavel z", pch=16)
    lines(smooth.spline(predict(fit.model), z, df=2))
    hatv <- ((influence.measures(m3o)$infmat[,15]))
    dicook <- ((influence.measures(m3o)$infmat[,14]))
    plot(hatv, ylab="Alavancagem", type="n")
    grid()
    lines(x=1:nrow(da), y=hatv, type="h", ylab="Alavancagem")
    points(x=1:nrow(da), y=hatv, pch=20)
    text(x = 45 - ni, hatv[45], "45", cex=0.95)
    ##
    plot(dicook, ylab="Dist. de Cook", type="n")
    grid()
    lines(x=1:nrow(da), y=dicook, type="h")
    points(x=1:nrow(da), y=dicook, pch=20)
    text(x = 69 - ni, dicook[69], "69", cex=0.95)
    text(x = 45 - ni, dicook[45], "45", cex=0.95)
    par(mfrow=c(1,1))
}

##======================================================================
## Lattice functions por Walmes Zeviane
##  Funções desenvolvidas para matérias do blog 'Ridículas <-
##  dicas curtas sobre R' www.ridiculas.wordpress.com, por
##  Walmes Zeviani envie sugestões para walmes<at>ufpr.br
##======================================================================

library(lattice)
library(latticeExtra)

##----------------------------------------------------------------------
##  painel para gráficos com contornos
panel.3d.contour <- function(x, y, z, rot.mat, distance, type = "on",
    nlevels = 20, zlim.scaled, col.contour = 1, ...) {
    clines <- contourLines(x, y, matrix(z, nrow = length(x),
        byrow = TRUE), nlevels = nlevels)
    if (any(type %in% c("bottom"))) {
        for (ll in clines) {
            n <- ltransform3dto3d(rbind(ll$x, ll$y, zlim.scaled[1]),
                rot.mat, distance)
            panel.lines(n[1, ], n[2, ], col = col.contour, lty = 1,
                lwd = 1)
        }
    }
    panel.3dwire(x, y, z, rot.mat, distance, zlim.scaled = zlim.scaled,
        ...)
    if (any(type %in% c("on"))) {
        for (ll in clines) {
            n <- ltransform3dto3d(rbind(ll$x, ll$y, ll$level),
                rot.mat, distance)
            panel.lines(n[1, ], n[2, ], col = col.contour, lty = 1,
                lwd = 1)
        }
    }
    if (any(type %in% c("top"))) {
        for (ll in clines) {
            n <- ltransform3dto3d(rbind(ll$x, ll$y, zlim.scaled[2]),
                rot.mat, distance)
            panel.lines(n[1, ], n[2, ], col = col.contour, lty = 1,
                lwd = 1)
        }
    }
}

## Exemplo de uso
da <- expand.grid(x = seq(0, 1, l = 20), y = seq(0, 1, l = 20))
da$z <- with(da, 10 + x + 0.5 * y + 2 * x * y - 1 * x^2 - 1 * y^2)
wireframe(z ~ x + y, da, zlim = c(8, 12), drape = TRUE, col = "white",
    col.contour = "red", panel.3d.wireframe = "panel.3d.contour",
    type = c("on", "bottom"), screen = list(z = 40, x = -75))

##----------------------------------------------------------------------
##  painéis para fazer curvas de regressão com bandas de
##  confiança

prepanel.ciH <- function(x, y, ly, uy, subscripts, ...) {
    x <- as.numeric(x)
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    list(ylim = range(uy, ly, finite = TRUE), xlim = range(x))
}

panel.ciH <- function(x, y, ly, uy, subscripts, ...) {
    y <- as.numeric(y)
    x <- as.numeric(x)
    or <- order(x)
    ly <- as.numeric(ly[subscripts])
    uy <- as.numeric(uy[subscripts])
    panel.polygon(c(x[or], x[rev(or)]), c(ly[or], uy[rev(or)]),
        col = 1, alpha = 0.15, border = NA)
    panel.lines(x[or], ly[or], lty = 3, lwd = 0.5, col = 1)
    panel.lines(x[or], uy[or], lty = 3, lwd = 0.5, col = 1)
    panel.xyplot(x, y, subscripts = subscripts, ...)
}

## ## Exemplo de uso
## da <- expand.grid(A = gl(2, 4), B = 0:10)
## da$y <- with(da, rnorm(A, mean = as.numeric(A) * B, sd = 2))
## m0 <- lm(y ~ A * B, da)
## da <- cbind(da, predict(m0, interval = "confidence"))
## str(da)
##
## ## Sem grupos
## xyplot(y ~ B | A, da) +
##     as.layer(
##         xyplot(fit ~ B | A,
##                da,
##                type = "l",
##                ly = da$lwr,
##                uy = da$upr,
##                prepanel = prepanel.ciH,
##                panel = panel.ciH))
##
## ## Com grupos
## xyplot(y ~ B, groups = A, da) +
##     as.layer(
##         xyplot(fit ~ B,
##                groups = A,
##                da,
##                type = "l",
##                ly = da$lwr,
##                uy = da$upr,
##                prepanel = prepanel.ciH,
##                panel = panel.superpose,
##                panel.groups = panel.ciH))
