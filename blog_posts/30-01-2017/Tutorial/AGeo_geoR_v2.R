# Rotina para análise de geoestatística utilizando o pacote geoR ----------
# Rodolfo Marcondes Silva Souza - rodolfomssouza@gmail.com


# Configurações globais ---------------------------------------------------
#options(OutDec='.')


# Carregando pacotes ------------------------------------------------------
require('geoR') 
require('Hmisc')
require('calibrate')
require('fields')
require('fBasics')


# Alterando diretório de trabalho -----------------------------------------
setwd(dir = '~/')
setwd('Programming/R/Geoestatística/')


# Carregando os dados -----------------------------------------------------
dgeo=read.table('Dados_geo.txt', h=T)
attach(dgeo); names(dgeo)

nome_analise = 'Bacia_Pajeu_Erosividade'


# Carregando os dados no módulo de geoestatística (geoR) ------------------
dGeo=read.geodata('Dados_geo.txt', h=T, coords.col=1:2, data.col=3)


#  Análise estatística ----------------------------------------------------
estb=basicStats(Z, ci=0.95); estb
tks=ks.test(Z, 'pnorm', mean=mean(Z), sd=sd(Z)); tks
shapiro.test(Z)
jarqueberaTest(Z)
boxplot(Z)


# Iniciando geoestatística ------------------------------------------------
# Máxima distância
mx=max(X)
my=max(Y)
mdta=sqrt(mx^2+my^2);mdta; mdta/3


# Gerando e plotando o semivariograma -------------------------------------
plot(variog(dGeo))
mdt=130; nlg=6 # máxima distância e número de logs
svgteorico=variog(dGeo, max.dist=mdt, uvec=nlg)
svgteorico1=variog(dGeo, max.dist=mdt, uvec=nlg, direction=pi/8)
svgteorico2=variog(dGeo, max.dist=mdt, uvec=nlg, direction=pi/4)
svgteorico3=variog(dGeo, max.dist=mdt, uvec=nlg, direction=pi/2)
par(mfrow=c(2,2))
plot(svgteorico); plot(svgteorico1, main=expression(pi/8)); plot(svgteorico2, main=expression(pi/4)); plot(svgteorico3, main=expression(pi/2))
layout(1)
h=svgteorico$u
v=svgteorico$v
npar=svgteorico$n
ltmx=(max(h)+0.4*max(h))
ltmy=(max(v)+0.4*max(v))
ic=(1.96*(sqrt(2*v)/sqrt(npar)))
is=abs(v+ic)
ii=abs(v-ic)
plot(svgteorico, las=1, xaxs='i', yaxs='i',pch=16, col='red', ylim=c(0,ltmy), xlim=c(0,ltmx))
textxy(h,v,npar, cex=0.7)
arrows(x0=h,x1=h,y0=ii,y1=is, code=3, angle=90, length=0.05)
sm=data.frame(h,v,npar, row.names=NULL);sm


# Fazendo ajuste do semivariograma teórico --------------------------------
esferico=variofit(svgteorico, cov.model='sph', max.dist=max(h), messages=F)
exponencial=variofit(svgteorico, cov.model='exponential', max.dist=max(h), messages=F)
gaussiano=variofit(svgteorico, cov.model='gaussian', max.dist=max(h),messages=F)
x11(); sentimento=eyefit(svgteorico, silent=F)
stp = unlist(sentimento)
sigmasq=as.numeric(stp[2]); phi=as.numeric(stp[3]); tausq=as.numeric(stp[4]);


# Plotando todos o semivariograma e seus ajustes --------------------------
par(mfrow=c(2,2),mar=c(5,5,2,2))
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Esférico')
lines.variomodel(esferico, col='red' , lwd=2, lty=1)
arrows(x0=h,x1=h,y0=ii,y1=is, code=3, angle=90, length=0.05)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Exponencial')
lines.variomodel(exponencial, col='red' , lwd=2, lty=1)
arrows(x0=h,x1=h,y0=ii,y1=is, code=3, angle=90, length=0.05)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Gaussiano')
lines.variomodel(gaussiano, col='red' , lwd=2, lty=1)
arrows(x0=h,x1=h,y0=ii,y1=is, code=3, angle=90, length=0.05)
plot(svgteorico, las=1, type='p',pch=19, cex=1.4, col='black' , xlab='', ylim=c(0,ltmy), xlim=c(0,ltmx),ylab='', xaxs='i', yaxs='i', cex.lab=1.4, cex.axis=1.2, main='Sentimento')
lines(sentimento, col='red' , lwd=2, lty=1)
arrows(x0=h,x1=h,y0=ii,y1=is, code=3, angle=90, length=0.05)

dev.copy(pdf, paste(nome_analise, '_Semivariogramas.pdf', sep = ''),  width=8, height=6)
dev.off()
layout(1)


# Parâmetros ajustados para os modelos ------------------------------------ 
c0.esf=esferico$nugget; c.esf=esferico$cov.pars[1]; a.esf=esferico$cov.pars[2]
c0.exp=exponencial$nugget; c.exp=exponencial$cov.pars[1]; a.exp=exponencial$cov.pars[2]
c0.gau=gaussiano$nugget; c.gau=gaussiano$cov.pars[1]; a.gau=gaussiano$cov.pars[2]
C0.snt = sigmasq; C.snt = tausq; a.snt = phi


# Validação cruzada dos ajustes -------------------------------------------
cv.esf=xvalid(dGeo, model=esferico); zsco.esf=cv.esf$std.error
cv.exp=xvalid(dGeo, model=exponencial); zsco.exp=cv.exp$std.error
cv.gau=xvalid(dGeo, model=gaussiano); zsco.gau=cv.gau$std.error
cv.sent=xvalid(dGeo, model=sentimento); zsco.sent=cv.sent$std.error


# Média e variância do erro reduzido --------------------------------------
jkmed.esf=round(mean(zsco.esf),5); jkvar.esf=round(var(zsco.esf),5)
jkmed.exp=round(mean(zsco.exp),5); jkvar.exp=round(var(zsco.exp),5)
jkmed.gau=round(mean(zsco.gau),5); jkvar.gau=round(var(zsco.gau),5)
jkmed.sent=round(mean(zsco.sent),5); jkvar.sent=round(var(zsco.sent),5)


# Resumo das análises para escolha do melhor modelo -----------------------
modelos=c('Esf', 'Exp','Gau','Sent')
m.jk=rbind(jkmed.esf,jkmed.exp, jkmed.gau, jkmed.sent)
v.jk=rbind(jkvar.esf, jkvar.exp, jkvar.gau, jkvar.sent)
c0.smfit=rbind(c0.esf, c0.exp, c0.gau, tausq)
c.smfit=rbind(c.esf, c.exp, c.gau, sigmasq)
a.smfit=rbind(a.esf, a.exp, a.gau, phi)

resumo=data.frame(row.names=modelos, c0.smfit, c.smfit, a.smfit, m.jk, v.jk)
resumo
write.table(x = resumo, file = paste(nome_analise, '_Resumo_Geo.txt', sep = ''))


# Escolha do melhor modelo ------------------------------------------------
# esferico; exponencial; gaussiano; sentimento
smfit = esferico


# Gerando o grid para interpolações ---------------------------------------
ndvi = 200 # Tamanho do intervalo para interpolação
x.range <- as.integer(range(X))
y.range <- as.integer(range(Y))
grid.map=expand.grid(x=seq(from=x.range[1], to=x.range[2], by= (x.range[2] - x.range[1])/ndvi),
                     y=seq(from=y.range[1], to=y.range[2], by=(y.range[2] - y.range[1])/ndvi))
plot(grid.map)


# Carregando limites da borda ---------------------------------------------
lmt=read.table('Dados_Contorno.txt', h=T)
dlmt=read.geodata('Dados_Contorno.txt', h=T, coords.col=1:2, data.col=NULL)


# Fazendo Krigagem --------------------------------------------------------
#krg=krige.conv(dGeo, locations=grid.map, krige=krige.control(obj.model=smfit))
krg=krige.conv(dGeo, locations=grid.map, krige=krige.control(obj.model=smfit), borders=dlmt)


# Definindo parâmetros para os mapas --------------------------------------
numlevel=100
escala=matrix(krg$predict, ncol=1)


# Gráfico de contornos ----------------------------------------------------
contour(krg, f=T, col=terrain.colors(12), nlevels=10)

# Configurações para salvar o semivariograma ------------------------------
#png(filename = paste(nome_analise, '_smg.png', sep = ''), units='cm', res=600, width=12, height=10, pointsize=12)
par(family=('Times'), las=1, xaxs='i', yaxs='i', mar=c(5,5,2,2), cex.axis = 1.0, cex.lab = 1.2)
plot(svgteorico, type='p',pch=19,col='black' , ylim=c(0,(max(sm$v)+max(sm$v)*0.15)), xlim=c(0,(max(sm$h)+max(sm$h)*0.15)), ylab='', xlab='h (m)')#, yaxt='n')
vr=var(dgeo$Z);
abline(v=NULL, h=vr, lty=2, lwd=1, untf=3)
lines(smfit, col='black' , lwd=2, lty=1)
mtext(text=expression(gamma(h)), side=2, line=3, las=3, cex=1.2)
dev.copy(pdf, paste(nome_analise, '_smg.pdf', sep = ''), width=8, height=6)
dev.off()


# Configurações para o mapa -----------------------------------------------
# Redefinindo número da escala
numlevel=10
# Sequências para a legenda
sl=seq(min(krg$predict), max(krg$predict),by=(max(krg$predict)-min(krg$predict))/numlevel)
sll=formatC(sl, digits=2, format='f', decimal.mark = ".")

# Mapa
par(family=('Times'), las=1, xaxs='i', yaxs='i', mar=c(5,5,2,6.5), cex.axis = 1.0, cex.lab = 1.2)
plot(dgeo$X,dgeo$Y, type='n', xlab='X (UTM) x 1000', ylab='Y (UTM) x 1000', yaxt = 'n',
     xlim = c(min(dgeo$X)*0.95, max(dgeo$X)*1.05), ylim = c(min(dgeo$Y)*0.999, max(dgeo$Y)*1.001))
axis(side = 2, at = seq(min(Y), max(Y), (max(Y)-min(Y))/5), las = 3)
image(krg,  add=T, col=gray(seq(1,0,l=numlevel)), zlim=c(min(escala),max(escala)))
contour(krg, add=T,  levels=sl, labcex=0.5, lwd=1, labels=sll, drawlabels=F)
lines(lmt, lwd=2)
box()
points(dgeo$X,dgeo$Y, pch='+', col='red')
image.plot(escala, nlevel = numlevel, zlim=c(min(escala), max(escala)), col=gray(seq(1,0,l=numlevel)), legend.only = T, breaks=sl, lab.breaks=sll)
#image.plot(escala, nlevel = numlevel, zlim=c(min(escala), max(escala)), col=rev(terrain.colors(numlevel)), legend.only = T, breaks=sl, lab.breaks=sll)#, legend.lab=expression('CAD,'~g~g^{-1}), legend.line=3, legend.mar=6)
dev.copy(pdf, paste(nome_analise, '_mapa.pdf', sep = ''), width=8, height=8)
dev.off()

# Mapa do desvio padrão
esc.var=matrix(krg$krige.var, ncol=1)
esc.desv=matrix(sqrt(krg$krige.var), ncol=1)
#numlevel=4

# Sequências para a legenda
slsd=seq(min(esc.desv), max(esc.desv),by=(max(esc.desv)-min(esc.desv))/numlevel)
sllsd=formatC(slsd, digits=2, format='f', decimal.mark = ".")

# Mapa
par(family=('Times'), las=1, xaxs='i', yaxs='i', mar=c(5,5,2,6.5), cex.axis = 1.0, cex.lab = 1.2)
plot(dgeo$X,dgeo$Y, type='n', xlab='X (UTM) x 1000', ylab='Y (UTM) x 1000', yaxt = 'n',
       xlim = c(min(dgeo$X)*0.95, max(dgeo$X)*1.05), ylim = c(min(dgeo$Y)*0.999, max(dgeo$Y)*1.001))
axis(side = 2, at = seq(min(Y), max(Y), (max(Y)-min(Y))/5), las = 3)
image(krg, val=(krg$krige.var)^0.5, add=T, col=gray(seq(1,0,l=numlevel)), , zlim=c(min(esc.desv),max(esc.desv)))
lines(lmt, lwd=2)
box()
points(dgeo$X,dgeo$Y, pch='+', col='red')
image.plot(esc.desv, nlevel = numlevel, zlim=c(min(esc.desv), max(esc.desv)), col=gray(seq(1,0,l=numlevel)), legend.only = T, breaks=slsd, lab.breaks=sllsd)
dev.copy(pdf, paste(nome_analise, '_mapa_desvio_padrao.pdf', sep = ''), width=8, height=8)
dev.off()


# Salvando sessão do R ----------------------------------------------------
save.image(paste(nome_analise, '_.RData', sep = ''))

# Fim
