library(vegan)
data(dune)
data(dune.env)
adonis(dune ~ Management*A1, data=dune.env, permutations=99)


### Example of use with strata, for nested (e.g., block) designs.

dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y)
plot(mod)
### Hulls show treatment
with(dat, ordihull(mod, group=NO3, show="0"))
with(dat, ordihull(mod, group=NO3, show="10", col=3))
### Spider shows fields
with(dat, ordispider(mod, group=field, lty=3, col="red"))

### Correct hypothesis test (with strata)
adonis(Y ~ NO3, data=dat, strata=dat$field, perm=999)

### Incorrect (no strata)
adonis(Y ~ NO3, data=dat, perm=999)

load(system.file("data_test", "robjects_600.Rdata", package="ExploreMetabar"))
otable = otu_table(data)
mdata = data.frame(sample_data(data))
mdata$Depth <- sample_sums(data)

dd <<- vegdist(t(otable), distance='bray')
dd
form <- 'dd ~ SampleType*Primer'
res.adonis = adonis(as.formula(form), data = mdata, permutations = 1000)
res.adonis
