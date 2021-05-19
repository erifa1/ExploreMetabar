# Deploy a Prod-Ready, Robust Shiny Application.
#
# 4. Test my package

devtools::test()
rhub::check_for_cran()

# 5. Deployment elements

## 5.1 If you want to deploy on RStudio related platforms
golem::add_rstudioconnect_file()
golem::add_shinyappsio_file()
golem::add_shinyserver_file()

## 5.2 If you want to deploy via a generic Dockerfile
golem::add_dockerfile()

## 5.2 If you want to deploy to ShinyProxy
golem::add_dockerfile_shinyproxy()

## 5.2 If you want to deploy to Heroku
golem::add_dockerfile_heroku()


#Document update
devtools::document(roclets = c('rd', 'collate', 'namespace'))

# shinyappIO

options(repos = c(BiocManager::repositories()))

options()$repos


devtools::install_github("erifa1/ranomaly", force = TRUE)
devtools::install_github("mahendra-mariadassou/phyloseq-extended", force = TRUE, ref = "dev")
#github::mahendra-mariadassou/phyloseq-extended
## list all dependencies
rsconnect::appDependencies()

### installed package need to be same as deployed!!!
# phangorn github
# ranomaly github
#test appli
attachment::att_amend_desc()
options(repos = c(BiocManager::repositories()))
rsconnect::configureApp("test_xplrMeta", size="xlarge")
rsconnect::deployApp(appName="test_xplrMeta")

# deploy master
options(repos = c(BiocManager::repositories()))
rsconnect::configureApp("APPNAME", size="xlarge")
rsconnect::deployApp(appName="exploremetabar")
