require(plotly)
require(htmlwidgets)
require(reticulate)

np <- import ('numpy')
biomass <- np$load("github_repos/dingo/VBOF/samples_biomass_max.npy")
vbof <- np$load("github_repos/dingo/VBOF/samples_vbof_max.npy")

# Find reaction indexes of interest; for example when the differenc in the mean flux is greater than an order of magnitude
for (i in 1: 3394)
{
  x = mean(biomass[i,])
  y = mean(vbof[i,])
  if (x > 0 && y >0) 
  {
    min_val = min(x,y)
    dif = abs(x-y)
    if (dif > 10*min_val)
    {
      print(i)
    }    
  }
}


biomass_index <- 3393
vbof_index    <- 3394




dgk_1_index   <- 1200





## plot style ##
m = 8     # 5x5 copula
axx <- list(
  title = '% max flux 1',
  ticketmode = 'array',
  ticktext = c("0", "0.5", "1"),
  tickvals = c(0, m/2, m-1)
)
axy <- list(
  title = '% max flux 2',
  ticketmode = 'array',
  ticktext = c("0", "0.5", "1"),
  tickvals = c(0,m/2,m-1)
)
axz <- list(
  title = 'probability',
  ticketmode = 'array',
  ticktext = c(" "),
  tickvals = c(0)
)


# Biomass ~ VBOF for the maximizing biomass model
#-------------------------------------------------
rets = biomass[biomass_index,]
vols = biomass[vbof_index,]
max(rets)
max(vols)
cop = compute_copula(rets, vols, m)

title_char = "Maximizing for biomass, biomass ~ VBOF"   # you could remove this if you don't want any title
fig = plotly::plot_ly(z = ~cop) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title = title_char, scene = list(xaxis=axx, yaxis=axy, zaxis=axz))
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/VBOF/",title_char,".html"))

# Biomass ~ dgk_1 for the maximizing biomass model
#-------------------------------------------------
rets = biomass[biomass_index,]
vols = biomass[dgk_1_index,]
max(rets)
max(vols)
cop = compute_copula(rets, vols, m)

title_char = "Maximizing for biomass, biomass ~ GK1"   # you could remove this if you don't want any title
fig = plotly::plot_ly(z = ~cop) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title = title_char, scene = list(xaxis=axx, yaxis=axy, zaxis=axz))


# VBOF ~ dgk_1 for the maximizing biomass model
#-------------------------------------------------
rets = biomass[vbof_index,]
vols = biomass[dgk_1_index,]
cop = compute_copula(rets, vols, m)

title_char = "Maximizing for biomass, VBOF ~ GK1"   # you could remove this if you don't want any title
fig = plotly::plot_ly(z = ~cop) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title = title_char, scene = list(xaxis=axx, yaxis=axy, zaxis=axz))
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/VBOF/",title_char,".html"))

#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

# Biomass ~ VBOF for the maximizing vbof model
#-------------------------------------------------
rets = vbof[biomass_index,]
vols = vbof[vbof_index,]
cop = compute_copula(rets, vols, m)

title_char = "Maximizing for vbof, biomass ~ VBOF"   # you could remove this if you don't want any title
fig = plotly::plot_ly(z = ~cop) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title = title_char, scene = list(xaxis=axx, yaxis=axy, zaxis=axz))
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/VBOF/",title_char,".html"))


# Biomass ~ dgk_1 for the maximizing vbof model
#-------------------------------------------------
rets = vbof[biomass_index,]
vols = vbof[dgk_1_index,]
cop = compute_copula(rets, vols, m)

title_char = "Maximizing for vbof, biomass ~ GK1"   # you could remove this if you don't want any title
fig = plotly::plot_ly(z = ~cop) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title = title_char, scene = list(xaxis=axx, yaxis=axy, zaxis=axz))
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/VBOF/",title_char,".html"))


# VBOF ~ dgk_1 for the maximizing vbof model
#-------------------------------------------------
rets = vbof[vbof_index,]
vols = vbof[dgk_1_index,]
cop = compute_copula(rets, vols, m)

title_char = "Maximizing for vbof, VBOF ~ GK1"   # you could remove this if you don't want any title
fig = plotly::plot_ly(z = ~cop) %>% add_surface(showscale=FALSE)
fig = fig %>% layout(title = title_char, scene = list(xaxis=axx, yaxis=axy, zaxis=axz))
htmlwidgets::saveWidget(as_widget(fig), paste0("~/github_repos/dingo/VBOF/",title_char,".html"))





