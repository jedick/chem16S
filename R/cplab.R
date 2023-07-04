# chem16S/cplab.R
# Function to format figure labels
# Adapted from canprot's cplab.R on 20230704

cplab <- function(varname) {
  switch(varname,
    nH2O = quote(italic(n)[H[2]*O]),
    DnH2O = quote(Delta*italic(n)[H[2]*O]),
    nO2 = quote(italic(n)[O[2]]),
    DnO2 = quote(Delta*italic(n)[O[2]]),
    Zc = quote(italic(Z)[C]),
    DZc = quote(Delta*italic(Z)[C]),
    logfO2 = quote(log~italic("f")[O[2]]),
    logaH2O = quote(log~italic("a")[H[2]*O]),
    nC = quote(italic(n)[C] * "/AA"),
    nN = quote(italic(n)[N] * "/AA"),
    nS = quote(italic(n)[S] * "/AA"),
    DnC = quote(Delta*italic(n)[C] * "/AA"),
    DnN = quote(Delta*italic(n)[N] * "/AA"),
    DnS = quote(Delta*italic(n)[S] * "/AA"),
    V0 = quote(list(italic("V") * degree, "cm" ^ 3 ~ "mol" ^ -1)),
    DV0 = quote(list(Delta * italic("V") * degree, "cm" ^ 3 ~ "mol" ^ -1)),
    nAA = quote(italic(n)[AA]),
    DnAA = quote(Delta*italic(n)[AA]),
    GRAVY = "GRAVY",
    DGRAVY = quote(Delta*"GRAVY"),
    pI = "pI",
    DpI = quote(Delta*"pI"),
    MW = quote("MW"),
    DMW = quote(Delta * "MW"),
    varname
  )
}

