
# ================================================================
# MODELO COMPLETO 2001–2026 — 4 variables
# Hoja "Datos" (mensual, 4 variables)
# MCO en primeras diferencias + diagnósticos + robustez
# ================================================================

# ---------------- PASO 1: PAQUETES ----------------
library(readxl)
library(dplyr)
library(zoo)
library(lmtest)
library(sandwich)
library(car)
library(tseries)
library(strucchange)
library(urca)
library(lubridate)
library(ggplot2)

# ---------------- PASO 2: CARGA ----------------
ruta <- "RUTA DE LA BASE DE DATOS: BASE DE DATOS DEL MODELO.xlsm"
datos_raw <- read_excel(ruta, sheet = "Data")

# ---------------- PASO 3: LIMPIEZA ----------------
datos <- datos_raw %>%
  rename(
    PERIODO    = 1,
    CARTERA    = 2,
    FBKF       = 3,
    CONCENT    = 4,
    MARGEN_PAR = 5
  ) %>%
  mutate(
    PERIODO = as.Date(PERIODO),
    across(-PERIODO, ~ as.numeric(.))
  ) %>%
  filter(!is.na(PERIODO)) %>%
  distinct(PERIODO, .keep_all = TRUE) %>%       # elimina duplicado 2026-08 con FORECAST
  arrange(PERIODO) %>%
  filter(PERIODO >= as.Date("2001-01-01")) %>%  # inicio del análisis
  filter(!is.na(FBKF), !is.na(CARTERA),
         !is.na(CONCENT), !is.na(MARGEN_PAR))

cat("\n=== RANGO FINAL ===\n")
cat("Desde:", format(min(datos$PERIODO)), "\n")
cat("Hasta:", format(max(datos$PERIODO)), "\n")
cat("N observaciones:", nrow(datos), "\n")

# ---------------- PASO 4: LOGS Y DIFERENCIAS ----------------
datos <- datos %>%
  mutate(
    ln_FBKF       = log(FBKF),
    ln_CARTERA    = log(CARTERA),
    ln_CONCENT    = log(CONCENT),
    ln_MARGEN_PAR = log(MARGEN_PAR),
    d_ln_FBKF       = c(NA, diff(ln_FBKF)),
    d_ln_CARTERA    = c(NA, diff(ln_CARTERA)),
    d_ln_CONCENT    = c(NA, diff(ln_CONCENT)),
    d_ln_MARGEN_PAR = c(NA, diff(ln_MARGEN_PAR))
  )

# ---------------- PASO 5: MODELO PRINCIPAL (2001–2026) ----------------
modelo <- lm(d_ln_FBKF ~ d_ln_CARTERA + d_ln_CONCENT + d_ln_MARGEN_PAR,
             data = datos)

cat("\n=== MODELO PRINCIPAL: MCO EN DIFERENCIAS 2001–2026 ===\n")
print(summary(modelo))

# ---------------- PASO 6: DIAGNÓSTICOS ----------------
cat("\n=== DIAGNÓSTICOS ===\n")
cat("N obs:", nobs(modelo), "\n")
cat("R2 ajustado:", round(summary(modelo)$adj.r.squared, 4), "\n")

cat("\nBreusch-Pagan (heteroscedasticidad):\n");  print(bptest(modelo))
cat("\nDurbin-Watson (autocorrelación):\n");       print(dwtest(modelo))
cat("\nBreusch-Godfrey (4 rezagos):\n");           print(bgtest(modelo, order = 4))
cat("\nShapiro-Wilk (normalidad):\n");             print(shapiro.test(residuals(modelo)))
cat("\nVIF (multicolinealidad):\n");               print(vif(modelo))

# ---------------- PASO 7: NEWEY-WEST ----------------
cat("\n=== NEWEY-WEST (errores robustos) ===\n")
print(coeftest(modelo, vcov = NeweyWest(modelo, lag = 4, prewhite = FALSE)))

# ---------------- PASO 8: ADF (orden de integración) ----------------
cat("\n=== ADF NIVELES ===\n")
for (v in c("ln_FBKF","ln_CARTERA","ln_CONCENT","ln_MARGEN_PAR")) {
  cat("\n---", v, "---\n")
  print(ur.df(datos[[v]], type = "drift", selectlags = "AIC"))
}
cat("\n=== ADF PRIMERAS DIFERENCIAS ===\n")
for (v in c("d_ln_FBKF","d_ln_CARTERA","d_ln_CONCENT","d_ln_MARGEN_PAR")) {
  cat("\n---", v, "---\n")
  print(ur.df(na.omit(datos[[v]]), type = "drift", selectlags = "AIC"))
}

# ---------------- PASO 9: COINTEGRACIÓN (Engle-Granger) ----------------
cat("\n=== ENGLE-GRANGER (residuos del MCO en niveles) ===\n")
modelo_niv <- lm(ln_FBKF ~ ln_CARTERA + ln_CONCENT + ln_MARGEN_PAR, data = datos)
print(ur.df(residuals(modelo_niv), type = "none", selectlags = "AIC"))

# ---------------- PASO 10: CHOW (quiebre COVID 2020-03) ----------------
idx_covid <- which(datos$PERIODO == as.Date("2020-03-01"))
cat("\n=== CHOW 2020-03 (COVID) ===\n")
print(sctest(d_ln_FBKF ~ d_ln_CARTERA + d_ln_CONCENT + d_ln_MARGEN_PAR,
             data = datos, type = "Chow", point = idx_covid))

# ---------------- PASO 11: ROBUSTEZ SIN COVID (2001–2019 + 2022–2026) ----------------
datos_sin_covid <- datos %>%
  filter(PERIODO < as.Date("2020-01-01") | PERIODO > as.Date("2021-12-31"))

modelo_sin_covid <- lm(d_ln_FBKF ~ d_ln_CARTERA + d_ln_CONCENT + d_ln_MARGEN_PAR,
                       data = datos_sin_covid)
cat("\n=== MODELO ROBUSTEZ: SIN COVID ===\n")
print(summary(modelo_sin_covid))
cat("\nNewey-West sin COVID:\n")
print(coeftest(modelo_sin_covid, vcov = NeweyWest(modelo_sin_covid, lag = 4, prewhite = FALSE)))





# ---------------- PASO 12: ROBUSTEZ ALTERNATIVA (solo hoja "Data") ----------------
# Nota: el desglose margen financiero puro vs. comisiones NO es posible con esta hoja
# porque solo contiene MARGEN_PAR agregado. Se testean tres especificaciones
# alternativas para evaluar robustez de los coeficientes principales.

# 12A. Modelo dinámico: MARGEN_PAR rezagado (¿hay efecto con retardo?)
datos <- datos %>%
  mutate(
    L_d_ln_MARGEN_PAR = dplyr::lag(d_ln_MARGEN_PAR, 1),
    L_d_ln_CONCENT    = dplyr::lag(d_ln_CONCENT, 1),
    L_d_ln_CARTERA    = dplyr::lag(d_ln_CARTERA, 1)
  )

modelo_rezago <- lm(d_ln_FBKF ~ d_ln_CARTERA + d_ln_CONCENT + d_ln_MARGEN_PAR +
                      L_d_ln_CARTERA + L_d_ln_CONCENT + L_d_ln_MARGEN_PAR,
                    data = datos)
cat("\n=== 12A. MODELO CON REZAGOS (dinámica) ===\n")
print(summary(modelo_rezago))
cat("\nNewey-West con rezagos:\n")
print(coeftest(modelo_rezago, vcov = NeweyWest(modelo_rezago, lag = 4, prewhite = FALSE)))

# 12B. Modelo con tendencia determinista (controla drift)
datos <- datos %>% mutate(t = row_number())
modelo_trend <- lm(d_ln_FBKF ~ d_ln_CARTERA + d_ln_CONCENT + d_ln_MARGEN_PAR + t,
                   data = datos)
cat("\n=== 12B. MODELO CON TENDENCIA DETERMINISTA ===\n")
print(summary(modelo_trend))
cat("\nNewey-West con tendencia:\n")
print(coeftest(modelo_trend, vcov = NeweyWest(modelo_trend, lag = 4, prewhite = FALSE)))

# 12C. Modelo con dummies estacionales (mes)
datos <- datos %>% mutate(mes = factor(month(PERIODO)))
modelo_estac <- lm(d_ln_FBKF ~ d_ln_CARTERA + d_ln_CONCENT + d_ln_MARGEN_PAR + mes,
                   data = datos)
cat("\n=== 12C. MODELO CON ESTACIONALIDAD (dummies mensuales) ===\n")
print(summary(modelo_estac))
cat("\nNewey-West con estacionalidad:\n")
print(coeftest(modelo_estac, vcov = NeweyWest(modelo_estac, lag = 4, prewhite = FALSE)))

# 12D. Modelo con interacción COVID (¿cambió la relación después de 2020?)
datos <- datos %>% mutate(D_covid = ifelse(PERIODO >= as.Date("2020-03-01"), 1, 0))
modelo_covid <- lm(d_ln_FBKF ~ d_ln_CARTERA + d_ln_CONCENT + d_ln_MARGEN_PAR +
                     D_covid + d_ln_CONCENT:D_covid + d_ln_MARGEN_PAR:D_covid,
                   data = datos)
cat("\n=== 12D. MODELO CON INTERACCIÓN COVID ===\n")
print(summary(modelo_covid))
cat("\nNewey-West con interacción COVID:\n")
print(coeftest(modelo_covid, vcov = NeweyWest(modelo_covid, lag = 4, prewhite = FALSE)))


# ¿La concentración causa el margen? (eslabón 2 → 3)
modelo_margen <- lm(d_ln_MARGEN_PAR ~ d_ln_CONCENT, data = datos)
summary(modelo_margen)
coeftest(modelo_margen, vcov = NeweyWest(modelo_margen, lag = 4, prewhite = FALSE))

# ---------------- PASO 13: ACF/PACF ----------------
par(mfrow = c(2, 2))
acf(residuals(modelo),          main = "ACF — Modelo principal 2001–2026")
pacf(residuals(modelo),         main = "PACF — Modelo principal 2001–2026")
acf(residuals(modelo_sin_covid),main = "ACF — Sin COVID")
pacf(residuals(modelo_sin_covid),main = "PACF — Sin COVID")
par(mfrow = c(1, 1))

# ---------------- PASO 14: GRÁFICOS ----------------
# (a) Series en logaritmos
plot(datos$PERIODO, datos$ln_FBKF, type = "l", col = "blue",
     xlab = "Fecha", ylab = "Logaritmo",
     main = "Evolución de las variables (logaritmos, 2001–2026)")
lines(datos$PERIODO, datos$ln_CARTERA,    col = "red")
lines(datos$PERIODO, datos$ln_CONCENT,    col = "green")
lines(datos$PERIODO, datos$ln_MARGEN_PAR, col = "orange")
legend("topleft",
       legend = c("ln(FBKF)","ln(CARTERA)","ln(CONCENT)","ln(MARGEN_PAR)"),
       col = c("blue","red","green","orange"), lty = 0.3)

# (b) Residuos del modelo principal
par(mfrow = c(2, 2))
plot(modelo)
par(mfrow = c(1, 1))

# ---------------- FIN ----------------
cat("\n============================================================\n")
cat("MODELO PRINCIPAL: modelo (2001–2026)\n")
cat("ROBUSTEZ: modelo_sin_covid, modelo_desg (2001–2016)\n")
cat("============================================================\n")
