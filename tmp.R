gg <- "Tcf23"

dt1 <- dmu[Geneid == gg, ]
dt1$log2mu <- log2(dt1$mu)
dt1
dt1[, log2W2ConUVB := log2mu[trt == "CON"] - log2mu[trt == "UVB"],
    by = "time"]
dt1[, log2W2SfnUVB := log2mu[trt == "SFN"] - log2mu[trt == "UVB"],
    by = "time"]
dt1

res_con_uvb_week2[res_con_uvb_week2@rownames == gg, ]
res_sfn_uvb_week2[res_sfn_uvb_week2@rownames == gg, ]

res_con_uvb_week15[res_con_uvb_week15@rownames == gg, ]
res_sfn_uvb_week15[res_sfn_uvb_week15@rownames == gg, ]