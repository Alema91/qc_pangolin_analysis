#!/usr/bin/env Rscript

################################################
## LOAD LIBRARIES                             ##
################################################

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(readxl, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(forcats, quietly = TRUE, warn.conflicts = FALSE)
library(reshape2, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)

### Not in o in ----

`%notin%` <- Negate(`%in%`)

# grid plots
library(magick, quietly = TRUE, warn.conflicts = FALSE)
library(kableExtra, quietly = TRUE, warn.conflicts = FALSE)
library(knitr, quietly = TRUE, warn.conflicts = FALSE)
library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
library(ggpubr, quietly = TRUE, warn.conflicts = FALSE)

################################################
## DATA          ###############################
################################################

# data_pangolin<- read.csv2("output_pangolin/lineages_16_lab.csv", sep = ",")
# data_pangolin[data_pangolin$taxon == "COD_2122_1", ]

## data

data_ss <- read_xlsx("files/estadisticas_consensos.xlsx", sheet = 6, col_types = c(rep("text", 4), rep("numeric", 9)))

## data sequencing

data_seq <- read_xlsx("files/estadisticas_consensos.xlsx", sheet = 5)
data_seq$var_sequencing_platforms <- factor(data_seq$var_sequencing_platforms, levels = c("Illumina MiSeq", "Illumina iSeq", "Illumina NextSeq"))

bar_plataformas <- data_seq %>%
    count(var_sequencing_platforms) %>%
    ggbarplot(., "var_sequencing_platforms", "n",
        fill = "#0073C299",
        ylab = "Laboratorios",
        label = TRUE, label.pos = "out", lab.pos = c("in"),
        xlab = FALSE,
        ggtheme = theme_minimal()
    ) + rotate_x_text(angle = 45)

ggexport(bar_plataformas, filename = "plots/multi_page_sequencing_instruments.pdf")

data_ct <- melt(data_seq[, c(1, 2, 8:11)], id.var = c("ID", "Sample ID"))
data_ct$value <- as.numeric(data_ct$value)
# Box plot Cts
bxp_ct <- ggboxplot(data_ct,
    x = "Sample ID", y = "value",
    title = "Cts",
    fill = "#0073C299",
    add = "jitter",
    ylab = FALSE,
    xlab = FALSE,
    ggtheme = theme_minimal()
) +
    font("xy.text", size = 12, color = "gray", face = "bold")

data_seq_parsed <- read_xlsx("files/estadisticas_consensos.xlsx", sheet = 6)
nombres_columnas_t2 <- c("Plataformas de secuenciación", "Diagnóstico SARS-CoV-2", "Kit de libreria")

data_seq_parsed %>%
    kbl(align = "rrr", col.names = nombres_columnas_t2, "html") %>%
    kable_classic(full_width = F, font_size = 15) %>%
    pack_rows("Illumina iSeq (27.3%)", 1, 3) %>%
    pack_rows("Illumina MiSeq (68.2%)", 4, 9) %>%
    pack_rows("Illumina NextSeq (4.5%)", 10, 10)


# Arrange over multiple pages

multipage <- ggarrange(ggarrange(coverage_f, vio_ns_f, ncol = 2, labels = c("A", "B"), legend = F),
    ggarrange(vio_variants_f, labels = "C", common.legend = TRUE, legend = "bottom"),
    nrow = 2
)
ggexport(multipage, filename = "plots/multi_page_sequencing.pdf")


# total
df_melt <- melt(data_ss[, c(2, 4:length(data_ss))], id.var = c("group", "type"))

# Box plot Total reads
bxp_reads <- ggboxplot(df_melt[df_melt$variable == "totalreads", ],
    x = "variable", y = "value",
    color = "type",
    title = "Lecturas totales",
    palette = "jco",
    add = "jitter",
    ylab = FALSE,
    xlab = FALSE,
    ggtheme = theme_minimal()
) +
    font("xy.text", size = 12, color = "gray", face = "bold") +
    rremove("x.text")

bxp_reads_f <- facet(bxp_reads,
    facet.by = c("group"),
    labeller = "label_value",
    short.panel.labs = FALSE, # Allow long labels in panels
    panel.labs.background = list(fill = "steelblue", color = "black")
)

# violin plot Read host %
vio_host <- ggviolin(df_melt[df_melt$variable == "readshost", ],
    x = "variable", y = "value",
    col = "type",
    title = "Lecturas hospedador (%)",
    palette = "jco",
    add = "jitter",
    ylab = FALSE,
    xlab = FALSE,
    ylim = c(0, 100),
    ggtheme = theme_minimal()
) +
    font("xy.text", size = 12, color = "gray", face = "bold") +
    rremove("x.text")

vio_host_f <- facet(vio_ns,
    facet.by = c("group"),
    labeller = "label_value",
    short.panel.labs = FALSE,
    panel.labs.background = list(fill = "steelblue", color = "black")
)

# violin plot Read virus %
vio_virus <- ggviolin(df_melt[df_melt$variable == "readsvirus", ],
    x = "variable", y = "value",
    col = "type",
    title = "Lecturas SARS-CoV-2 (%)",
    palette = "jco",
    add = "jitter",
    ylab = FALSE,
    xlab = FALSE,
    ylim = c(0, 100),
    ggtheme = theme_minimal()
) +
    font("xy.text", size = 12, color = "gray", face = "bold") +
    rremove("x.text")

vio_virus_f <- facet(vio_virus,
    facet.by = c("group"),
    labeller = "label_value",
    short.panel.labs = FALSE,
    panel.labs.background = list(fill = "steelblue", color = "black")
)

# violin plot Read unmapped %
vio_unmapped <- ggviolin(df_melt[df_melt$variable == "unmapedreads", ],
    x = "variable", y = "value",
    col = "type",
    title = "Lecturas no mapeadas",
    palette = "jco",
    add = "jitter",
    ylab = FALSE,
    xlab = FALSE,
    ylim = c(0, 100),
    ggtheme = theme_minimal()
) +
    font("xy.text", size = 12, color = "gray", face = "bold") +
    rremove("x.text")

vio_unmapped_f <- facet(vio_unmapped,
    facet.by = c("group"),
    labeller = "label_value",
    short.panel.labs = FALSE,
    panel.labs.background = list(fill = "steelblue", color = "black")
)

# plot coverage

coverage <- df_melt[df_melt$variable == "medianDPcoveragevirus", ]
coverage$value <- log10(coverage$value)

p_coverage <- ggboxplot(coverage,
    y = "value",
    palette = "futurama",
    ylab = "log(mediana de la cobertura)",
    xlab = FALSE,
    ggtheme = theme_minimal()
) +
    geom_jitter(aes(color = group)) +
    rremove("x.text")
coverage_f <- facet(p_coverage,
    facet.by = c("type"),
    labeller = "label_value",
    short.panel.labs = FALSE,
    panel.labs.background = list(fill = "steelblue", color = "black")
)

# violin plot Ns
vio_ns <- ggviolin(df_melt[df_melt$variable == "Ns10x", ],
    y = "value",
    palette = "futurama",
    ylab = "% Ns",
    xlab = FALSE,
    ylim = c(0, 100),
    ggtheme = theme_minimal()
) +
    geom_jitter(aes(color = group)) +
    font("xy.text", size = 12, color = "gray", face = "bold") +
    rremove("x.text")

vio_ns_f <- facet(vio_ns,
    facet.by = c("type"),
    labeller = "label_value",
    short.panel.labs = FALSE,
    panel.labs.background = list(fill = "steelblue", color = "black")
)

# violin plot variants
vio_variants <- ggviolin(df_melt[df_melt$variable == "Variantsinconsensusx10" | df_melt$variable == "MissenseVariants", ],
    y = "value",
    palette = "futurama",
    ylab = "Mutaciones",
    xlab = FALSE,
    ggtheme = theme_minimal()
) +
    geom_jitter(aes(color = group)) +
    font("xy.text", size = 12, color = "gray", face = "bold") +
    rremove("x.text")

vio_variants_f <- facet(vio_variants,
    facet.by = c("type", "variable"),
    labeller = "label_value",
    short.panel.labs = FALSE,
    panel.labs.background = list(fill = "steelblue", color = "black")
)

# Arrange over multiple pages

multipage <- ggarrange(ggarrange(coverage_f, vio_ns_f, ncol = 2, labels = c("A", "B"), legend = F),
    ggarrange(vio_variants_f, labels = "C", common.legend = TRUE, legend = "bottom"),
    nrow = 2
)
ggexport(multipage, filename = "plots/multi_page_sequencing.pdf")

# Table

data_bioinfo <- read_xlsx("files/estadisticas_consensos.xlsx", sheet = 3)
f_d_bioinfo <- data_bioinfo[, c(2, 4:9)]
f_d_bioinfo$plataformas <- factor(f_d_bioinfo$plataformas, levels = c("Illumina MiSeq", "Illumina iSeq", "Illumina NextSeq"))
f_d_bioinfo$bioinformatic_protocol <- factor(f_d_bioinfo$bioinformatic_protocol, levels = c("DRAGEN COVID", "In-house protocol", "SeqCOVID", "DeepChek pipeline"))
f_d_bioinfo$Preprocessing <- factor(f_d_bioinfo$Preprocessing, levels = c("fastp", "DRAGEN COVID trimmer", "NA", "artic protocol", "trimmomatic"))
f_d_bioinfo$Mapping <- factor(f_d_bioinfo$Mapping, levels = c("bwa mem", "NA", "DRAGEN COVID mapping"))
f_d_bioinfo$Assembly <- factor(f_d_bioinfo$Assembly, levels = c("NA", "DRAGEN COVID", "None"))
f_d_bioinfo$Variant_Calling <- factor(f_d_bioinfo$Variant_Calling, levels = c("ivar variants", "DRAGEN COVID Variant Caller", "NA", "bcftools", "freebayes"))
f_d_bioinfo$Consensus <- factor(f_d_bioinfo$Consensus, levels = c("ivar consensus", "DRAGEN COVID", "bcftools consensus", "NA"))

bar_plataformas <- f_d_bioinfo %>%
    count(plataformas) %>%
    ggbarplot(., "plataformas", "n",
        fill = "#0073C299",
        ylab = "Laboratorios",
        label = TRUE, label.pos = "out", lab.pos = c("in"),
        xlab = FALSE,
        ggtheme = theme_minimal()
    ) + rotate_x_text(angle = 45)

bar_bioinfopro <- f_d_bioinfo %>%
    count(bioinformatic_protocol) %>%
    ggbarplot(., "bioinformatic_protocol", "n",
        fill = "#0073C299",
        ylab = "Laboratorios",
        label = TRUE, label.pos = "out", lab.pos = c("in"),
        xlab = FALSE,
        ggtheme = theme_minimal()
    ) + rotate_x_text(angle = 45)

bar_pre <- f_d_bioinfo %>%
    count(Preprocessing) %>%
    ggbarplot(., "Preprocessing", "n",
        fill = "#0073C299",
        ylab = "Laboratorios",
        label = TRUE, label.pos = "out", lab.pos = c("in"),
        xlab = FALSE,
        ggtheme = theme_minimal()
    ) + rotate_x_text(angle = 45)

bar_mapping <- f_d_bioinfo %>%
    count(Mapping) %>%
    ggbarplot(., "Mapping", "n",
        fill = "#0073C299",
        ylab = "Laboratorios",
        label = TRUE, label.pos = "out", lab.pos = c("in"),
        xlab = FALSE,
        ggtheme = theme_minimal()
    ) + rotate_x_text(angle = 45)

bar_assembly <- f_d_bioinfo %>%
    count(Assembly) %>%
    ggbarplot(., "Assembly", "n",
        fill = "#0073C299",
        ylab = "Laboratorios",
        label = TRUE, label.pos = "out", lab.pos = c("in"),
        xlab = FALSE,
        ggtheme = theme_minimal()
    ) + rotate_x_text(angle = 45)

bar_vcaller <- f_d_bioinfo %>%
    count(Variant_Calling) %>%
    ggbarplot(., "Variant_Calling", "n",
        fill = "#0073C299",
        ylab = "Laboratorios",
        label = TRUE, label.pos = "out", lab.pos = c("in"),
        xlab = FALSE,
        ggtheme = theme_minimal()
    ) + rotate_x_text(angle = 45)

bar_consensus <- f_d_bioinfo %>%
    count(Consensus) %>%
    ggbarplot(., "Consensus", "n",
        fill = "#0073C299",
        ylab = "Laboratorios",
        label = TRUE, label.pos = "out", lab.pos = c("in"),
        xlab = FALSE,
        ggtheme = theme_minimal()
    ) + rotate_x_text(angle = 45)

ggexport(bar_plataformas, filename = "plots/plat_sequencing.png", ncol = 1, nrow = 1)
multipage <- ggarrange(bar_pre, bar_mapping, bar_vcaller, bar_consensus, labels = c("A", "B", "C", "D"), heights = c(0.5, 0.5), widths = c(1, 1), ncol = 2, nrow = 2)
ggexport(multipage, filename = "plots/multi_page_bioinfo.pdf")

# Tabla bioinfo

software_1 <- c(rep("ivar", 3))
params <- c(
    "minimum quality for consensus calling = 20, minimum frequency to consider fixed a SNP = 0.8, minimum position depth = 30 (ambiguous base otherwise)",
    "samtools mpileup -aa -A -d 0 -B -Q 0 X.trim.sort.bam | ivar consensus -p X -q 20 -t 0.8 -m 30 -n N",
    "-q 20 -t 0.03 -m 15"
)
df_tabla_bioinfo_1 <- data.frame(software_1, params)

software_2 <- c(rep("bcftools consensus", 2), rep("DRAGEN COVID", 2), rep("ivar", 1))
QC_metrics <- c(
    "consensus with 75% calling variants",
    "> 3000 N's",
    "Virus amplicons <5; Human control genes <4; Pangolin: Max Ambiguous Rate=0.5",
    ">90% viral amplicons detected; Coverage threshold = 30",
    "COV >= 90% 30X; --max-ambig 0.3"
)
df_tabla_bioinfo_2 <- data.frame(software_2, QC_metrics)

# main_title <- "Generación de consensos"
# subtitle <- paste0(
#    "Los consensos se han generado por lo softwares ivar (ivar consensus, 50%), bcftools (bcftools consensus, 16.7%) y DRAGEN (16.7%).",
#    " Un 16.7% de los laboratorios han dejado vacío este campo."
# ) %>%
#    strwrap(width = 80) %>%
#    paste(collapse = "\n")

nombres_columnas_t1_bioinfo <- c("Software", "Parámetros")
df_tabla_bioinfo_1 %>%
    kbl(align = "rl", col.names = nombres_columnas_t1_bioinfo, "html") %>%
    row_spec(0, bold = T) %>%
    kable_classic(full_width = F, font_size = 15)

nombres_columnas_t2_bioinfo <- c("Software", "Métricas control de calidad")
df_tabla_bioinfo_2 %>%
    kbl(align = "rl", col.names = nombres_columnas_t2_bioinfo, "html") %>%
    row_spec(0, bold = T) %>%
    kable_classic(full_width = F, font_size = 15)

# tabla resumen

data_ss <- read_xlsx("files/estadisticas_consensos.xlsx", sheet = 1)

nombres_columnas <- c("Grupo", "Muestra", "v.3.1.16", "v3.1.17", "v3.1.18", "v3.1.19", "v3.1.16", "v3.1.17", "v3.1.18", "v3.1.19", "Problemáticas", "Causas")
nombres_columnas_t1 <- c("Grupo", "Muestra", "v.3.1.16", "v3.1.17", "v3.1.18", "v3.1.19", "v3.1.16", "v3.1.17", "v3.1.18", "v3.1.19", "Problemáticas")
# filter(., Causas == "Version Pangolin")

data_ss %>%
    arrange(., Causas) %>%
    select(., !Causas) %>%
    kbl(align = "rrrrrrrrrrll", col.names = nombres_columnas_t1, "html") %>%
    kable_classic(full_width = F, font_size = 15) %>%
    row_spec(0, bold = T) %>%
    add_header_above(c(" " = 2, "Laboratorios" = 4, "Viralrecon" = 4, " " = 1), bold = T) %>%
    pack_rows("Protocolo bioinformático", 1, 14) %>%
    pack_rows("Secuenciación", 16, 17) %>%
    pack_rows("Versión Pangolin", 18, 23)

data_ss %>%
    arrange(., Causas) %>%
    filter(., Causas == "Version Pangolin") %>%
    kbl(align = "rrrrrrrrrrll", col.names = nombres_columnas, "html") %>%
    kable_classic(full_width = F, font_size = 15) %>%
    add_header_above(c(" " = 2, "Laboratorios" = 4, "Viralrecon" = 4, " " = 2))

data_ss %>%
    arrange(., Causas) %>%
    filter(., Causas == "Protocolo bioinfo") %>%
    kbl(align = "rrrrrrrrrrll", col.names = nombres_columnas, "html") %>%
    kable_classic(full_width = F, font_size = 15) %>%
    add_header_above(c(" " = 2, "Laboratorios" = 4, "Viralrecon" = 4, " " = 2))

data_ss %>%
    arrange(., Causas) %>%
    filter(., Causas == "Secuenciación") %>%
    kbl(align = "rrrrrrrrrrll", col.names = nombres_columnas, "html") %>%
    kable_classic(full_width = F, font_size = 15) %>%
    add_header_above(c(" " = 2, "Laboratorios" = 4, "Viralrecon" = 4, " " = 2))


# Ejemplos

# Ejemplos secuenciación

columnas <- c("Laboratorio", "Muestra", "Lecturas mapeadas", "% mediana de la cobertura", "% cobertura > 10x", "SNPs", "INDELs", "Variantes Missense", "Ns por 100kb consensus", "%Ns", "Pangolin linaje")
values_2140 <- c("COD_2140", "muestra 6", "28954", "39", "67", "17", "NA", "10", "42597.73", "42", "None")

data_ejemplo_seq <- data.frame(t(values_2140))
colnames(data_ejemplo_seq) <- columnas

data_ejemplo_seq %>%
    kbl(align = "rrlllllllll", "html") %>%
    kable_classic(full_width = F, font_size = 15)

## DATA resumen
niveles_samples <- c(
    "muestra_1",
    "muestra_3",
    "muestra_4",
    "muestra_5",
    "muestra_6",
    "muestra_8",
    "muestra_9",
    "muestra_10"
)

data_resumen <- read.csv2("files/df_resumen.csv", sep = "\t")
data_resumen$muestra <- factor(data_resumen$muestra, levels = niveles_samples)

# geom_text(stat = "count", aes(label = ..count..)) +

ggplot(data_resumen, aes(muestra, fill = Problemáticas)) +
    geom_bar() +
    guides(fill = guide_legend(title = "")) +
    labs(y = "Causas", x = "", title = "") +
    # geom_text(aes(label = ..count..), stat = "count", position = "fill", vjust = -5) +
    geom_text(stat = "count", aes(label = after_stat(count)), size = 6) +
    theme(
        text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(),
        legend.title = element_text(),
        legend.text = element_text()
    )
ggsave("plots/causas_pangolin.png", width = 60, height = 40, units = "cm")

ggplot(data_resumen, aes(grupo, fill = Problemáticas)) +
    geom_bar() +
    guides(fill = guide_legend(title = "")) +
    labs(y = "Causas", x = "", title = "", size = 12) +
    geom_text(stat = "count", aes(label = ..count..), size = 6.5, vjust = -0.1) +
    theme(
        text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(),
        legend.title = element_text(),
        legend.text = element_text()
    )
ggsave("plots/grupo_causas_pangolin.png", width = 60, height = 40, units = "cm")

ggplot(data_resumen, aes(Causas, fill = Problemáticas)) +
    geom_bar() +
    guides(fill = guide_legend(title = "")) +
    labs(y = "Causas", x = "", title = "", size = 12) +
    geom_text(stat = "count", aes(label = ..count..), size = 6.5, vjust = -0.1) +
    theme(
        text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_text(),
        legend.title = element_text(),
        legend.text = element_text()
    )
ggsave("plots/resumen_causas.png", width = 60, height = 40, units = "cm")
