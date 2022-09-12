## Versiones de Pangolin

Hemos usado las siguientes versiones de Pangolin para discernir si existen diferencias a nivel de versión, análisis bioinformático y secuenciación:

- v3.1.16
- v3.1.17
- v3.1.18
- v3.1.19
- v3.1.20

Dividiremos el análisis por muestras para 30 laboratorios (Illumina y Nanopore) analizados con las versiones anteriores de Pangolin con sus análisis y el análisis propio de la Unidad de Bioinformática (Viralrecon).

### Muestra 1

Los resultados son homogeneos para todos los casos, con el linaje B.1.1.7 menos para los siguientes casos:

- **COD_2122, v3.1.17**: Q.4
Segun [COGUK](https://cov-lineages.org/lineage_list.html) es un sublinaje de la B.1.1.7 con una mutación en la spike (spike:H681R). El issue relacionado [Pango_designation: 176](https://github.com/cov-lineages/pango-designation/issues/176) y añadido a la versión v1.2.57 de pango-designation [Pango_designation v1.2.57](https://github.com/cov-lineages/pango-designation/releases/tag/v1.2.57).

Parece que es un problema del PANGOLearn (Tabla X).

| taxon    | sample | class        | pangolin_version | lineage | conflict | ambiguity_score | scorpio_call         | scorpio_support | scorpio_conflict | version         | pangoLEARN_version | pango_version | status    | note                                                                      |
|----------|--------|--------------|------------------|---------|----------|-----------------|----------------------|-----------------|------------------|-----------------|--------------------|---------------|-----------|---------------------------------------------------------------------------|
| COD_2122 |      1 | laboratories | 3.1.16           | B.1.1.7 | 0        | 0.8388268881    | Alpha (B.1.1.7-like) | 0.8261          | 0.0435           | PLEARN-v1.2.97  | 2021-11-18         | v1.2.97       | passed_qc | scorpio call: Alt alleles 19; Ref alleles 1; Amb alleles 3; Oth alleles 0 |
| COD_2122 |      1 | viralrecon   | 3.1.16           | B.1.1.7 | 0        | 0.9135502671    | Alpha (B.1.1.7-like) | 0.8261          | 0.1304           | PLEARN-v1.2.97  | 2021-11-18         | v1.2.97       | passed_qc | scorpio call: Alt alleles 19; Ref alleles 3; Amb alleles 1; Oth alleles 0 |
| COD_2122 |      1 | laboratories | 3.1.17           | Q.4     |        0 |    0.8523217457 | Alpha (B.1.1.7-like) |           0.913 |           0.0435 | PLEARN-v1.2.101 |         2021-11-25 | v1.2.101      | passed_qc | scorpio call: Alt alleles 21; Ref alleles 1; Amb alleles 1; Oth alleles 0 |
| COD_2122 |      1 | viralrecon   | 3.1.17           | Q.4     |        0 |    0.9174246051 | Alpha (B.1.1.7-like) |          0.8696 |           0.1304 | PLEARN-v1.2.101 |         2021-11-25 | v1.2.101      | passed_qc | scorpio call: Alt alleles 20; Ref alleles 3; Amb alleles 0; Oth alleles 0 |

Las mutaciones observadas para el sublinaje según el issue son las siguientes:

```
- ORF1a:A2320V (nsp3:A1502V, nucleotide C7224T)
- S:H681R (A23604G)
```

Aquí se encuentran las mutaciones compartirdas entre B.1.1.7 y Q.4:

![image](q.4vsb1.1.7.png)

He buscado estas dos mutaciones en la longtable de variantes de la muestra 1 de COD_2122 y la he comparado con los resultados de la muestra 1 de COD_2107 (un laboratorio que lo ha hecho todo muy bien). Aparece la mutacion a nivel de nucleótido C23604T (misssense effect) en el gen S (aunque el cambio recogido es una C>T) pero no la mutacion C7224T en el gen ORF1a.

[longtable_variantes_muestra_1](https://docs.google.com/spreadsheets/d/1eq7tpPi9YkRsEP3AQXTR_suEM46nU4i9sURXi6IpCFk/edit#gid=2059206970)
[longtable_outputpangolin_muestra_1](https://docs.google.com/spreadsheets/d/1eq7tpPi9YkRsEP3AQXTR_suEM46nU4i9sURXi6IpCFk/edit#gid=39736565)

Además, he realizado el alineamiento con [clustalo](https://github.com/hybsearch/clustalo) de los fasta de la muestra 1 tanto de los laboratorios (COD_2122 y COD_2107 como control) como los de viralrecon y la referencia (NC_045512.2). Se puede visualizar en el siguiente enlace:

[NCBI visual MSA](https://www.ncbi.nlm.nih.gov/projects/msaviewer/?key=ZtX8DHrXpQ4J-esJKujd941PP-BgkW6UYpJKhF6ATK7doDffpeWn95hvrRKglurru_Pm5_jDo8bk3PDR9tf6zMjl9evZ1_M,O4ihUSeK-FNUpLZUd7WAqtARFeJKk0SWSJBghnSCZqz3oh3dj-eXsPhEzcXPACF9cGUtcTNVaFAvSjtHPUExWgNzPn0SQTg)

- **COD_2131, todas las versiones en Viralrecon**

EL fastq de la muestra 1 (COD_2131_1_R2.fastq.gz) que han subido al sFTP esta corrupto y por lo tanto el resultado de Viralrecon es NA.

- **COD_2137, todas las versiones en Viralrecon**

El output de Pangolin para esas muestras nos dice que aunque pangoLEARN si asigna correctamente el linaje, este no es soportado por scorpio. 

Resultado de Pangolin:

```
COD_2137_1_NC_045512.2,None,,,,,,PLEARN-v1.2.123,3.1.20,2022-02-28,v1.2.123,passed_qc,pangoLEARN lineage assignment B.1.1.7 was not supported by scorpio
```

Analizamos los vcf en busca de la posible diferencia entre los fastas procedentes de los laboratorios y los procedentes de Viralrecon:

- COD_2137_1: 153 mutaciones
- Viralrecon_2137_1: 27 mutaciones

He realizado el alineamiento como el caso de la 2122 y se puede visualizar en el siguiente enlace:

[NCBI visual MSA](https://www.ncbi.nlm.nih.gov/projects/msaviewer/?key=qBsywrQZa8DHNyXH5CYTOUOA5hC5Ybdku2KTdIdwlV4EUO10yE63zXu9TpiT41GeAIZdkkO2GLNfqUukTaJBuXOQTp5iokg,QvHYKF7zgSot3c8tDsz506lqDPhTiV2MUYp5nG2Yf7buuAecIqb1hIzWufWTo_veqsb30un2svP16eHk5-Lr-dnQ5N7I4uI)

LLama la atención la región comprendida entre la posición 5000 a 5283 con una región con N. La misma muestra 1 del COD_2137 analizada con viralrecon recupera información de la región 5070 a la 5160. Esta región se encuentra el orf1ab. Sin embargo, tenemos estos resultados de viralrecon:

| taxon    | sample | class      | pangolin_version | lineage | conflict | ambiguity_score | scorpio_call | scorpio_support | scorpio_conflict | version         | pangoLEARN_version | pango_version | status    | note                                                               | Reason         |
|----------|--------|------------|------------------|---------|----------|-----------------|--------------|-----------------|------------------|-----------------|--------------------|---------------|-----------|--------------------------------------------------------------------|----------------|
| COD_2137 |      1 | viralrecon | 3.1.16           | B.1.1   |        0 |    0.9923423423 |              |                 |                  | PLEARN-v1.2.97  |         2021-11-18 | v1.2.97       | passed_qc |                                                                    | 2137_1 problem |
| COD_2137 |      1 | viralrecon | 3.1.17           | None    |          |                 |              |                 |                  | PLEARN-v1.2.101 |         2021-11-25 | v1.2.101      | passed_qc | pangoLEARN lineage assignment B.1.1.7 was not supported by scorpio | 2137_1 problem |

El mismo resultado de la versión 3.1.17 lo tendríamos para la version 3.1.18 y 3.1.19.

