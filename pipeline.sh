#!/usr/bin/env bash
set -eo pipefail

source "Config.txt"
source /opt/conda/etc/profile.d/conda.sh
mkdir -p "${WORKDIR}"/{raw_data,ref,qc/raw,qc/trimmed,qc/multiqc,trimmed,aligned,qualimap,coverage,variants}
cd "${WORKDIR}"

REF_FNA="ref/reference.fna"
REF_GTF="ref/reference.gtf"

# ==============================================================================
# 1. СКАЧИВАНИЕ РЕФЕРЕНСА И SRA
# ==============================================================================
if [[ ! -f "$REF_FNA" ]]; then
    echo "[ШАГ 1] Скачивание референсного генома"
    wget -qO- "$REF_FNA_URL" | gunzip -c > "$REF_FNA"
fi

if [[ ! -f "$REF_GTF" ]]; then
    echo "[ШАГ 1] Скачивание аннотации"
    wget -qO- "$REF_GTF_URL" | gunzip -c > "$REF_GTF"
fi

for SRA_ID in "${SRA_IDS[@]}"; do
    if [[ ! -f "raw_data/${SRA_ID}_1.fastq.gz" || ! -f "raw_data/${SRA_ID}_2.fastq.gz" ]]; then
        echo "-> Скачивание прочтений для образца: ${SRA_ID}"
        "$FASTERQDUMP_BIN" "$SRA_ID" -O raw_data/ -e "$THREADS" --split-files
        gzip -f "raw_data/${SRA_ID}_1.fastq" "raw_data/${SRA_ID}_2.fastq"
    fi
done

# ==============================================================================
# 2. FASTQC -> ТРИММИНГ -> FASTQC -> MULTIQC
# ==============================================================================
conda activate "$ENV_MAIN"

for SRA_ID in "${SRA_IDS[@]}"; do
    RAW_R1="raw_data/${SRA_ID}_1.fastq.gz"
    RAW_R2="raw_data/${SRA_ID}_2.fastq.gz"
    TRIM_R1="trimmed/${SRA_ID}_${TRIMMER}_1.trimmed.fastq.gz"
    TRIM_R2="trimmed/${SRA_ID}_${TRIMMER}_2.trimmed.fastq.gz"
    if [[ ! -f "qc/raw/${SRA_ID}_1_fastqc.html" || ! -f "qc/raw/${SRA_ID}_2_fastqc.html" ]]; then
        echo "-> [FastQC] Анализ сырых ридов образца: ${SRA_ID}"
        fastqc -t "$THREADS" -o qc/raw "$RAW_R1" "$RAW_R2"
    fi
    if [[ ! -f "$TRIM_R1" || ! -f "$TRIM_R2" ]]; then
        echo "-> [Тримминг] Обработка образца ${SRA_ID} через ${TRIMMER}"
        if [[ "$TRIMMER" == "fastp" ]]; then
            fastp -i "$RAW_R1" -I "$RAW_R2" -o "$TRIM_R1" -O "$TRIM_R2" \
                --detect_adapter_for_pe --thread "$THREADS" \
                --html "qc/trimmed/${SRA_ID}_fastp.html" --json "qc/trimmed/${SRA_ID}_fastp.json"
        elif [[ "$TRIMMER" == "trimmomatic" ]]; then
            trimmomatic PE -threads "$THREADS" -phred33 "$RAW_R1" "$RAW_R2" \
                "$TRIM_R1" "trimmed/${SRA_ID}_${TRIMMER}_1.unpaired.fastq.gz" \
                "$TRIM_R2" "trimmed/${SRA_ID}_${TRIMMER}_2.unpaired.fastq.gz" \
                ILLUMINACLIP:"${TRIMMOMATIC_ADAPTERS}":2:30:10 \
                LEADING:"$TRIM_LEADING" TRAILING:"$TRIM_TRAILING" \
                SLIDINGWINDOW:"$TRIM_SLIDINGWINDOW" MINLEN:"$TRIM_MINLEN"
        fi
    fi
    if [[ ! -f "qc/trimmed/${SRA_ID}_${TRIMMER}_1.trimmed_fastqc.html" || ! -f "qc/trimmed/${SRA_ID}_${TRIMMER}_2.trimmed_fastqc.html" ]]; then
        echo "-> [FastQC] Анализ оттримленных ридов образца: ${SRA_ID}"
        fastqc -t "$THREADS" -o qc/trimmed "$TRIM_R1" "$TRIM_R2"
    fi
done
echo "-> [MultiQC] Сборка общего отчёта по тримированным файлам"
multiqc qc/ -o qc/multiqc -n "all_samples_multiqc_report.html" --force --dirs


# ==============================================================================
# 3. ИНДЕКСАЦИЯ РЕФЕРЕНСА И ПОИСК КООРДИНАТ ФАКТОРА ПАТОГЕННОСТИ
# ==============================================================================
echo "[ШАГ 3] Индексация референса (${ALIGNER}) и поиск гена ${VIRULENCE_GENE}"
if [[ "$ALIGNER" == "bwa" ]]; then
    if [[ ! -f "${REF_FNA}.amb" || ! -f "${REF_FNA}.ann" || ! -f "${REF_FNA}.bwt" || ! -f "${REF_FNA}.pac" || ! -f "${REF_FNA}.sa" ]]; then
        bwa index "$REF_FNA"
    fi
elif [[ "$ALIGNER" == "bowtie2" ]]; then
    if [[ ( ! -f "ref/bowtie2_index.1.bt2" || ! -f "ref/bowtie2_index.2.bt2" || ! -f "ref/bowtie2_index.3.bt2" || ! -f "ref/bowtie2_index.4.bt2" || ! -f "ref/bowtie2_index.rev.1.bt2" || ! -f "ref/bowtie2_index.rev.2.bt2" ) &&
          ( ! -f "ref/bowtie2_index.1.bt2l" || ! -f "ref/bowtie2_index.2.bt2l" || ! -f "ref/bowtie2_index.3.bt2l" || ! -f "ref/bowtie2_index.4.bt2l" || ! -f "ref/bowtie2_index.rev.1.bt2l" || ! -f "ref/bowtie2_index.rev.2.bt2l" ) ]]; then
        bowtie2-build --threads "$THREADS" "$REF_FNA" "ref/bowtie2_index"
    fi
fi

if [[ ! -f "${REF_FNA}.fai" ]]; then
    samtools faidx "$REF_FNA"
fi
if [[ ! -f "ref/reference.dict" ]]; then
    "$GATK_BIN" CreateSequenceDictionary -R "$REF_FNA" -O "ref/reference.dict"
fi

REGION_STR=$(grep -i -w "$VIRULENCE_GENE" "$REF_GTF" | awk -F'\t' '$3=="gene" {print $1":"$4"-"$5; exit}')
echo "Ген ${VIRULENCE_GENE}, найден на координатах ${REGION_STR}"

# ==============================================================================
# 4. КАРТИРОВАНИЕ, QUALIMAP, ПОКРЫТИЕ, MARKDUPLICATES И КОЛЛИНГ
# ==============================================================================
ALL_BAMS=()
echo "[ШАГ 4] Картирование, покрытие, дубликаты, Qualimap и коллинг (${CALLER})"
for SRA_ID in "${SRA_IDS[@]}"; do
    conda activate "$ENV_MAIN"
    echo "Обработка образца: ${SRA_ID}"
    TRIM_R1="trimmed/${SRA_ID}_${TRIMMER}_1.trimmed.fastq.gz"
    TRIM_R2="trimmed/${SRA_ID}_${TRIMMER}_2.trimmed.fastq.gz"
    SORTED_BAM="aligned/${SRA_ID}_${TRIMMER}.${ALIGNER}.sorted.bam"
    RG_BAM="aligned/${SRA_ID}_${TRIMMER}.${ALIGNER}.rg.bam"
    DEDUP_BAM="aligned/${SRA_ID}_${TRIMMER}.${ALIGNER}.dedup.bam"

    # 4.1 Картирование с учетом кастомных параметров из конфига и сортировка
    if [[ ! -f "$SORTED_BAM" ]]; then
        if [[ "$ALIGNER" == "bwa" ]]; then
            bwa mem $BWA_PARAMS -t "$THREADS" "$REF_FNA" "$TRIM_R1" "$TRIM_R2" | \
                samtools sort -@ "$THREADS" -o "$SORTED_BAM"
        elif [[ "$ALIGNER" == "bowtie2" ]]; then
            bowtie2 $BOWTIE2_PARAMS -p "$THREADS" -x "ref/bowtie2_index" -1 "$TRIM_R1" -2 "$TRIM_R2" | \
                samtools sort -@ "$THREADS" -o "$SORTED_BAM"
        fi
    fi

    # 4.2 Запись заголовков @RG через AddOrReplaceReadGroups
    if [[ ! -f "$RG_BAM" ]]; then
        "$GATK_BIN" AddOrReplaceReadGroups -I "$SORTED_BAM" -O "$RG_BAM" \
            -RGID "${RG_ID}_${SRA_ID}" \
            -RGLB "${RG_LB}_${SRA_ID}" \
            -RGPL "$RG_PL" \
            -RGPU "$RG_PU" \
            -RGSM "${RG_SM}_${SRA_ID}" \
            --CREATE_INDEX true
    fi

    # 4.3 Расчёт покрытия участка гена вирулентности
    if [[ ! -f "coverage/${SRA_ID}_${TRIMMER}.${ALIGNER}.${VIRULENCE_GENE}.coverage.txt" ]]; then
        samtools coverage -r "$REGION_STR" "$RG_BAM" > "coverage/${SRA_ID}_${TRIMMER}.${ALIGNER}.${VIRULENCE_GENE}.coverage.txt"
    fi
    if [[ ! -f "coverage/${SRA_ID}_${TRIMMER}.${ALIGNER}.${VIRULENCE_GENE}.depth.tsv" ]]; then
        samtools depth -a -r "$REGION_STR" "$RG_BAM" > "coverage/${SRA_ID}_${TRIMMER}.${ALIGNER}.${VIRULENCE_GENE}.depth.tsv"
    fi

    # 4.4 Маркировка дубликатов
    if [[ ! -f "$DEDUP_BAM" || ! -f "aligned/${SRA_ID}_${TRIMMER}.${ALIGNER}.dup_metrics.txt" ]]; then
        "$GATK_BIN" MarkDuplicates -I "$RG_BAM" -O "$DEDUP_BAM" -M "aligned/${SRA_ID}_${TRIMMER}.${ALIGNER}.dup_metrics.txt" --CREATE_INDEX true
    fi
    if [[ ! -f "${DEDUP_BAM}.bai" ]]; then
        samtools index "$DEDUP_BAM"
    fi
    ALL_BAMS+=("$DEDUP_BAM")

    # 4.5 Отчёт Qualimap
    conda activate "$ENV_QUALIMAP"
    if [[ ! -f "qualimap/${SRA_ID}_${TRIMMER}_${ALIGNER}_bamqc/qualimapReport.html" ]]; then
        qualimap bamqc -bam "$RG_BAM" -gff "$REF_GTF" -outdir "qualimap/${SRA_ID}_${TRIMMER}_${ALIGNER}_bamqc" -nt "$THREADS" --java-mem-size=4G
    fi
    if [[ ! -f "qualimap/${SRA_ID}_${TRIMMER}_${ALIGNER}_dedup__bamqc/qualimapReport.html" ]]; then
        qualimap bamqc -bam "$DEDUP_BAM" -gff "$REF_GTF" -outdir "qualimap/${SRA_ID}_${TRIMMER}_${ALIGNER}_dedup__bamqc" -nt "$THREADS" --java-mem-size=4G
    fi
    echo "-> [Variant Calling] Коллинг вариантов через ${CALLER}"
    # 4.6 Коллинг вариантов (gatk или freebayes)
    conda activate "$ENV_MAIN"
    if [[ "$CALLER" == "gatk" ]]; then
        if [[ ! -f "variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.gatk.vcf.gz" ]]; then
            "$GATK_BIN" HaplotypeCaller -R "$REF_FNA" -I "$DEDUP_BAM" \
                -O "variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.gatk.vcf.gz" \
                -ploidy "$PLOIDY" --native-pair-hmm-threads "$THREADS"
        fi
    elif [[ "$CALLER" == "freebayes" ]]; then
        if [[ ! -f "variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.freebayes.vcf.gz" ]]; then
            conda activate "$ENV_FREEBAYES"
            if [[ ! -f "variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.freebayes.vcf" ]]; then
                freebayes -f "$REF_FNA" -p "$PLOIDY" "$DEDUP_BAM" > "variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.freebayes.vcf"
            fi
            conda activate "$ENV_MAIN"
            bgzip "variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.freebayes.vcf"
        fi
        if [[ ! -f "variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.freebayes.vcf.gz.tbi" ]]; then
            tabix -p vcf "variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.freebayes.vcf.gz"
        fi
    fi
done
conda activate "$ENV_MAIN"
echo "->Создание общего отчёта для картированных файлов"
multiqc qualimap/ -o qualimap/multiqc -n "all_samples_bamqc_multiqc_report.html" --force --dirs

# График deepTools plotCoverage
conda activate "$ENV_MAIN"
if [[ ! -f "coverage/${VIRULENCE_GENE}_${TRIMMER}_${ALIGNER}_plotCoverage.png" ]]; then
    echo "Построение графика deepTools plotCoverage по участку фактора вирулентности"
    plotCoverage -b "${ALL_BAMS[@]}" -n "$NUMOFSAMPL" --region "${REGION_STR/-/:}" --labels "${SRA_IDS[@]}" \
        --plotFile "coverage/${VIRULENCE_GENE}_${TRIMMER}_${ALIGNER}_plotCoverage.png" -p "$THREADS"
fi
if [[ ! -f "coverage/whole_${TRIMMER}_${ALIGNER}_plotCoverage.png" ]]; then
    echo "Построение общего графика deepTools plotCoverage"
    plotCoverage -b "${ALL_BAMS[@]}" -n "$NUMOFSAMPL" --labels "${SRA_IDS[@]}" \
        --plotFile "coverage/whole_${TRIMMER}_${ALIGNER}_plotCoverage.png" -p "$THREADS"
fi
# ==============================================================================
# 5. ОБЪЕДИНЕНИЕ И ИНДЕКСАЦИЯ VCF
# ==============================================================================
conda activate "$ENV_BCF"
VCF_FILES=()
for SRA_ID in "${SRA_IDS[@]}"; do
    VCF="variants/${SRA_ID}_${TRIMMER}.${ALIGNER}.${CALLER}.vcf.gz"
    if [[ ! -f "$VCF" ]]; then
        echo "Ошибка: не найден файл $VCF" >&2
        exit 1
    fi
    echo "-> Подготовка VCF: $VCF"
    bcftools index -f "$VCF"
    VCF_FILES+=("$VCF")
done
MERGED_VCF="variants/merged_${TRIMMER}_${ALIGNER}_${CALLER}.vcf.gz"
bcftools merge "${VCF_FILES[@]}" \
    -Oz -o "$MERGED_VCF" \
    --threads "$THREADS"
bcftools index -f "$MERGED_VCF"
# ==============================================================================
# 6. ФИЛЬТРАЦИЯ ОБЩЕГО VCF
# ==============================================================================
conda activate "$ENV_BCF"
PREFIX="variants/merged_${TRIMMER}_${ALIGNER}_${CALLER}"
FILTERED_VCF="${PREFIX}.filtered.vcf.gz"
ANNOTATED_RAW="${PREFIX}.filtered.ann.vcf"
ANNOTATED_VCF="${ANNOTATED_RAW}.gz"
WITH_ANN_VCF="${PREFIX}.with_annotation.vcf.gz"
WITHOUT_ANN_VCF="${PREFIX}.without_annotation.vcf.gz"
REPORT_DIR="variants/reports_${TRIMMER}_${ALIGNER}_${CALLER}"
mkdir -p "$REPORT_DIR"
echo "-> Фильтрация общего VCF"
bcftools view \
    -i "QUAL!=\".\" && QUAL>=${MIN_QUAL}" \
    "$MERGED_VCF" -Ou |
bcftools filter \
    -S . \
    -e "FMT/DP=\".\" | FMT/DP<${MIN_DP}" \
    -Ou |
bcftools view \
    -i 'GT="alt"' \
    -Oz -o "$FILTERED_VCF"
bcftools index -f "$FILTERED_VCF"
# ==============================================================================
# 7. АННОТАЦИЯ SNPEFF
# ==============================================================================
bcftools annotate --rename-chrs <(printf 'NC_002516.2\tChromosome\n') \
    -Oz -o "${PREFIX}.snpeff_input.vcf.gz" "$FILTERED_VCF"
conda activate "$ENV_SNPEFF"
echo "-> Аннотация SnpEff: ${SNPEFF_DB}"
snpEff ann \
    -stats "$REPORT_DIR/snpeff.html" \
    "$SNPEFF_DB" "${PREFIX}.snpeff_input.vcf.gz" \
    > "$ANNOTATED_RAW" \
    2> "$REPORT_DIR/snpeff.log"
conda activate "$ENV_BCF"

bgzip -c "$ANNOTATED_RAW" > "$ANNOTATED_VCF"
bcftools index -f "$ANNOTATED_VCF"
# ==============================================================================
# 8. РАЗДЕЛЕНИЕ ЗАПИСЕЙ С АННОТАЦИЕЙ И БЕЗ
# ==============================================================================
echo "-> Разделение по наличию INFO/ANN"
bcftools view -i 'INFO/ANN!="."' "$ANNOTATED_VCF" \
    -Oz -o "$WITH_ANN_VCF"

bcftools view -i 'INFO/ANN="."' "$ANNOTATED_VCF" \
    -Oz -o "$WITHOUT_ANN_VCF"

bcftools index -f "$WITH_ANN_VCF"
bcftools index -f "$WITHOUT_ANN_VCF"
# ==============================================================================
# 9. ПОДСЧЁТ ЗАПИСЕЙ
# ==============================================================================
N_BEFORE=$(bcftools view -H "$MERGED_VCF" | wc -l)
N_FILTERED=$(bcftools view -H "$FILTERED_VCF" | wc -l)
N_WITH_ANN=$(bcftools view -H "$WITH_ANN_VCF" | wc -l)
N_WITHOUT_ANN=$(bcftools view -H "$WITHOUT_ANN_VCF" | wc -l)

{
    printf 'Категория\tКоличество записей\n'
    printf 'До фильтрации\t%s\n' "$N_BEFORE"
    printf 'После фильтрации\t%s\n' "$N_FILTERED"
    printf 'С аннотацией ANN\t%s\n' "$N_WITH_ANN"
    printf 'Без аннотации ANN\t%s\n' "$N_WITHOUT_ANN"
} > "$REPORT_DIR/variant_counts.tsv"

cat "$REPORT_DIR/variant_counts.tsv"
# ==============================================================================
# 10. СТАТИСТИКА И ГРАФИКИ
# ==============================================================================
echo "-> Статистика до и после фильтрации"
conda activate "$ENV_BCF"
bcftools stats -s - "$MERGED_VCF" \
    > "$REPORT_DIR/before_filtering.stats"

bcftools stats -s - "$FILTERED_VCF" \
    > "$REPORT_DIR/after_filtering.stats"

echo "-> Построение графиков"

plot-vcfstats -P \
    -p "$REPORT_DIR/plots_before/" \
    "$REPORT_DIR/before_filtering.stats"

plot-vcfstats -P \
    -p "$REPORT_DIR/plots_after/" \
    "$REPORT_DIR/after_filtering.stats"

echo "-> Готово. Отчёты: $REPORT_DIR"