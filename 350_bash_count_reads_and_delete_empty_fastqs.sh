













# as fastest solution taken from:
# https://bioinformatics.stackexchange.com/questions/935/fast-way-to-count-
#   number-of-reads-and-number-of-bases-in-a-fastq-file

fix_base_count() {
    local counts=($(cat))
    echo "${counts[0]} $((${counts[1]} - ${counts[0]}))"
}

gzip -dc ERR047740_1.filt.fastq.gz \
    | awk 'NR % 4 == 2' \
    | wc -cl \
    | fix_base_count
