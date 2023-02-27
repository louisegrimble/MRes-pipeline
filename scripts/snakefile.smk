(CLUSTER,) = glob_wildcards("data/{CLUSTER}.fasta")

rule all:
  input:
    expand("results/prokka_out_test/{CLUSTER}/{CLUSTER}.faa", CLUSTER = CLUSTER),
    expand("results/eggnog_out_test/{CLUSTER}.emapper.annotations", CLUSTER = CLUSTER),
    expand("results/eggnog_out_test/{CLUSTER}.in.emapper.annotations", CLUSTER = CLUSTER),
    expand("results/kofamscan_out_test/kofam_out_{CLUSTER}.tsv", CLUSTER = CLUSTER),
    expand("results/kofamscan_out_test/kofam_out_{CLUSTER}.txt", CLUSTER = CLUSTER),
    expand("results/eggnog_out_test/r2_egg_{CLUSTER}", CLUSTER = CLUSTER),
    expand("results/eggnog_out_test/r2_egg_{CLUSTER}.txt", CLUSTER = CLUSTER),
    expand("results/eggnog_out_test/head_r2_egg_{CLUSTER}.txt", CLUSTER = CLUSTER),
    expand("results/eggnog_out_test/{CLUSTER}_egg_for_decoder.txt", CLUSTER = CLUSTER),
    expand("results/kofamscan_out_test/kofam_decoder_{CLUSTER}.txt", CLUSTER = CLUSTER),
    expand("results/decoder_out_test/concat_kofam_{CLUSTER}.txt", CLUSTER = CLUSTER),
    expand("results/decoder_out_test/decoder_out_{CLUSTER}", CLUSTER = CLUSTER),
    expand("results/python_out_{CLUSTER}/stats_df", CLUSTER = CLUSTER),
    expand("results/python_out_{CLUSTER}/low_qual_df", CLUSTER = CLUSTER),
    expand("results/krakren_test/kraken2_{CLUSTER}_out", CLUSTER = CLUSTER),
    expand("results/interproscan_test/{CLUSTER}_interpro_out", CLUSTER = CLUSTER),
    expand("results/blast_out_test/blast_out_{CLUSTER}", CLUSTER = CLUSTER),
    expand("results/seqkit_test/{CLUSTER}_seqkit_length.txt", CLUSTER = CLUSTER)

rule prokka:
  output: "results/prokka_out_test/{CLUSTER}/{CLUSTER}.faa"
  input:  "data/{CLUSTER}.fasta"
  params: dir= "results/prokka_out_test/{CLUSTER}",
          prefix= CLUSTER
  conda: "envs/prokka_env.yaml"
  shell:
    "prokka {input} --outdir {params.dir} --prefix {params.prefix} --force"

rule eggnog:
  output: "results/eggnog_out_test/{CLUSTER}.emapper.annotations"
  input:  expand("results/prokka_out_test/{{CLUSTER}}/{{CLUSTER}}.faa")
  params: dir1= "results/eggnog_out_test",
          name= "{CLUSTER}"
  conda: "envs/pip_egg.yaml"
  shell:
      """
      mkdir -p {params.dir1}
      emapper.py --dmnd_ignore_warnings --no_file_comments --data_dir /mnt/lustre/groups/biol-chong-2019/databases/eggnog/ --itype proteins -i {input} --output_dir {params.dir1} -o {params.name}
      """

rule egg_zero:
  output: "results/eggnog_out_test/{CLUSTER}.in.emapper.annotations"
  input: "results/eggnog_out_test/{CLUSTER}.emapper.annotations"
  shell:
    "sed -i 's/\t-/\t0/g' {input} > {output}"

rule kofam_scan:
  output:
    tsv = "results/kofamscan_out_test/{CLUSTER}.tsv",
    txt = "results/kofamscan_out_test/{CLUSTER}.txt"
  input: expand("results/prokka_out_test/{{CLUSTER}}/{{CLUSTER}}.faa")
  conda: "envs/kofamscan.yaml"
  shell:
     """
     ../kofam_test/kofam_scan-1.3.0/exec_annotation -f detail-tsv -o {output.tsv} {input} --profile=../kofam_test/profiles
     ../kofam_test/kofam_scan-1.3.0/exec_annotation -f mapper -o {output.txt} {input} --profile=../kofam_test/profiles
     """

rule egg_for_decoder:
  output:  "results/eggnog_out_test/{CLUSTER}_regex_egg"
  input:  "results/eggnog_out_test/{CLUSTER}.in.emapper.annotations"
  script: "scripts/egg_regex.py"

rule egg_format_decoder:
  output:
    one = "results/eggnog_out_test/{CLUSTER}_egg_1.txt",
    two = "results/eggnog_out_test/{CLUSTER}_egg_2.txt"
  input: "results/eggnog_out_test/{CLUSTER}_regex_egg"
  shell:
    """
    sed 's/,/\t/g' r2_egg_{CLUSTER} > {output.one} #removes commas
    awk 'NR!= 1' r2_egg_{CLUSTER}.txt > {output.two} #removes header
    """

rule rename_egg_decoder:
  output:  "results/eggnog_out_test/{CLUSTER}_egg_for_decoder.txt"
  input:  "results/eggnog_out_test/{CLUSTER}_egg_2.txt"
  shell: "sed 's/^/egg_/g' {input} > {output}"

rule rename_kofam_decoder:
  output:  "results/kofamscan_out_test/{CLUSTER}_kofam.txt"
  input:  "results/kofamscan_out_test/{CLUSTER}.txt"
  shell: "sed 's/^/kofam_/g'  {input} > {output}"

rule concat_for_decoder:
  output:  "results/decoder_out_test/{CLUSTER}_concat.txt"
  input:
    egg = "results/eggnog_out_test/{CLUSTER}_egg_for_decoder.txt",
    kofam = "results/kofamscan_out_test/{CLUSTER}_kofam.txt"
  shell: "cat {input.egg} {input.kofam} > {output}.txt"

rule keggdecoder:
  output:  "results/decoder_out_test/{CLUSTER}_decoder_out"
  input:  "results/decoder_out_test/{CLUSTER}_concat.txt"
  conda: "envs/keggdecoder.yaml"
  shell:
    "KEGG-decoder --input {input} --output {output} --vizoption static"

rule python_info:
  output: "results/python_out_{CLUSTER}/stats_df",
          "results/python_out_{CLUSTER}/egg_tax_df",
          "results/python_out_{CLUSTER}/egg_module_df",
          "results/python_out_{CLUSTER}/egg_cog_df",
          "results/python_out_{CLUSTER}/subset",
          "results/python_out_{CLUSTER}/low_qual_df"
  input: "results/eggnog_out_test/{CLUSTER}.emapper.annotations",
         "results/kofamscan_out_test/kofam_out_{CLUSTER}.tsv"
  script: "scripts/stats_python.py"

rule low_qual:
   output: "results/seqkit_test/{CLUSTER}_low_qual.faa"
   input:
    txt = "results/python_out_{CLUSTER}/low_qual_df",
    seq = "results/prokka_out_test/{CLUSTER}/{CLUSTER}.faa"
   conda: "envs/seqkit_env.yaml"
   shell:
     "seqkit grep -f {input.txt} {input.seq} > {output}"

rule kraken2:
  output: "results/krakren_test/{CLUSTER}_kraken2_out"
  input:  "results/seqkit_test/{CLUSTER}_low_qual.faa"
  conda: "envs/kraken_env.yaml"
  shell:
    "kraken2 --db /mnt/lustre/groups/biol-chong-2019/databases/krakendb/kraken2_samstudio8 {input} --output {output}"

rule interproscan:
  output: "results/interproscan_test/{CLUSTER}_interpro_out"
  input:  "results/seqkit_test/{CLUSTER}_low_qual.faa"
  conda: "envs/interproscan_env.yaml"
  shell:
    "../interproscan_test/interproscan-5.59-91.0/interproscan.sh -i {input} -goterms -b {output}"

rule blast_p:
  output: "results/blast_out_test/blast_out_{CLUSTER}"
  input: "results/seqkit_test/{CLUSTER}_low_qual.faa"
  conda: "envs/blast_env.yaml"
  shell:
    """
    module load bio/BLAST+/2.11.0-gompi-2021a
    blastp -query {input} -db /mnt/lustre/groups/biol-database-2020/blastdb/2022-02-21/nr -outfmt 6 -out {output}
    """

rule hypo_list:
  output: "results/seqkit_test/{CLUSTER}_hypothetical_list.txt"
  input:  "results/seqkit_test/{CLUSTER}_low_qual.faa"
  shell:
    "grep hypothetical {input} > {output}"

rule seqkit_length:
  output: "results/seqkit_test/{CLUSTER}_seqkit_length.txt"
  input:  "results/seqkit_test/{CLUSTER}_low_qual.faa"
  conda: "envs/seqkit_env.yaml"
  shell:
    "seqkit fx2tab -l {input} | cut -f 1,4 > {output}"
