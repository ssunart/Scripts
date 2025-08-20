-P 위치에 " https"링크의 파일 다운로드 &  download_smallbfd.log 에 로그 기록 & 백그라운드 작동
nohup wget -P ../alphafold_data/small_bfd/ "https://storage.googleapis.com/alphafold-databases/reduced_dbs/bfd-first_non_consensus_sequences.fasta.gz" > download_smallbfd.log 2>&1 &
