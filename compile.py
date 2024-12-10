###압축###
import shutil                                           
                                                        
# 압축할 디렉토리와 출력 파일 경로 지정                 
directory_to_compress = "/home/soowon/project/02_CACHE/CACHE3_4/sample_CACHE_model"

# 압축할 디렉토리 경로                                                   
output_file = "/home/soowon/project/02_CACHE/CACHE3_4/main_trainingdata"  

# 결과 파일 경로 (확장자 제외)                                                
                                                        
# tar.gz 형식으로 압축                                  
shutil.make_archive(output_file, 'gztar', directory_to_compress)                              
                                                        
print(f"압축 완료: {output_file}.tar.gz")   
