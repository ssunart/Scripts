import os

# 현재 디렉토리의 모든 하위 디렉토리 리스트 얻기
current_directory = os.getcwd()
directories = [d for d in os.listdir(current_directory) if os.path.isdir(os.path.join(current_directory, d))]

# 기본 세션 파일 내용
session_template = """
open {dir}/*.pdb {dir}/*.sdf
set bgColor white
select protein & H
del sel
select ::name="UNL"
style sel ball
select sel :< 6
color (#!1 & sel) cyan
color (#!1 & sel) byhetero
show sel atoms
show sel surfaces
transparency (#!1 & sel) 60
label (#!1 & sel) text "{{0.name}}{{0.number}}{{0.insertion_code}}"
ui tool show H-Bonds
hbonds reveal true
ui tool show Clashes
clashes ignoreHiddenModels true
~select
"""

# 디렉토리별 세션 파일 생성
for dir_name in directories:
    session_content = session_template.format(dir=dir_name)
    session_filename = f"{dir_name}.cxc"
    
    with open(session_filename, "w") as session_file:
        session_file.write(session_content)
        
    print(f"Created session file: {session_filename}")
