import os

# 설정
#target_data_path = "/home/soowon/project/07_PPAR/docking/6c5t_PPAR/6c5t_prep.maps.fld"
target_data_path = "/home/soowon/project/07_PPAR/docking/6c5t_PPAR/6c5t_f/6c5t_prep.maps.fld"
ligand_dir = "/data1/DB/Molecule/Enamine/Screening/PDBQT/39/"
output_dir = "/home/soowon/project/07_PPAR/docking/6c5t_PPAR/gpuR/39/"

# 모든 pdbqt 파일을 찾아 job file 생성
with open("job_file39_f", "w") as job_file:
    for ligand in os.listdir(ligand_dir):
        if ligand.endswith(".pdbqt"):
            # 각 파일에 대한 경로 생성
            ligand_path = os.path.join(ligand_dir, ligand)
            output_prefix = ligand.split('.')[0]
            result_file = os.path.join(output_dir, output_prefix + ".dlg")  # 결과 파일 경로 확인용

            # 결과 파일(.dlg)이 이미 존재하는지 확인
            if not os.path.exists(result_file):
                # job file에 내용 작성 (결과 파일이 없는 경우에만)
                job_file.write(f"{target_data_path}\n{ligand_path}\n{output_dir}{output_prefix}\n")

print("Job file 생성 완료. 중복 계산을 방지했습니다.")
