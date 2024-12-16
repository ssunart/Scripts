import numpy as np

def get_center_coordinates_mol2_direct(filename):
    coordinates = []
    read_coordinates = False

    with open(filename, 'r') as f:
        for line in f:
            # @<TRIPOS>ATOM 섹션을 찾음
            if line.startswith("@<TRIPOS>ATOM"):
                read_coordinates = True
                continue
            # @<TRIPOS>BOND 섹션을 찾으면 좌표 읽기를 중단
            if line.startswith("@<TRIPOS>BOND"):
                read_coordinates = False

            # 좌표 읽기
            if read_coordinates:
                parts = line.split()
                if len(parts) < 6:  # 유효한 줄인지 확인
                    continue
                try:
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                except ValueError:  # 숫자로 변환할 수 없는 값이면 건너뜀
                    continue
                coordinates.append([x, y, z])

    # numpy 배열로 변환 후 중심 좌표를 계산
    if len(coordinates) == 0:
        return "No coordinates found"

    coordinates = np.array(coordinates)
    center = np.mean(coordinates, axis=0)
    return center

mol2_file = '1G6S_Nlig.mol2'
center_coordinates = get_center_coordinates_mol2_direct(mol2_file)

# 결과를 nativeLig_center.txt에 저장
with open('nativeLig_center.txt', 'w') as f:
    if isinstance(center_coordinates, str):
        f.write(center_coordinates + '\n')
    else:
        f.write(f"The center coordinates are: {center_coordinates}\n")

print(f"The center coordinates are: {center_coordinates}")
