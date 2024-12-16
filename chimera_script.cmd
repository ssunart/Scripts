###use chimera without GUI
$chimera --nogui

###use script for autoprocess
$chimera --nogui script.cmd


-----script example(script.cmd)----

open 4LAE.pdb   #현재 위치에 있는 4LAE.pdb 파일 열기
select          #전체 선택
~select ::1VM   #'1VM' 분자만 선택해제
delete sel      # 선택한 구조 제거
select :.X      # chain X 선택
select invert   # 선택 반전(chain X 빼고 전부 선택)
delete sel      # 선택한 구조 제거
addh            # 수소 추가
addcharge std   # standard residue charge 할당
addcharge nonstd #0 0 method gas      #0번 model nonstandard residue gasteiger charge 할당
write format mol2 #0 4LAE_Nlig.mol2  #0번 model mol2 format으로 '4LAE_Nlig.mol2' 이름으로 저장

----------------------

<Other commands>
delete :/isHet  ###HETATM delete


##Command test##
1. Chimera를 GUI로 열어 상단의 'Favorites-Command Line'을 클릭하면 화면 하단에 Command 활성화
2. Test할 명령어 입력 후 제대로 작동하는 지 확인



##아래 링크에서 자주 쓰는 명령어들 찾아서 정리중##
https://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/framecommand.html
