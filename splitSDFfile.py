###다중 리간드 SDF 파일을 개별 리간드 SDF 파일로 분리

import os, sys

def readSDFfile(sdf_f):
    return open(sdf_f).readlines()


def splitMultiSDFtoSingleSDF(sdf_f):
    multi_lig = readSDFfile(sdf_f)
    idxList = []
    name = {}
    data = []
    chk = []
    for idx, i in enumerate(multi_lig):
        if i.startswith('$$$$'):
            idxList.append(idx+1)

    
    for idx, i in enumerate(idxList):
        if idx == 0:
            f = open(f"1onp_Lig_{idx}.sdf", "w")
            for j in multi_lig[0:i]:
                f.write(f"{j}")
        else:
            f = open(f"1onp_Lig_{idx}.sdf", "w")
            for j in multi_lig[idxList[idx-1]:i]:
                f.write(f"{j}")



if __name__ == '__main__':
    splitMultiSDFtoSingleSDF('1onp_Ligandsource.sdf')
