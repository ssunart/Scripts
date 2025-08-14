import matplotlib.pyplot as plt

# 데이터 정의
stages = ["Heating", "Equilibration", "Production"]
case1 = [-43906, -43861, -44279]
case2 = [-43274, -44317, -44176]

# 그래프 생성
plt.figure(figsize=(8, 5))
plt.plot(stages, case1, marker='o', label='Case 1 (0K → 300K)')
plt.plot(stages, case2, marker='o', label='Case 2 (0K → 500K → 300K)')

# 그래프 꾸미기
plt.title("MD Simulation Energy Changes", fontsize=14)
plt.xlabel("Simulation Stage", fontsize=12)
plt.ylabel("Energy (kcal/mol)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

plt.show()
