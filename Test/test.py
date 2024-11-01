import subprocess
import filecmp
import os
def Compile_LF_LDG():
	subprocess.call('make -f Makefile_lf_ldg', shell=True)

def Compile_ROE_BR2():
	subprocess.call('make -f Makefile_roe_br2', shell=True)

def SolvePDE(num):
	# solve_call = 'demo_solve min_p'+str(num)+'.pde'
	solve_call = 'ngs min_p'+str(num)+'.pde'
	subprocess.call(solve_call, shell=True)
def TestResult(golden_file, compare_file, num):
	if(filecmp.cmp(golden_file, compare_file)):
		os.remove(compare_file)
		return 1
	else:
		return 0

TotalSize = 12
TestResultsList = TotalSize * [0] 

Compile_LF_LDG()
# Running test 4
golden_file = 'golden_solution_p4.out'
compare_file = 'solution-32-3.out'
num = 4
SolvePDE(num)
val = TestResult(golden_file, compare_file, num)
TestResultsList[num-1] = val

# Running test 6
golden_file = 'golden_solution_p6.out'
compare_file = 'solution-32-5.out'
num = 6
SolvePDE(num)
val = TestResult(golden_file, compare_file, num)
TestResultsList[num-1] = val

# # Running test 10
golden_file = 'golden_solution_p10.out'
compare_file = 'solution-96-3.out'
num = 10
SolvePDE(num)
val = TestResult(golden_file, compare_file, num)
TestResultsList[num-1] = val

# Running test 12
golden_file = 'golden_solution_p12.out'
compare_file = 'solution-96-5.out'
num = 12
SolvePDE(num)
val = TestResult(golden_file, compare_file, num)
TestResultsList[num-1] = val

Compile_ROE_BR2()
# Running test 3
golden_file = 'golden_solution_p3.out'
compare_file = 'solution-32-2.out'
num = 3
SolvePDE(num)
val = TestResult(golden_file, compare_file, num)
TestResultsList[num-1] = val

# Running test 9
golden_file = 'golden_solution_p9.out'
compare_file = 'solution-96-2.out'
num = 9
SolvePDE(num)
val = TestResult(golden_file, compare_file, num)
TestResultsList[num-1] = val

# print(TestResultsList)

# for i in range(0, len(TestResultsList)):
# 	if(TestResultsList[i]):
# 		print ("Passed test "+str(i+1))
# 	else:
# 		print ("Failed test "+str(i+1))