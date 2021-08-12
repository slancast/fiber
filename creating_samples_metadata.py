import os 
import numpy as np
from shutil import copyfile
import shutil

cwd=os.getcwd()
counter = 0
array = [["sample","participant","fiber","week","run"]]

for root, dirs, files in os.walk(cwd):
	for dir in dirs:
		shutil.rmtree(dir)

for root, dirs, files in os.walk(cwd):
	print(files)
	try:
		files.remove("samples.txt")
		files.remove(".DS_Store")
	except ValueError:
		print("Error")
	for file in files:
		print(file)
		z=file.split("_")
		print(z)
		if z[1] == "1" or z[1] == "2":
			z[2] = z[1] + z[2]
			del z[1]
		if z[1] == "LC" or z[1] == "SC" or z[1] == "1LC" or z[1] == "2LC" or z[1] == "1SC" or z[1] == "2SC":
			z[1] = z[1] + z[2]
			del z[2]
		if z[2] == "Washout":
			z[2] = z[2] + z[3]
			del z[3]
		if z[2] == "WashoutWeek" or z[2] == "Washoutweek" or z[2] == "Washoutweek":
			z[2] = z[2] + z[3]
			del z[3]
		if z[2] == "End":
			del z[3]
		if z[2] == "":
			z[2] = "?"
		print(z)
		if z[3] == "rsem.genes.results.txt":
			counter = counter + 1
			a = [counter]
			a.append(z[0])
			a.append(z[1])
			a.append(z[2])
			b=z[0]+"_"+z[1]+"_"+z[2]
			a.append(b)
			if a[2].startswith("1") or a[2].startswith("2"):
				a[2] = a[2][1:]
			a.append(z[3])
			array.append(a)
		print(b)
		try:
			os.mkdir(cwd+"/"+b)
			break
		except OSError:
			print(cwd+"/"+b+" already made")
			break
		copyfile(cwd+"/"+file, cwd+"/"+b+"/"+b+".genes.results")
	
array = np.array(array)
with open('samples.txt','wb') as f:
    np.savetxt(f,array,fmt="%s",delimiter="\t")