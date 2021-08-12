import os 
import numpy as np
from shutil import copyfile
import shutil

cwd=os.getcwd()
counter = 0
array = [["sample","participant","fiber","week","run","batch"]]

for root, dirs, files in os.walk(cwd):
	for dir in dirs:
		shutil.rmtree(dir)

counter= 0 
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
		z[2] = z[2].replace(" ","")
		z[4] = z[4].replace("\r","")
		print(z)
		if z[4] == "rsem.genes.results.txt":
			counter = counter + 1
			a = [counter]
			a.append(z[0])
			a.append(z[1])
			a.append(z[2])
			b=z[0]+"_"+z[1]+"_"+z[2]+"_"+z[3]
			a.append(b)
			if a[2].startswith("1") or a[2].startswith("2"):
				a[2] = a[2][1:]
			a.append(z[3])
			array.append(a)
		
		try:
			os.mkdir(cwd+"/"+b)
		except OSError:
			print(cwd+"/"+b+" already made")
		copyfile(cwd+"/"+file, cwd+"/"+b+"/"+b+".genes.results")

print(counter)
	
array = np.array(array)
with open('samples.txt','wb') as f:
    np.savetxt(f,array,fmt="%s",delimiter="\t")