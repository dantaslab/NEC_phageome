#python3
import os

samples=[]

for sample in os.listdir():
	sample = sample.split("_ShortBRED")[0]
	samples.append(sample)

joint_output = []
joint_output.append("Sample\tFamily\tCount\tHits\tTotMarkerLength")

for sample in samples:
	path = os.getcwd() + "/" + sample + "_ShortBRED.txt"
	for line in open(path, "r"):
		if not line.startswith("Family"):
			entry = sample + "\t" + line.strip()
			joint_output.append(entry)

joint_output_file = open(os.getcwd() + "/joint_output.txt", "w")

for line in joint_output:
	print(line, file = joint_output_file)

joint_output_file.close()