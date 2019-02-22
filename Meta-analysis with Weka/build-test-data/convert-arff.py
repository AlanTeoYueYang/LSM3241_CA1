import csv
import os

EMT = [1226581,1226582,1226583]
Normal = [1226584,1226585,1226586]

arff = open("test_data.arff","w+")
arff.write('@relation genes_expression\n\n')

attributes_w = False

attributes = {}
test = True
for file in os.listdir():
	if file.endswith(".csv"):
		with open(file) as csv_file:
			csv_reader = csv.reader(csv_file)
			header = True
			for line in csv_reader:
				if header:
					samples = list(map(lambda x: int(x[3:-4]) if x.endswith('.CEL') else int(x[3:]), line[1:]))
					data = [''] * len(samples)
					header = False
				else:
					if not attributes_w:
						attribute = line[0]
						arff.write('@attribute ' + attribute + ' numeric\n')
					temp_data = line[1:]
					for i in range(len(data)):
						data[i] += (temp_data[i] + ',')
			for j in range(len(samples)):
				if samples[j] in Normal:
					data[j] += ("Normal")
				elif samples[j] in EMT:
					data[j] += ("EMT")
				else:
					print("Error: sample is not in EMT/Normal samples")
			if not attributes_w:
				arff.write('@attribute subtype{Normal,EMT}\n\n@data\n')
			attributes_w = True
			for d in data:
				arff.write(d + '\n')
