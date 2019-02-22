import csv
import os

normal = [315430,315431,315432,359925,359926,368860,368863,368866,368867,451814,451815,451816,451820,451821,451822,1063981]
EMT = [315443,315444,315445,359927,359929,359930,359931,368868,368869,451817,451818,451819,1063982,1063983,1063984]

arff = open("train_data.arff","w+")
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
				if samples[j] in normal:
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
