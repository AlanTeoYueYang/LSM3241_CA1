import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.integrate as integrate

f = open('SMO_classification_train_details.txt','r')

lst = []
for line in f.readlines():
	if line.endswith('at\n'):
		temp = line.split(" ")
		for i in temp:
			try:
				weight = float(i)
				if weight < 0:
					weight *= -1
				lst.append(weight)
				break
			except:
				continue

lam = 1/np.mean(lst)
f = lambda x:lam*math.exp(-lam*x)
plt.figure(figsize=(10,5))

# plot histogram 
ax1=plt.subplot(1, 2, 1)
ax1.hist(lst)
ax1.set_xlabel('Absolute Feature Weights')
ax1.set_ylabel('Number of Features')

# plot exponential
ax2=plt.subplot(1, 2, 2)
x = np.arange(0.0, 0.01, 0.0001)
y = np.asarray(list(map(lambda x: f(x),x)))
ax2.plot(x,y)
ax2.set_xlabel('x')
ax2.set_ylabel('Exp(x)')
plt.show()

# find p-value
print(integrate.quad(lambda x: f(x), 0.002,1)[0])
