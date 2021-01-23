import numpy as np
from sklearn.linear_model import LinearRegression

def getDeliveries():
        X= []
        Y= []
        f  = open("deliveries.txt")
        s = f.readlines()
        for line in s:
                line = line.strip()
                g = ""
                for (_, c) in enumerate(line):
                        if c.isdigit():
                                g+=c
                        else:
                                g+=" "
                nums = list(map(int, g.split()))
                X.append(nums[:4])
                Y.append(nums[4])

        f.close()
        print(X)
        print(Y)
        x, y = np.array(X), np.array(Y)
        model = LinearRegression().fit(x, y)
        print('slope:', model.coef_)



getDeliveries()