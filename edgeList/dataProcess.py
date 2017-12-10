#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 20:26:32 2017

@author: jingqiu
"""

import numpy as np
import matplotlib.pyplot as plt
time=np.loadtxt("time4P.txt");
time5=time[:,1]
print(np.median(time5))
plt.hist(time5)
time=np.loadtxt("time1P.txt");
time5=time[:,1]
plt.hist(time5)
print(np.median(time5))#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 20:26:32 2017

@author: jingqiu
"""

'''import numpy as np
import matplotlib.pyplot as plt
time=np.loadtxt("time4P.txt");
time_sum=np.sum(time[:,1:2],axis=1);
print(np.median(time_sum))
plt.hist(time_sum)
time=np.loadtxt("time1P.txt");
time_sum=np.sum(time[:,1:2],axis=1);
plt.hist(time_sum)
print(np.median(time_sum))'''



    
    



    
    