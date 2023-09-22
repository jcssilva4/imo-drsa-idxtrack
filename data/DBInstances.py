import math
import numpy as np
import os

class Instance: 
    def __init__(self, instance_n, sourceName):
        if sourceName.lower() == "beasley":
            self.readBeasleyData(instance_n)

    def readBeasleyData(self, instance_n):
        '''
        Data Description:
        nStocks nWeeks (nWeeks = nPrices)
        IndexPrice1 IndexPrice2 ... 
        IndexPrice_nWeeks
        Stock1Price1 Stock1Price2...
        Stock1Price_nWeeks
        .
        .
        .
        Stock_nStocksPrice1 Stock_nStocksPrice2
        Stock_nStocksPrice_nWeeks
        '''
        print('reading Beasley_indtrack' + str(instance_n) + " instance data")

        # get file object
        filelines = []
        filehandle = open("data/indtrack" + str(instance_n)+ ".txt" , "r")
        filelines.extend(filehandle.readlines())
        allData = []

        # get prices
        for line in filelines:
            allData.extend(line.split(" "))
        allData = [float(i) for i in allData if i != '\n' and i != ''] # filter some elements
        n = int(allData[0]) # number of stocks
        T = int(allData[1]) # nWeeks
        prices_Index = allData[2:2 + (T+1)]

        # get index's returns
        self.returns_Index = [
            (prices_Index[t] - prices_Index[t-1])/prices_Index[t-1] for t in range(1, T + 1)]

        '''
        # get stocks' returns
        self.returns_i = []
        for stock in range(n):
            low_ = 2 + (T+1)*(stock + 1) #T + 1 because we have 291 prices (290 returns)
            upper_ = low_ + (T+1) #T+1 because we have 291 prices (290 returns)
            prices_t = allData[low_:upper_]
            self.returns_i.append(
                [(prices_t[t] - prices_t[t-1])/prices_t[t-1] for t in range(1, T + 1)]
            )
        '''
        # Here the return are calculate as Beasley2003 returns (R_t = log_e(P_t/P_(t-1)))
        self.returns_i = []
        for stock in range(n):
            low_ = 2 + (T+1)*(stock + 1) # T + 1 because we have 291 prices (290 returns)
            upper_ = low_ + (T+1) # T+1 because we have 291 prices (290 returns)
            prices_t = allData[low_:upper_]
            self.returns_i.append(
                [np.log(prices_t[t]/prices_t[t-1]) for t in range(1, T + 1)]
            )        

        filehandle.close()

 