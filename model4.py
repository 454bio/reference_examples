# copyright 2023 - Mel Davey

import numpy as np

#
# The 'Model' class approximates the physical system of DNA strand incorporation and UV cleavage
# here, we account for the following per dye:
# incomplete extention (ie) where not all available binding sites incorporate a nucleotide 
# carry-forward (cf) where further incorporation can occur after the tagged nucleotide is cleaved
# droop (dr) where some DNA strands can no longer support new incorporations
#

class Model:
    def __init__(self, strandLen:int = 100):
        self.strandLen = strandLen

        self.state = np.zeros(strandLen)
        self.prevState = np.zeros(strandLen)

        self.state[0] = 1.0
        self.dr = 0.01
        self.ie = np.zeros(4)
        self.ie[:] = 0.05
        self.cf = np.zeros(4)
        self.cf[:] = 0.07
        self.cfDefault = np.mean(self.cf)
        self.extraBucket = 0
        self.bases = ['A', 'C', 'G', 'T']

    def SetParams(self, ie, cf, dr:float = 0):
        self.ie = ie
        self.cf = cf
        self.cfDefault = np.mean(self.cf)
        self.dr = dr

    def ApplyUV(self, dnaTemplate:str, maxLen:int = 0):
        # save current state
        self.prevState = self.state[:]
        self.extraBucket = 0
        numExtensions = 3 # technically this goes on forever, but after 3 rounds there is not much left
        extendAmount = np.zeros(numExtensions)
        totalChange = np.zeros(self.strandLen+numExtensions)
        templateLen = len(dnaTemplate)

        # apply incompletion and carry-forward effects to our state
        for i in range(0, self.strandLen):
            # amount of product available to extend
            extendAmount[0] = self.state[i] * (1.0 - self.ie[0])
            changeAmount = extendAmount[0]
            if extendAmount[0] == 0 or (maxLen > 0 and i >= (maxLen-1)):
                continue

            # continue to extend (carry-forward)
            for s in range(1,numExtensions):
                templateIndex = i + s - 1
                if templateIndex < templateLen:
                    cf = self.cf[self.bases.index(dnaTemplate[templateIndex])]
                else:
                    cf = self.cfDefault
                cfAmount = extendAmount[s-1] * cf
                extendAmount[s] = cfAmount
                extendAmount[s-1] -= extendAmount[s]

            # update the total change from this strand's position
            totalChange[i+1:i+1+numExtensions] += extendAmount
            totalChange[i] -= changeAmount

        # apply the change to the current state
        self.state += totalChange[:self.strandLen]

        # droop is applied across all states
        self.state *= (1.0 - self.dr)


    def GetSignal(self, dnaTemplate):
        signal = np.zeros(6) # 4 bases plus unknown & extra
        templateLen = len(dnaTemplate)
        # each position within our state, sum the signal, binned by known DNA bases
        for i in range(self.strandLen):
            if i < templateLen:
                signal[self.bases.index(dnaTemplate[i])] += self.state[i]
            else:
                signal[4] += self.state[i]
        signal[5] = self.extraBucket
        return signal

    def Revert(self):
        # revert to previous state
        self.state = self.prevState[:]

    def GetState(self):
        return self.state

    def Reset(self):
        self.state.fill(0)
        self.state[0] = 1.0

