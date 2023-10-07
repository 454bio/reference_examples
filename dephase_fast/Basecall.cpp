#include <stdio.h>
#include "Model4.h"


struct PhaseParams {
    double ie;
    double cf;
    double dr;
    double err;
};

void template2bases(char *dnaTemplate, char *basecalls)
{
    char bases[] = {'G', 'C', 'A', 'T'};
    int len = strlen(dnaTemplate);
    int i = 0;
    for(i=0;i<len;i++)
        basecalls[i] = bases[dnaTemplate[i]-1];
    basecalls[i] = 0;
}

std::vector<Signal> LoadCycleIntensities(const char *filename)
{
    std::vector<Signal> vals;
    Signal s;

    FILE *fp = fopen(filename, "r");
    if (fp) {
        char line[1024];
        while (fgets(line, sizeof(line), fp)) {
            sscanf(line, "%lf,%lf,%lf,%lf", &s.v[0], &s.v[1], &s.v[2], &s.v[3]);
            s.v[0] /= 100.0;
            s.v[1] /= 100.0;
            s.v[2] /= 100.0;
            s.v[3] /= 100.0;
            vals.push_back(s);
        }
        fclose(fp);
    }

    return vals;
}

void LoadSpotData(const char *filename, std::vector<std::vector<Signal>> &spotData)
{
    std::vector<Signal> vals;
    int spotId;
    char name[256];
    int cycle;
    Signal s;
    int lastSpotId = 0;

    FILE *fp = fopen(filename, "r");
    if (fp) {
        char line[1024];
        // read and ignore the header
        fgets(line, sizeof(line), fp);
        while (fgets(line, sizeof(line), fp)) {
            // skip to the 4 values
            char *ptr = line;
            for(int i=0;i<3;i++) {
                // find a comma, then skip over it
                if (ptr)
                    ptr = strchr(ptr, ',');
                if (ptr)
                    ptr++;
            }
            sscanf(ptr, "%lf,%lf,%lf,%lf", &s.v[0], &s.v[1], &s.v[2], &s.v[3]);
            sscanf(line, "%d", &spotId);
            if (spotId != lastSpotId) {
                if (lastSpotId != 0) {
                    std::vector<Signal> v;
                    v = vals;
                    spotData.push_back(v);
                }
                lastSpotId = spotId;
                vals.clear();
            }
            vals.push_back(s);
        }
        fclose(fp);
        spotData.push_back(vals);
    }
}

double CallBases(char *dnaTemplate, std::vector<Signal> &measuredSignal, double ie, double cf, double dr)
{
    int numCycles = measuredSignal.size();

    Model4 m;
    m.Init(100);
    m.SetParams(ie, cf, dr);

    std::vector<Signal> dyeIntensities(numCycles);
    std::vector<double> totalSignal(numCycles, 0);
    std::vector<double> errorPerCycle(numCycles, 0);

    memset(dnaTemplate, 0, 1024); // TODO really need to include size
    char bases[4] = {'A', 'C', 'G', 'T'};
    double cumulativeError;

    int numIterations = 3;
    for(int iteration=0;iteration<numIterations;iteration++) {
        m.Reset();
        cumulativeError = 0.0;
        for(int cycle=0;cycle<numCycles;cycle++) {
            double best_error = 0.0;
            int best_base = -1;
            Signal best_signal;
            for(int base=0;base<4;base++) {
                // insert a "what-if" base at the current position, and predict what the signal looks like
                dnaTemplate[cycle] = base+1;
                Signal signal = m.GetSignal(dnaTemplate);
                double signalSum = 0.0;
                for(int i=0;i<4;i++)
                    signalSum += signal.v[i];

                double error = 0.0;
                for(int i=0;i<4;i++) {
                    double delta = (measuredSignal[cycle].v[i] - signal.v[i])/signalSum;
                    error += delta*delta;
                }

                // keep track of the lowest error, this is the best predition
                if (error < best_error || best_base == -1) {
                    best_base = base;
                    best_error = error;
                    best_signal = signal;
                }
            }

            // append/replace with best base at current position (cycle)
            dnaTemplate[cycle] = best_base+1;
            for(int i=0;i<5;i++)
                dyeIntensities[cycle].v[i] = best_signal.v[i];
            totalSignal[cycle] = 0.0;
            for(int i=0;i<4;i++)
                totalSignal[cycle] += best_signal.v[i];
            errorPerCycle[cycle] = best_error;

            // update the model - note that we do this after getting the measured signals, because this matches the physical
            // system where the first base is exposed to nucleotides prior to UV cleavage
            m.ApplyUV(dnaTemplate, numCycles);

            cumulativeError += errorPerCycle[cycle];
        }
    }

    return cumulativeError;
}

PhaseParams GridSearch(char *dnaTemplate, std::vector<Signal> &measuredSignal)
{
    // grid-search first, then call based on lowest error
    double drmin = 0.01;
    double drmax = 0.03;
    int    drnum = 10;
    double iemin = 0.07;
    double iemax = 0.12;
    int    ienum = 20;
    double cfmin = 0.07;
    double cfmax = 0.12;
    int    cfnum = 20;

    int numCycles = measuredSignal.size();
    PhaseParams params;
    double minerr = 99999.0;

    for(int dri=0;dri<drnum;dri++) {
        double drtest = drmin + (dri/(double)(drnum-1.0)) * (drmax-drmin);
        for(int cfi=0;cfi<cfnum;cfi++) {
            double cftest = cfmin + (cfi/(double)(cfnum-1.0)) * (cfmax-cfmin);
            for(int iei=0;iei<ienum;iei++) {
                double ietest = iemin + (iei/(double)(ienum-1.0)) * (iemax-iemin);
                double err = CallBases(dnaTemplate, measuredSignal, ietest, cftest, drtest);
                if (err < minerr) {
                    minerr = err;
                    params.ie = ietest;
                    params.cf = cftest;
                    params.dr = drtest;
                }
            }
        }
    }
    params.err = CallBases(dnaTemplate, measuredSignal, params.ie, params.cf, params.dr);
    return params;
}

int main(int argc, char *argv[])
{
    // load spots
    // std::vector<Signal> measuredSignal = LoadCycleIntensities("jmrdata1.csv");
    std::vector<std::vector<Signal>> spotData;
    LoadSpotData("color_transformed_spots.csv", spotData);
    int numSpots = spotData.size();
    printf("Loaded %d spots\n", numSpots);

    // grid-search phase-correct each spot
    char dnaTemplate[1024];
    char basecalls[1024];
    for(int i=0;i<numSpots;i++) {
        PhaseParams params = GridSearch(dnaTemplate, spotData[i]);
        template2bases(dnaTemplate, basecalls);
        printf("spot: %d ie/cf/dr: %.3lf/%.3lf/%.3lf err: %.3lf template: %s\n", i, params.ie, params.cf, params.dr, params.err, basecalls);
    }

    // generate some interesting csv outputs
}

