#include <stdio.h>
#include <stdlib.h>
#include "Model4.h"
#include <string>
#include <map>
#include "EditDist.h"
#include <math.h>

struct PhaseParams {
    double ie;
    double cf;
    double dr;
    double err;
};

struct SpotData {
    std::vector<Signal> vals;
    std::string spotName;
};

std::map<std::string,std::string> spotSeq = {
    {"355", "ACGTGACTAGTGCATCACGTGACTAGTGCATC"},
    {"357", "ATGCAGTCGACGTACTATGCAGTCGACGTACT"},
    {"358", "CGTATCGACTATGCAGCGTATCGACTATGCAG"},
    {"360", "GACTCGATGCTCAGTAGACTCGATGCTCAGTA"},
    {"364", "TCAGTACGATGACTGCTCAGTACGATGACTGC"},
    {"370", "ACGTGACTAGTGCATCACGTGACTAGTGCATC"},
    {"372", "ATGCAGTCGACGTACTATGCAGTCGACGTACT"},
    {"373", "CGTATCGACTATGCAGCGTATCGACTATGCAG"},
    {"375", "GACTCGATGCTCAGTAGACTCGATGCTCAGTA"},
    {"377", "GTCAGCTACGACTGATGTCAGCTACGACTGAT"},
    {"379", "TCAGTACGATGACTGCTCAGTACGATGACTGC"},
    {"574", "GGGGGGGGGGTAAGAA"},
    {"575", "AAAAAAAAAATAAGAA"},
    {"576", "CCCCCCCCCCTAAGAA"},
    {"577", "TTTTTTTTTTTAAGAA"},
    {"632", "AAATGCAGTCGACGTACTATGCAGTC"},
    {"633", "CCCGTATCGACTATGCAGCGTATCGA"},
    {"634", "GGGACTCGATGCTCAGTAGACTCGAT"},
    {"635", "TTTCAGTACGATGACTGCTCAGTACG"},
    {"648", "TTTGCATTAAGAAATTAAAAAAGCTAAAAAAAAAA"},
    {"649", "AAAGCATTAAGAAATTAAAAAAGCTAAAAAAAAAA"},
    {"650", "GGGGCATTAAGAAATTAAAAAAGCTAAAAAAAAAA"},
    {"651", "CCCGCATTAAGAAATTAAAAAAGCTAAAAAAAAAA"},
    {"657", "GGGCATCTCGTATGCC"},
    {"662", "ACTGATCTCGTATGCC"},
    {"663", "GCTGATCTCGTATGCC"},
    {"664", "CAGCATCTCGTATGCC"},
    {"665", "TCTGATCTCGTATGCC"}
};

std::vector<std::string> ref = {"662", "663", "664", "665"};

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

std::vector<std::string> Parse(char *str, char delim)
{
    std::vector<std::string> tokens;
    char *ptr = str;
    while (ptr && *ptr != 0) {
        char *delimPtr = strchr(ptr, delim);
        if (delimPtr) {
            *delimPtr = 0;
            delimPtr++;
        }
        tokens.push_back(ptr);
        ptr = delimPtr;
    }

    return tokens;
}

void LoadSpotData(const char *filename, std::vector<SpotData> &spotData)
{
    std::vector<Signal> vals;
    int spotId;
    char name[256];
    int cycle;
    Signal s;
    int lastSpotId = 0;
    std::string spotName;

    FILE *fp = fopen(filename, "r");
    if (fp) {
        char line[1024];
        // read and ignore the header
        fgets(line, sizeof(line), fp);
        while (fgets(line, sizeof(line), fp)) {
            std::vector<std::string> tokens = Parse(line, ',');
            spotId = atoi(tokens[0].c_str());
            // if spot id is different, save the last spot data and clear things out
            if (spotId != lastSpotId) {
                if (lastSpotId != 0) {
                    SpotData d;
                    d.vals = vals;
                    d.spotName = spotName;
                    spotData.push_back(d);
                }
                lastSpotId = spotId;
                vals.clear();
            }
            // now read the new spot data and store
            spotName = tokens[1];
            for(int i=0;i<4;i++)
                s.v[i] = atof(tokens[3+i].c_str());
            vals.push_back(s);
        }
        fclose(fp);

        // save the last spot and vals
        SpotData d;
        d.vals = vals;
        d.spotName = spotName;
        spotData.push_back(d);
    }
}

double CallBases(char *dnaTemplate, std::vector<Signal> &measuredSignal, std::vector<double> &errorPerCycle, double ie, double cf, double dr)
{
    int numCycles = measuredSignal.size();

    Model4 m;
    m.Init(100);
    m.SetParams(ie, cf, dr);

    std::vector<Signal> dyeIntensities(numCycles);
    std::vector<double> totalSignal(numCycles, 0);

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

PhaseParams GridSearch(char *dnaTemplate, std::vector<Signal> &measuredSignal, std::vector<double> &errorPerCycle)
{
    // grid-search first, then call based on lowest error
    double drmin = 0.005;
    double drmax = 0.025;
    int    drnum = 5; // 10;
    double iemin = 0.13;
    double iemax = 0.19;
    int    ienum = 12; // 20;
    double cfmin = 0.07;
    double cfmax = 0.12;
    int    cfnum = 10; // 20;

    int numCycles = measuredSignal.size();
    PhaseParams params;
    double minerr = 99999.0;

    for(int dri=0;dri<drnum;dri++) {
        double drtest = drmin + (dri/(double)(drnum-1.0)) * (drmax-drmin);
        for(int cfi=0;cfi<cfnum;cfi++) {
            double cftest = cfmin + (cfi/(double)(cfnum-1.0)) * (cfmax-cfmin);
            for(int iei=0;iei<ienum;iei++) {
                double ietest = iemin + (iei/(double)(ienum-1.0)) * (iemax-iemin);
                double err = CallBases(dnaTemplate, measuredSignal, errorPerCycle, ietest, cftest, drtest);
                if (err < minerr) {
                    minerr = err;
                    params.ie = ietest;
                    params.cf = cftest;
                    params.dr = drtest;
                }
            }
        }
    }
    params.err = CallBases(dnaTemplate, measuredSignal, errorPerCycle, params.ie, params.cf, params.dr);
    return params;
}

int main(int argc, char *argv[])
{
    const char *spotFile = "color_transformed_spots.csv";
    bool gridsearch = true;
    double ie = 0.0, cf = 0.0, dr = 0.0;
    const char *fastQFileName = "out.fastq";

    int argcc = 1;
    while (argcc < argc) {
        switch (argv[argcc][1]) {
            case 'f':
                argcc++;
                spotFile = argv[argcc];
            break;

            case 'o':
                argcc++;
                fastQFileName = argv[argcc];
            break;

            case 'i':
                gridsearch = false;
                argcc++;
                ie = atof(argv[argcc]);
            break;

            case 'c':
                gridsearch = false;
                argcc++;
                cf = atof(argv[argcc]);
            break;

            case 'd':
                gridsearch = false;
                argcc++;
                dr = atof(argv[argcc]);
            break;
        }
        argcc++;
    }

    // load spots
    // std::vector<Signal> measuredSignal = LoadCycleIntensities("jmrdata1.csv");
    std::vector<SpotData> spotData;
    LoadSpotData(spotFile, spotData);
    int numSpots = spotData.size();
    printf("Loaded %d spots\n", numSpots);
    int numCycles = spotData[0].vals.size();

    // grid-search phase-correct each spot
    int numReadsAll = 0; // needs to be min 4 bases correct
    double qualScoreAll = 0.0;

    int numReadsHQ = 0;
    double qualScoreHQ = 0.0;
    double readLenHQ = 0;

    int num6Q7 = 0;

    FILE *fastQFile = fopen(fastQFileName, "w");

    char dnaTemplate[1024];
    char basecalls[1024];
    std::vector<double> errorPerCycle(numCycles, 0);
    std::vector<double> qualScoreList;
    for(int i=0;i<numSpots;i++) {
        PhaseParams params;
        if (gridsearch)
            params = GridSearch(dnaTemplate, spotData[i].vals, errorPerCycle);
        else {
            params.ie = ie;
            params.cf = cf;
            params.dr = dr;
            params.err = CallBases(dnaTemplate, spotData[i].vals, errorPerCycle, params.ie, params.cf, params.dr);
        }
        template2bases(dnaTemplate, basecalls);
        // std::string spotSequence = spotSeq[spotData[i].spotName];
        // map against reference genome
        int len = strlen(basecalls);
        int bestDist = -1;
        int bestRef = -1;
        for(int r=0;r<ref.size();r++) {
            std::string spotSequence = spotSeq[ref[r]];
            int dist = distance(spotSequence.c_str(), len, basecalls, len);
            if (bestDist == -1 || dist < bestDist) {
                bestDist = dist;
                bestRef = r;
            }
        }
        std::string spotSequence = spotSeq[ref[bestRef]];

        // assign a phred score to each base, using the error in measured vs predicted at each position
        std::string qscores;
        for(int c=0;c<numCycles;c++) {
            double err = errorPerCycle[c] * 2.0; // range is 0.0 to 1.0, 0.0 is very high quality Q30, 1.0 is Q3 (50% prob)
            if (err < 0.0) err = 0.0;
            if (err > 1.0) err = 1.0;
            double prob = 0.001 + err*(0.5-0.001);
            int QScore = int(-10.0 * log10(prob)); // truncate down to be fair
            char QChar = QScore + 33; // phred-33 encoding
            qscores += QChar;
        }

        int perfectLen = 0;
        while (perfectLen < len && spotSequence[perfectLen] == basecalls[perfectLen]) perfectLen++;

        double qualScore = -10.0 * log10((double)bestDist/(double)len);
        printf("spot: %d ie/cf/dr: %.3lf/%.3lf/%.3lf err: %.3lf template: %s qual: %s spot seq: %s dist: %d Q: %.3lf perfect: %d\n",
            i, params.ie, params.cf, params.dr, params.err, basecalls, qscores.c_str(), spotSequence.c_str(), bestDist, qualScore, perfectLen);

        if (fastQFile) {
            if (perfectLen > 4) {
                fprintf(fastQFile, "@spot_%06d\n%s\n+\n%s\n", i, basecalls, qscores.c_str());
            }
        }

        // see how long the read was until we make a 2nd error
        int numErrors = 0;
        int goodLen = 0;
        int goodEditDist = 0;
        for(int l=0;l<len && numErrors<2;l++) {
            bool isError = (basecalls[l] != spotSequence.c_str()[l]);
            if (isError)
                numErrors++;
            if (!isError && numErrors < 2) { // need to be careful here not to count an error towards the good readlength
                goodLen = l;
                goodEditDist = numErrors;
            }
        }

        if (goodLen >= 4) {
            if (goodEditDist == 0) { // I don't like this, but if our read ends with all errors, goodEditDist is 0 which puts our Q-score at +inf
                goodEditDist++;
                goodLen++;
            }
            double goodQualScore = -10.0 * log10((double)goodEditDist/(double)goodLen);
            numReadsHQ++;
            qualScoreHQ += goodQualScore;
            readLenHQ += goodLen;
        }

        if (goodLen >= 6) {
            num6Q7++;
        }

        // calc phred score and save some stats
        // this is simply counting a 8 cycles and tracking reads at lease 50% correct
        if (bestDist <= 4) {
            numReadsAll++;
            qualScoreAll += qualScore;
        }
    }

    if (numReadsAll > 0)
        qualScoreAll /= numReadsAll;
    if (numReadsHQ > 0) {
        qualScoreHQ /= numReadsHQ;
        readLenHQ /= numReadsHQ;
    }

    if (fastQFile)
        fclose(fastQFile);

    printf("%d %d Q %.3lf\n", numReadsAll, numCycles, qualScoreAll);
    printf("%d %.3lf Q %.3lf\n", numReadsHQ, readLenHQ, qualScoreHQ);
    printf("%d 6 Q 7\n", num6Q7);
}

