//
//  SCEPTTrC 1.2
//
//  Created by Jeffrey Hartgerink, summer 2021
//
//  Reads parameters from parameters.txt.
//  Reads a series of sequences in from seq_input.txt.
//  Scores a series of 1, 2 or 3 peptides based on length, propensity and pairwise amino acid interactions.
//  Produces Tm scores for all canonical registers. Highlights the best, second best and specificity of the system.
//  Writes a summary to output.txt.


#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <thread>
#include <time.h>

using namespace std;

struct parameterType
{
    double axial[27][27]; // 0 is undefined. 1 is "A", 26 is "Z".
    double lateral[27][27];
    double propensityX[27];
    double propensityY[27];
    double  A, B, C; // A + Bx + Cx^2 for the length parameter. Starting values of -82.57, 7.549, -0.0853.
    double charge;
    double Nterm;
    double Cterm;
    
    // These are true if we want to optimize them.
    bool    optLength;
    bool    optPropX[27];
    bool    optPropY[27];
    bool    optAxial[27][27];
    bool    optLat[27][27];
    
    // These are the "initial" values set primarily from experimentation.
    double exAxial[27][27];
    double exLateral[27][27];
    double exPropensityX[27];
    double exPropensityY[27];
    double exA, exB, exC; // A + Bx + Cx^2 for the length parameter. Starting values of -82.57, 7.549, -0.0853.
    double exCharge;
    double exNterm;
    double exCterm;
};


// This is a recursive function.
// It returns a double containing the maximum possible change in Tm from pairwise
// interactions from that start point of the call through to the end of the peptide strand.
// The final return should be the total maximum pairwise interactions from that strand.
double PairWiseCalc (double XPW[], double LPW[], short currentPair, short lastPair, short previousPWType, double currentPWSum, double BestPWTotal)
{
    bool NoPairT = false;
    bool AxPairT = false;
    bool LatPairT = false;
        
    if (currentPair > lastPair)
    {
        return BestPWTotal;
    }
    
    if ((not NoPairT) && (previousPWType != 0))
    {
        currentPWSum += 0;
        if (currentPair == lastPair)
        {
            if (currentPWSum > BestPWTotal)
            {
                BestPWTotal = currentPWSum;
            }
            return BestPWTotal;
        }
        NoPairT = true;
                
        // "We need to go deeper!"
        BestPWTotal = PairWiseCalc(XPW, LPW, (currentPair +1), lastPair, 0 /* the previous is no pair or zero */, currentPWSum, BestPWTotal);
    }
        
    if ((not LatPairT) && (previousPWType != 2))
    {
        // only sum stabilizing interactions. All possible destabilizing intereactions will be accounted for elsewhere.
        if (LPW[currentPair] > 0) currentPWSum += LPW[currentPair];
        if (currentPair == lastPair)
         {
             if (currentPWSum > BestPWTotal)
             {
                 BestPWTotal = currentPWSum;
             }
             return BestPWTotal;
         }
         LatPairT = true;
                  
     // "We need to go deeper!"
     BestPWTotal = PairWiseCalc(XPW, LPW, (currentPair +1), lastPair, 1 /* the previous is lateral or one */, currentPWSum, BestPWTotal);
    }
            
    if (not AxPairT)
    {
        // only sum stabilizing interactions. All possible destabilizing intereactions will be accounted for elsewhere.
        if (XPW[currentPair] > 0) currentPWSum += XPW[currentPair];
        if (currentPair == lastPair)
        {
            if (currentPWSum > BestPWTotal)
            {
                BestPWTotal = currentPWSum;
            }
            return BestPWTotal;
        }
        AxPairT = true;
        
        // "We need to go deeper!"
        BestPWTotal = PairWiseCalc(XPW, LPW, (currentPair +1), lastPair, 2 /* the previous is axial or two */, currentPWSum, BestPWTotal);
    }
    //cout << "terminal return" << endl;
    return BestPWTotal;
}


struct TripleHelix
{
    short   numPep;
    short   numAA;
    char    sequences[3][100];
    string  Nterm, Cterm;
    
    double  expTm;
    double  CCTm; // Tm of the "correct" composition. For A2B systems this requires at least one of each peptide (no homotrimers). For ABC systems this requires one of each peptide (no homotrimers and no A2B systems).
    
    double  HighTm;
    double  deviation;
    double  secTm;
    double  specificity;
    
    short   bestRegister[4];
    short   secRegister[4];
    short   CCRegister[4];
    
    short   netCharge[3][3][3][9];
    short   totalCharge[3][3][3][9];
    // These arrays contain the following:
    // [0] the value of the peptide in the leading position (peptide numbers are 0-2)
    // [1] the value of the peptide in the middle position
    // [2] the value of the peptide in the trailing position
    // [3] This is the value of the Offset with 9 possible values as follows:
    // Values range from 1-9. Nine Offsets or Staggers can be considered.
    // 0) {012} This is the canonical offset
    // --------------
    // 1) {015} -1 triplet.     trailing strand offset by an additional 3 amino acids
    // 2) {042} -1 triplet.     middle strand offset by an additional 3 amino acids
    // 3) {045} -1 triplet.     middle and trailing strand offset by an additional 3 amino acids
    // --------------
    // 4) {018} -2 triplets.    trailing strand offset by an additional 6 amino acids
    // 5) {048} -2 triplets.    middle stand offset by 3, trailing by 6
    // 6) {072} -2 triplets.    middle strand offset by 6
    // 7) {075} -2 triplets.    middle strand offset by 6, trailing by 3
    // 8) {078} -2 triplets.    middle and trailing strand offset by 6
    // Scoring will be handled by trimming the helix down to a canonical triple helix and then scoring based on this new smaller helix.
        
    double  bestPropensity;
    double  bestPairwise;
    double  Propensity[3][3][3][9];
    double  PairWise[3][3][3][9];
    double  Tm[3][3][3][9];
    
    short   XaaPos; // The position of the first Xaa amino acid.
    
    void initializeAll(void)
    {
        short a,b,c,d;
        
        bestPropensity = 0;
        bestPairwise = 0;
        HighTm = 0;
        deviation = 0;
        secTm = 0;
        specificity = 0;
        expTm = 0;
        CCTm = 0;
        numPep = 0;
        numAA = 0;
        Nterm = "initial";
        Cterm = "initial";
        
        for (a=0;a<3;a++)for(b=0;b<100;b++) sequences[a][b] = '.';
        for (a=0;a<3;a++) for(b=0;b<3;b++) for(c=0;c<3;c++) for(d=0;d<9;d++)
        {
            Propensity[a][b][c][d] = 0;
            PairWise[a][b][c][d] = 0;
            Tm[a][b][c][d] = 0;
            netCharge[a][b][c][d] = 0;
            totalCharge[a][b][c][d] = 0;
        }
        for (a=0;a<4;a++)
        {
            bestRegister[a] = 0;
            secRegister[a] = 0;
            CCRegister[a] = 0;
            
        }
    };
    
    bool isXaa(short position)
    {
        if (abs(position + (3 - this->XaaPos)) % 3 == 0) return true;
        return false;
    };
    
    bool isYaa(short position)
    {
        if (abs(position + (3 - this->XaaPos)) % 3 == 1) return true;
        return false;
    };
    
    bool isGly(short position)
    {
        if (abs(position + (3 - this->XaaPos)) % 3 == 2) return true;
        return false;
    };
    
    void dissect(void)
    {
        short x, y;
        cout << "numPep = " << numPep << endl;
        cout << "numAA =  " << numAA << endl;
        for (x=0;x<numPep;x++)
        {
            for (y=0;y<numAA;y++)
            {
                cout << sequences[x][y];
            }
            cout << endl;
        }
        cout << "termination: " << Nterm << " " << Cterm << endl;
        cout << "XaaPos = " << XaaPos << endl;
        cout << "expTm = " << expTm << ". CCTm = " << CCTm << ". Deviation = " << deviation << endl;
        cout << "CCregister = " << CCRegister[0] << "," << CCRegister[1] << "," << CCRegister[2] << "." << CCRegister[3] << endl;
        cout << "High Tm = " << HighTm << " = " << bestPropensity << " + " << bestPairwise << endl;
        cout << "Best register = " << bestRegister[0] << "," << bestRegister[1] << "," << bestRegister[2] << "." << bestRegister[3] <<endl;
        cout << "Second highest Tm = " << secTm << endl;
        cout << "Second Best register = " << secRegister[0] << "," << secRegister[1] << "," << secRegister[2] << "." << secRegister[3] << endl;
        cout << "Specificity = " << specificity << "." << endl;
        cout << endl;
    };
    
    void determine_reptition(void)
    {
        short x, xCount, yCount, zCount;
        bool goodPeptide;
        
        xCount = 0;
        yCount = 0;
        zCount = 0;
        goodPeptide = false;
        
        for (x=0;x<numAA;x++)
        {
            if ((x%3 == 0) && (sequences[0][x] == 'G')) xCount++;
            if ((x%3 == 1) && (sequences[0][x] == 'G')) yCount++;
            if ((x%3 == 2) && (sequences[0][x] == 'G')) zCount++;
        }
        if (xCount >= (numAA)/3)
        {
            // we are ok
            // cout << n << " Gly was found in the first position and in every subsequent i+3 position. xCount = " << xCount << endl;
            goodPeptide = true;
            XaaPos = 1;
        }
        if (yCount >= (numAA)/3)
        {
            // we are ok
            // cout << n << " Gly was found in the second position and in every subsequent i+3 position. yCount = " << yCount << endl;
            goodPeptide = true;
            XaaPos = 2;
        }
        if (zCount >= (numAA)/3)
        {
            // we are ok
            // cout << n << " Gly was found in the third position and in every subsequent i+3 position. zCount = " << zCount << endl;
            goodPeptide = true;
            XaaPos = 0;
        }
        if (not goodPeptide)
        {
            cout << "This peptide does not appear to have a Gly every third residue!" << endl;
            dissect();
        }
    };
        
    void userOutput(void)
    {
        short a,b,c;
        short x,y;
        cout << "----------------------------------------------------------" << endl;
        cout << "Number of peptides: " << numPep << endl;
        cout << "Number of amino acids: " << numAA << endl;
        cout << Nterm << "...peptide..." << Cterm << endl;
        cout << "Experimental Tm = " << expTm << endl;
        cout << "Deviation (Tm(predicted) - Tm(experimental)) = " <<  CCTm - expTm << endl;
        cout << endl;
        cout << "The most stable register/composition is {" << bestRegister[0] <<  bestRegister[1]  << bestRegister[2] << "}. Tm = " << HighTm << "." << endl;
        cout << "Total charge on {" << bestRegister[0] <<  bestRegister[1]  << bestRegister[2] << "} = " << totalCharge[bestRegister[0]][bestRegister[1]][bestRegister[2]][bestRegister[3]] << endl;
        cout << "Net charge on {" << bestRegister[0] <<  bestRegister[1]  << bestRegister[2] << "} = " << netCharge[bestRegister[0]][bestRegister[1]][bestRegister[2]][bestRegister[3]] << endl;
        if (numPep == 2)
        {
            if ((bestRegister[0] != bestRegister[1]) || (bestRegister[0] != bestRegister[2]) || (bestRegister[1] != bestRegister[2]))
            {
                //cout << "This matches the input diversity of peptides (2)." << endl;
            }
            else
            {
                cout << "\x1b[1m\x1b[31mWARNING: The most stable register/composition does not include all the peptides you input.\x1b[0m" << endl;
            }
        }
        if (numPep == 3)
        {
            if ((bestRegister[0] != bestRegister[1]) && (bestRegister[0] != bestRegister[2]) && (bestRegister[1] != bestRegister[2]))
            {
                //cout << "This matches the input diversity of peptides (3)." << endl;
            }
            else
            {
                cout << "\x1b[1m\x1b[31mWARNING: The most stable register/composition does not include all the peptides you input.\x1b[0m" << endl;
            }
        }
        
        
        // 0) {012} This is the canonical offset
        // --------------
        // 1) {015} -1 triplet.     trailing strand offset by an additional 3 amino acids
        // 2) {042} -1 triplet.     middle strand offset by an additional 3 amino acids
        // 3) {045} -1 triplet.     middle and trailing strand offset by an additional 3 amino acids
        // --------------
        // 4) {018} -2 triplets.    trailing strand offset by an additional 6 amino acids
        // 5) {048} -2 triplets.    middle stand offset by 3, trailing by 6
        // 6) {072} -2 triplets.    middle strand offset by 6
        // 7) {075} -2 triplets.    middle strand offset by 6, trailing by 3
        // 8) {078} -2 triplets.    middle and trailing strand offset by 6
        
        
        
        x = bestRegister[0];
        cout << bestRegister[0] << ": ";
        for (y=0;y<numAA;y++)
        {
            if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
            if (sequences[x][y] == 'R') cout << "\x1b[34m";
            if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
            if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
            if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
            cout << sequences[x][y];
            cout << "\x1b[0m";
        }
        cout << endl;
        
        x = bestRegister[1];
        if ((bestRegister[3] == 0) || bestRegister[3] == 1) cout << bestRegister[1] << ":  ";
        if ((bestRegister[3] == 2) || bestRegister[3] == 3) cout << bestRegister[1] << ":     ";
        for (y=0;y<numAA;y++)
        {
            if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
            if (sequences[x][y] == 'R') cout << "\x1b[34m";
            if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
            if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
            if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
            cout << sequences[x][y];
            cout << "\x1b[0m";
        }
        cout << endl;
        
        x = bestRegister[2];
        if ((bestRegister[3] == 0) || bestRegister[3] == 2) cout << bestRegister[2] << ":   ";
        if ((bestRegister[3] == 1) || bestRegister[3] == 3) cout << bestRegister[2] << ":      ";
        for (y=0;y<numAA;y++)
        {
            if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
            if (sequences[x][y] == 'R') cout << "\x1b[34m";
            if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
            if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
            if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
            cout << sequences[x][y];
            cout << "\x1b[0m";
        }
        cout << endl;
        
        if (numPep != 1)
        {
            cout << endl;
            cout << "The second most stable register/composition is {" << secRegister[0] << secRegister[1] << secRegister[2] << "}. Tm = " << secTm << "." << endl;
            x = secRegister[0];
            cout << secRegister[0] << ": ";
            for (y=0;y<numAA;y++)
            {
                if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
                if (sequences[x][y] == 'R') cout << "\x1b[34m";
                if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
                if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
                if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
                cout << sequences[x][y];
                cout << "\x1b[0m";
            }
            cout << endl;
            x = secRegister[1];
            cout << secRegister[1] << ":  ";
            for (y=0;y<numAA;y++)
            {
                if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
                if (sequences[x][y] == 'R') cout << "\x1b[34m";
                if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
                if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
                if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
                cout << sequences[x][y];
                cout << "\x1b[0m";
            }
            cout << endl;
            x = secRegister[2];
            cout << secRegister[2] << ":   ";
            for (y=0;y<numAA;y++)
            {
                if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
                if (sequences[x][y] == 'R') cout << "\x1b[34m";
                if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
                if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
                if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
                cout << sequences[x][y];
                cout << "\x1b[0m";
            }
            cout << endl;
            cout << endl;
            cout << "The specificity is = " << specificity << "." << endl;
            cout << endl;
        }
        
        cout << "Melting temperatures of all canonical registers." << endl;
        cout << "Best in blue, second best in red." << endl;
        cout << "Tm < 10C faded to indicate experimentally unreliable (frequently will not fold)." << endl;
        for (a=0;a<numPep;a++)
        {
            for (b=0;b<numPep;b++)
            {
                for (c=0;c<numPep;c++)
                {
                    cout << "\x1b[0m";
                    if ((a == bestRegister[0]) && (b == bestRegister[1]) && (c == bestRegister[2])) cout << "\x1b[1m\x1b[34m";
                    if ((a == secRegister[0]) && (b == secRegister[1]) && (c == secRegister[2])) cout << "\x1b[1m\x1b[31m";
                    if (Tm[a][b][c][0] < 10) cout << "\x1b[2m";
                    cout << "{" << a << b << c << "} = " << Tm[a][b][c][0] << endl;
                    cout << "\x1b[0m";
                }
            }
            cout << endl;
        }
    };
    
};
    
parameterType ReadParameters(void)
{
    // Read parameters from file.
    parameterType parameters;
    short x, y;
    char SingleAminoAcid, parameterChar;
    short optValue;
    string StringLine;
    
    ifstream parameterFile("parameters.txt");
    if (!parameterFile.is_open()) cout << "We couldn't open the parameters.txt file." << endl;
    
    // Zero initial values.
    for (x=0;x<27;x++)
    {
        parameters.propensityX[x] = 0;
        parameters.propensityY[x] = 0;
        
        parameters.exPropensityX[x] = 0;
        parameters.exPropensityY[x] = 0;
        
        parameters.optPropX[x] = false;
        parameters.optPropY[x] = false;
        
        for (y=0;y<27;y++)
        {
            parameters.axial[x][y] = 0;
            parameters.lateral[x][y] = 0;
            
            parameters.exAxial[x][y] = 0;
            parameters.exLateral[x][y] = 0;
            
            parameters.optAxial[x][y] = false;
            parameters.optLat[x][y] = false;
        }
    }
    
    getline(parameterFile, StringLine);
    cout << "Parameter File: " << StringLine << endl;
  
    while (!parameterFile.eof())
    {
        getline(parameterFile, StringLine);
        if (StringLine == "Length")
        {
            //cout << "we found Length" << endl;
            parameterFile >> parameters.A;
            parameterFile >> parameters.B;
            parameterFile >> parameters.C;
            
        }
        if (StringLine == "XaaPropensity")
        {
            //cout << "we found XaaPropensity" << endl;
            for (x=1;x<27;x++)
            {
                parameterFile >> SingleAminoAcid;
                parameterFile >> parameters.propensityX[(short)SingleAminoAcid-64];
            }
        }
        if (StringLine == "YaaPropensity")
        {
            //cout << "we found YaaPropensity" << endl;
            for (x=1;x<27;x++)
            {
                parameterFile >> SingleAminoAcid;
                parameterFile >> parameters.propensityY[(short)SingleAminoAcid-64];
            }
        }
        if (StringLine == "PairwiseLateral")
        {
            //cout << "we found PairwiseLateral" << endl;
            getline(parameterFile, StringLine);
            for (y=1; y<27; y++)
            {
                parameterFile >> parameterChar;
                for (x=1; x<27; x++)
                {
                    parameterFile >> parameters.lateral[y][x];
                }
            }
        }
        
        if (StringLine == "PairwiseAxial")
        {
            //cout << "we found PairwiseAxial" << endl;
            getline(parameterFile, StringLine);
            for (y=1; y<27; y++)
            {
                parameterFile >> parameterChar;
                for (x=1; x<27; x++)
                {
                    parameterFile >> parameters.axial[y][x];
                }
            }
        }
        if (StringLine == "EOF")
        {
            //cout << "we found EOF" << endl;
        }
    }
    parameterFile.close();
    
    parameterFile.open("parameters_exp.txt");
    if (!parameterFile.is_open()) cout << "We couldn't open the parameters.txt file." << endl;
      
    getline(parameterFile, StringLine);
    cout << "Experimental Parameter File: " << StringLine << endl;
    
    while (!parameterFile.eof())
      {
          getline(parameterFile, StringLine);
          if (StringLine == "Length")
          {
              //cout << "we found Length" << endl;
              parameterFile >> parameters.exA;
              parameterFile >> parameters.exB;
              parameterFile >> parameters.exC;
              
          }
          if (StringLine == "XaaPropensity")
          {
              //cout << "we found XaaPropensity" << endl;
              for (x=1;x<27;x++)
              {
                  parameterFile >> SingleAminoAcid;
                  parameterFile >> parameters.exPropensityX[(short)SingleAminoAcid-64];
              }
          }
          if (StringLine == "YaaPropensity")
          {
              //cout << "we found YaaPropensity" << endl;
              for (x=1;x<27;x++)
              {
                  parameterFile >> SingleAminoAcid;
                  parameterFile >> parameters.exPropensityY[(short)SingleAminoAcid-64];
              }
          }
          if (StringLine == "PairwiseLateral")
          {
              //cout << "we found PairwiseLateral" << endl;
              getline(parameterFile, StringLine);
              for (y=1; y<27; y++)
              {
                  parameterFile >> parameterChar;
                  for (x=1; x<27; x++)
                  {
                      parameterFile >> parameters.exLateral[y][x];
                  }
              }
          }
          
          if (StringLine == "PairwiseAxial")
          {
              //cout << "we found PairwiseAxial" << endl;
              getline(parameterFile, StringLine);
              for (y=1; y<27; y++)
              {
                  parameterFile >> parameterChar;
                  for (x=1; x<27; x++)
                  {
                      parameterFile >> parameters.exAxial[y][x];
                  }
              }
          }
          if (StringLine == "EOF")
          {
              //cout << "we found EOF" << endl;
          }
      }
      parameterFile.close();
    
    // Optionally read in which values will be optimized.
    
    parameterFile.open("opt_list.txt");
    
    getline(parameterFile, StringLine);
    cout << "Optimization List: " << StringLine << endl;
    
    if (!parameterFile.is_open()) cout << "We couldn't open the parameters.txt file." << endl;
    
    while (!parameterFile.eof())
    {
        getline(parameterFile, StringLine);
        if (StringLine == "Length")
        {
            //cout << "we found Length" << endl;
            parameterFile >> optValue;
            parameterFile >> optValue;
            parameterFile >> optValue;
            if (optValue == 1) parameters.optLength = true; else parameters.optLength = false;
            
        }
        if (StringLine == "XaaPropensity")
        {
            //cout << "we found XaaPropensity" << endl;
            for (x=1;x<27;x++)
            {
                parameterFile >> SingleAminoAcid;
                parameterFile >> optValue;
                if (optValue == 1) parameters.optPropX[x] = true; else parameters.optPropX[x] = false;
            }
        }
        if (StringLine == "YaaPropensity")
        {
            //cout << "we found YaaPropensity" << endl;
            for (x=1;x<27;x++)
            {
                parameterFile >> SingleAminoAcid;
                parameterFile >> optValue;
                if (optValue == 1) parameters.optPropY[x] = true; else parameters.optPropY[x] = false;
            }
        }
        if (StringLine == "PairwiseLateral")
        {
            //cout << "we found PairwiseLateral" << endl;
            getline(parameterFile, StringLine);
            for (y=1; y<27; y++)
            {
                parameterFile >> parameterChar;
                for (x=1; x<27; x++)
                {
                    parameterFile >> optValue;
                    if (optValue == 1) parameters.optLat[y][x] = true; else parameters.optLat[y][x] = false;
                }
            }
        }
        
        if (StringLine == "PairwiseAxial")
        {
            //cout << "we found PairwiseAxial" << endl;
            getline(parameterFile, StringLine);
            for (y=1; y<27; y++)
            {
                parameterFile >> parameterChar;
                for (x=1; x<27; x++)
                {
                    parameterFile >> optValue;
                    if (optValue == 1) parameters.optAxial[y][x] = true; else parameters.optAxial[y][x] = false;
                }
            }
        }
        if (StringLine == "EOF")
        {
            //cout << "we found EOF" << endl;
        }
    }
    parameterFile.close();
    
    
    return parameters;
}

void DisplayParameters(parameterType parameters)
{
    short x, y;
  
    cout << endl;
    cout << "OptXaaPropensity" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.optPropX[x] << endl;
    cout << "OptYaaPropensity" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.optPropY[x] << endl;
    cout << "OptPairwiseLateral" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.optLat[x][y] << "\t";
        }
        cout << endl;
    }
    
    cout << "OptPairwiseAxial" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.optAxial[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "EOF" << endl;
    cout << endl;
    cout << "----------------------------------------" << endl;
    cout << endl;
    cout << "LengthEx" << endl;
    cout << parameters.exA << endl;
    cout << parameters.exB << endl;
    cout << parameters.exC << endl;
    cout << "XaaPropensityEx" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.exPropensityX[x] << endl;
    cout << "YaaPropensityEx" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.exPropensityY[x] << endl;
    cout << "PairwiseLateralEX" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.exLateral[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "PairwiseAxial" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.exAxial[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "EOF" << endl;
    cout << endl;
    cout << "----------------------------------------" << endl;
    cout << endl;
    cout << "Length" << endl;
    cout << parameters.A << endl;
    cout << parameters.B << endl;
    cout << parameters.C << endl;
    cout << "XaaPropensity" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.propensityX[x] << endl;
    cout << "YaaPropensity" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.propensityY[x] << endl;
    cout << "PairwiseLateral" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.lateral[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "PairwiseAxial" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.axial[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "EOF" << endl;
    cout << endl;
}

void WriteParameters(parameterType parameters)
{
    // Write parameters to file.
    ofstream newParameters("newParameters.txt");
    short x, y;
    newParameters << "Length" << endl;
    newParameters << parameters.A << endl;
    newParameters << parameters.B << endl;
    newParameters << parameters.C << endl;
    newParameters << "XaaPropensity" << endl;
    for (x=1;x<27;x++) newParameters << char(64+x) << "\t" << parameters.propensityX[x] << endl;
    newParameters << "YaaPropensity" << endl;
    for (x=1;x<27;x++) newParameters << char(64+x) << "\t" << parameters.propensityY[x] << endl;
    newParameters << "PairwiseLateral" << endl;
    for (x=1;x<27;x++) newParameters << "\t" << char(64+x);
    newParameters << endl;
    for (y=1;y<27;y++)
    {
        newParameters << char(64+y) << "\t";
        for (x=1;x<27;x++)
        {
            newParameters << parameters.lateral[y][x] << "\t";
        }
        newParameters << endl;
    }
    
    newParameters << "PairwiseAxial" << endl;
    for (x=1;x<27;x++) newParameters << "\t" << char(64+x);
    newParameters << endl;
    for (y=1;y<27;y++)
    {
        newParameters << char(64+y) << "\t";
        for (x=1;x<27;x++)
        {
            newParameters << parameters.axial[y][x] << "\t";
        }
        newParameters << endl;
    }
    newParameters << "EOF" << endl;
    
    newParameters.close();
    
    
    
}


// Will be better to pass a pointer to parameters rather than the entire struct, but this works for now.

void ScoreHelix (parameterType parameters, TripleHelix * theHelix)
{
    short a, b, c, d, x, y, i;
        
    for (a=0;a<3;a++)for(b=0;b<3;b++)for(c=0;c<3;c++)for(d=0;d<9;d++)
    {
        theHelix->Propensity[a][b][c][d] = 0;
        theHelix->PairWise[a][b][c][d] = 0;
        theHelix->Tm[a][b][c][d] = 0;
    }
    
    double XinteractionThread[20];
    double LinteractionThread[20];
    for (x=0; x<20; x++)
    {
        XinteractionThread[x] = 0;
        LinteractionThread[x] = 0;
    }
    
    // this triple helix will hold the temporary values for all trimmed non-canonical offsets.
    TripleHelix trimmedHelix;
    
    short numYaa = 0;
    numYaa = theHelix->numAA / 3;
    
    double maxTm, secondBestTm;
    double bestReg[4], secondBestReg[4];
    double lengthBasis = 0;
    secondBestTm = -2000;
    secondBestReg[0] = 5;
    secondBestReg[1] = 5;
    secondBestReg[2] = 5;
    secondBestReg[3] = 10;
    
    maxTm = -1000;
    bestReg[0] = 6;
    bestReg[1] = 6;
    bestReg[2] = 6;
    bestReg[3] = 11;
    
    theHelix->CCTm = -1500;
    theHelix->CCRegister[0] = 7;
    theHelix->CCRegister[1] = 7;
    theHelix->CCRegister[2] = 7;
    theHelix->CCRegister[3] = 12;
    
    // NOTE: d loop is being short circuited to only look at canonical registers here! //
    for (a=0; a<theHelix->numPep; a++) for (b=0; b<theHelix->numPep; b++) for (c=0; c<theHelix->numPep; c++) for (d=0;d<1;d++)
    {
        // // // // // // // //
        // Offset Switch       //
        // // // // // // // //
        trimmedHelix.numPep = theHelix->numPep;
        trimmedHelix.Nterm = theHelix->Nterm;
        trimmedHelix.Cterm = theHelix->Cterm;
        trimmedHelix.XaaPos = theHelix->XaaPos;
        // Trim triple helix for offset scoring.
        switch(d)
        {
            case 0:
                // 0 is {012}, canonical and therefore needs no trimming.
                // Still need to move thisHelix into trimmedHelix for scoring below
                for(y=0;y<theHelix->numAA;y++) trimmedHelix.sequences[a][y] = theHelix->sequences[a][y];
                for(y=0;y<theHelix->numAA;y++) trimmedHelix.sequences[b][y] = theHelix->sequences[b][y];
                for(y=0;y<theHelix->numAA;y++) trimmedHelix.sequences[c][y] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA;
                break;
                
            case 1:
                // 1 is {015}
                // trailing strand is offset 3 and therefore needs to be cut at C term while leading and middle strands need to lose first triplet
                for(y=3;y<theHelix->numAA;y++) trimmedHelix.sequences[a][y-3] = theHelix->sequences[a][y];
                for(y=3;y<theHelix->numAA;y++) trimmedHelix.sequences[b][y-3] = theHelix->sequences[b][y];
                for(y=0;y<theHelix->numAA-3;y++) trimmedHelix.sequences[c][y] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA - 3;
                break;

            case 2:
                // 3 is {042}
                // middle strand is offset 3, others are not shifted. Cut lead and trail at N-term, mid at C term
                for(y=3;y<theHelix->numAA-0;y++) trimmedHelix.sequences[a][y-3] = theHelix->sequences[a][y];
                for(y=0;y<theHelix->numAA-3;y++) trimmedHelix.sequences[b][y-0] = theHelix->sequences[b][y];
                for(y=3;y<theHelix->numAA-0;y++) trimmedHelix.sequences[c][y-3] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA - 3;
                break;
                
            case 3:
                // 4 is {045}
                // mid and trail are offset 3. Cut lead at N-term, mid and trail at C term
                for(y=3;y<theHelix->numAA-0;y++) trimmedHelix.sequences[a][y-3] = theHelix->sequences[a][y];
                for(y=0;y<theHelix->numAA-3;y++) trimmedHelix.sequences[b][y-0] = theHelix->sequences[b][y];
                for(y=0;y<theHelix->numAA-3;y++) trimmedHelix.sequences[c][y-0] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA - 3;
                break;
                
            case 4:
                // 2 is {018}
                for(y=6;y<theHelix->numAA;y++) trimmedHelix.sequences[a][y-6] = theHelix->sequences[a][y];
                for(y=6;y<theHelix->numAA;y++) trimmedHelix.sequences[b][y-6] = theHelix->sequences[b][y];
                for(y=0;y<theHelix->numAA-6;y++) trimmedHelix.sequences[c][y] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA - 6;
                break;
                
            case 5:
                // 5 is {048}
                // mid by 3, trail by 6.
                for(y=6;y<theHelix->numAA-0;y++) trimmedHelix.sequences[a][y-6] = theHelix->sequences[a][y];
                for(y=3;y<theHelix->numAA-3;y++) trimmedHelix.sequences[b][y-3] = theHelix->sequences[b][y];
                for(y=0;y<theHelix->numAA-6;y++) trimmedHelix.sequences[c][y-0] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA - 6;
                break;
            case 6:
                // 6 is {072}
                // mid by 6, trail by 3
                for(y=6;y<theHelix->numAA-0;y++) trimmedHelix.sequences[a][y-6] = theHelix->sequences[a][y];
                for(y=0;y<theHelix->numAA-6;y++) trimmedHelix.sequences[b][y-0] = theHelix->sequences[b][y];
                for(y=3;y<theHelix->numAA-6;y++) trimmedHelix.sequences[c][y-3] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA - 6;
                break;
            case 7:
                // 7 is {075}
                // mid by 6, trail by 3
                for(y=6;y<theHelix->numAA-0;y++) trimmedHelix.sequences[a][y-6] = theHelix->sequences[a][y];
                for(y=0;y<theHelix->numAA-6;y++) trimmedHelix.sequences[b][y-0] = theHelix->sequences[b][y];
                for(y=3;y<theHelix->numAA-3;y++) trimmedHelix.sequences[c][y-3] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA - 6;
                break;
            case 8:
                // 8 is {078}
                // mid by 6, trail by 6
                for(y=6;y<theHelix->numAA-0;y++) trimmedHelix.sequences[a][y-6] = theHelix->sequences[a][y];
                for(y=0;y<theHelix->numAA-6;y++) trimmedHelix.sequences[b][y-0] = theHelix->sequences[b][y];
                for(y=0;y<theHelix->numAA-6;y++) trimmedHelix.sequences[c][y-0] = theHelix->sequences[c][y];
                trimmedHelix.numAA = theHelix->numAA - 6;
                break;
            default:
                cout << "Value of d out of bounds. d = " << d << endl;
                break;
                
        }
        
        // // // // //
        // Length   //
        // // // // //
        if (trimmedHelix.numAA >50)
        {
            lengthBasis = parameters.A + (parameters.B*50) + (parameters.C*50*50);
        }
        else
        {
            if (d == 0) lengthBasis = parameters.A + (parameters.B*trimmedHelix.numAA) + (parameters.C*trimmedHelix.numAA*trimmedHelix.numAA);
            if ((d == 1) || (d == 2) || (d == 3)) lengthBasis = parameters.A + (parameters.B*(trimmedHelix.numAA+1)) + (parameters.C*(trimmedHelix.numAA+1)*(trimmedHelix.numAA+1));
        }
        //cout << "lengthBasis = " << lengthBasis << endl;
        theHelix->Propensity[a][b][c][d] = lengthBasis;
        
        // // // // // //
        // TERMINATION //
        // // // // // //
        if (trimmedHelix.Nterm == "n") theHelix->Propensity[a][b][c][d] -= 1.8;
        if (trimmedHelix.Cterm == "c") theHelix->Propensity[a][b][c][d] -= 1.8;
        
        // // // // // // // // //
        // Tyrosine/Tryptophan Termination //
        // // // // // // // // //
        if ((trimmedHelix.sequences[a][0] == 'Y') && (trimmedHelix.sequences[b][0] == 'Y') && (trimmedHelix.sequences[c][0] == 'Y'))
        {
            theHelix->Propensity[a][b][c][d] += 3;
        }
        if ((trimmedHelix.sequences[a][trimmedHelix.numAA-1] == 'Y') && (trimmedHelix.sequences[b][trimmedHelix.numAA-1] == 'Y') && (trimmedHelix.sequences[c][trimmedHelix.numAA-1] == 'Y'))
        {
            theHelix->Propensity[a][b][c][d] += 3;
        }
        if ((trimmedHelix.sequences[a][0] == 'W') && (trimmedHelix.sequences[b][0] == 'W') && (trimmedHelix.sequences[c][0] == 'W'))
        {
            theHelix->Propensity[a][b][c][d] += 3;
        }
        if ((trimmedHelix.sequences[a][trimmedHelix.numAA-1] == 'W') && (trimmedHelix.sequences[b][trimmedHelix.numAA-1] == 'W') && (trimmedHelix.sequences[c][trimmedHelix.numAA-1] == 'W'))
        {
            theHelix->Propensity[a][b][c][d] += 3;
        }
        
        //cout << "capping mod = " << Propensity[a][b][c] << endl;
        
        // // // // // // // // // // //
        // Terminal Hydrogen Bonding  //
        // // // // // // // // // // //
        // if (not theHelix->isXYG) Propensity[a][b][c][d] -= 3.6;
        if (not trimmedHelix.isXaa(0)) theHelix->Propensity[a][b][c][d] -= 1.8;
        if (not trimmedHelix.isGly(theHelix->numAA-1)) theHelix->Propensity[a][b][c][d] -= 1.8;
        
        //cout << "terminal H-bond mod = " << Propensity[a][b][c] << endl;
        
        // // // // // // // //
        // Single AA Score   //
        // // // // // // // //
            
        for (x=0;x<theHelix->numAA;x++)
        {
            if ((trimmedHelix.sequences[a][x] == 'K') || (trimmedHelix.sequences[a][x] == 'R'))
            {
                theHelix->netCharge[a][b][c][d]++;
                theHelix->totalCharge[a][b][c][d]++;
            }
            if ((trimmedHelix.sequences[a][x] == 'E') || (trimmedHelix.sequences[a][x] == 'D'))
            {
                theHelix->netCharge[a][b][c][d]--;
                theHelix->totalCharge[a][b][c][d]++;
            }
            if ((trimmedHelix.sequences[b][x] == 'K') || (trimmedHelix.sequences[b][x] == 'R'))
            {
                theHelix->netCharge[a][b][c][d]++;
                theHelix->totalCharge[a][b][c][d]++;
            }
            if ((trimmedHelix.sequences[b][x] == 'E') || (trimmedHelix.sequences[b][x] == 'D'))
            {
                theHelix->netCharge[a][b][c][d]--;
                theHelix->totalCharge[a][b][c][d]++;
            }
            if ((trimmedHelix.sequences[c][x] == 'K') || (trimmedHelix.sequences[c][x] == 'R'))
            {
                theHelix->netCharge[a][b][c][d]++;
                theHelix->totalCharge[a][b][c][d]++;
            }
            if ((trimmedHelix.sequences[c][x] == 'E') || (trimmedHelix.sequences[c][x] == 'D'))
            {
                theHelix->netCharge[a][b][c][d]--;
                theHelix->totalCharge[a][b][c][d]++;
            }
            
            if ((x>2) && (x<(theHelix->numAA-2))) // not the tips
            {
            // cout << trimmedHelix.sequences[a][x] << " " << (short)trimmedHelix.sequences[a][x]-64 << endl;
            if (trimmedHelix.isXaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[a][x]-64];
            if (trimmedHelix.isYaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[a][x]-64];
            
            if (trimmedHelix.isXaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[b][x]-64];
            if (trimmedHelix.isYaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[b][x]-64];
            
            if (trimmedHelix.isXaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[c][x]-64];
            if (trimmedHelix.isYaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[c][x]-64];
            }
            else // the tips
            {
                if (trimmedHelix.isXaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[a][x]-64]/3;
                if (trimmedHelix.isYaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[a][x]-64]/3;
                
                if (trimmedHelix.isXaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[b][x]-64]/3;
                if (trimmedHelix.isYaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[b][x]-64]/3;
                
                if (trimmedHelix.isXaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[c][x]-64]/3;
                if (trimmedHelix.isYaa(x)) theHelix->Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[c][x]-64]/3;
            }
        }
        
        // Charge Scoring
        if (abs(theHelix->netCharge[a][b][c][d]) > 6)
        {
            theHelix->Propensity[a][b][c][d] -= ((abs(theHelix->netCharge[a][b][c][d]) - 6)/3);
        }

        // // // // // // // //
        // set pairwise Tm   //
        // // // // // // // //
        
        // FIRST THREAD
        i = 0;
        for (x=0;x<theHelix->numAA;x++)
        {
            if (theHelix->isYaa(x))
            {
                // Populate all values for the First Interaction Thread (both axial and lateral options)
                // a, b & c are the peptide number of this particular composition / registration
                // i is tracking the number of Yaa's
                // x is tracking the amino acid position of the peptide
                if ((x+2) < trimmedHelix.numAA) XinteractionThread[i] = parameters.axial[(short)trimmedHelix.sequences[a][x]-64][(short)trimmedHelix.sequences[b][x+2]-64];
                    else XinteractionThread[i] = 0;
                if ((x-1) >= 0) LinteractionThread[i] = parameters.lateral[(short)trimmedHelix.sequences[a][x]-64][(short)trimmedHelix.sequences[b][x-1]-64];
                    else LinteractionThread[i] = 0;
                i++;
            }
        }
                    
        // find best combination of stabilizing interactions
        theHelix->PairWise[a][b][c][d] = PairWiseCalc (XinteractionThread, LinteractionThread, 0, numYaa, 9, 0, 0);
        
        // force *ALL* possible destabilizing interactions
        for (x=0;x<numYaa;x++)
        {
            //cout << "Axial Interaction Thread = " << XinteractionThread[i] << " ";
            if (XinteractionThread[x] < 0) theHelix->PairWise[a][b][c][d] += XinteractionThread[x];
            if (LinteractionThread[x] < 0) theHelix->PairWise[a][b][c][d] += LinteractionThread[x];
        }
        
        // SECOND THREAD
        i = 0;
        for (x=0;x<trimmedHelix.numAA;x++)
        {
            if (trimmedHelix.isYaa(x))
            {
                // Populate all values for the Second Interaction Thread (both axial and lateral options)
                if ((x+2) < theHelix->numAA)
                    XinteractionThread[i] = parameters.axial[(short)trimmedHelix.sequences[b][x]-64][(short)trimmedHelix.sequences[c][x+2]-64];
                else
                    XinteractionThread[i] = 0;
                if ((x-1) >= 0)
                    LinteractionThread[i] = parameters.lateral[(short)trimmedHelix.sequences[b][x]-64][(short)trimmedHelix.sequences[c][x-1]-64];
                else
                    LinteractionThread[i] = 0;
                i++;
            }
        }
        
        // find best combination of stabilizing interactions
        theHelix->PairWise[a][b][c][d] += PairWiseCalc (XinteractionThread, LinteractionThread, 0, numYaa, 9, 0, 0);
        
        // force *ALL* possible destabilizing interactions
        for (x=0;x<numYaa;x++)
        {
            if (XinteractionThread[x] < 0) theHelix->PairWise[a][b][c][d] += XinteractionThread[x];
            if (LinteractionThread[x] < 0) theHelix->PairWise[a][b][c][d] += LinteractionThread[x];
        }
        
        // THIRD THREAD
        i = 0;
        for (x=0;x<trimmedHelix.numAA;x++)
        {
            if (trimmedHelix.isYaa(x))
            {
                // Populate all values for the Third Interaction Thread (both axial and lateral options)
                if ((x+5) < theHelix->numAA)
                    XinteractionThread[i] = parameters.axial[(short)trimmedHelix.sequences[c][x]-64][(short)trimmedHelix.sequences[a][x+5]-64];
                else
                    XinteractionThread[i] = 0;
                if ((x+2) < theHelix->numAA)
                    LinteractionThread[i] = parameters.lateral[(short)trimmedHelix.sequences[c][x]-64][(short)trimmedHelix.sequences[a][x+2]-64];
                else
                    LinteractionThread[i] = 0;
                i++;
            }
        }
        
        // find best combination of stabilizing interactions
        theHelix->PairWise[a][b][c][d] += PairWiseCalc (XinteractionThread, LinteractionThread, 0, numYaa, 9, 0, 0);
        // force *ALL* possible destabilizing interactions
        for (x=0;x<numYaa;x++)
        {
            if (XinteractionThread[x] < 0) theHelix->PairWise[a][b][c][d] += XinteractionThread[x];
            if (LinteractionThread[x] < 0) theHelix->PairWise[a][b][c][d] += LinteractionThread[x];
        }
        
        //cout << "Pairwise Mod = " << PairWise[a][b][c] << endl;
        
        // Calculate the Tm for this composition / register determined.
        theHelix->Tm[a][b][c][d] = theHelix->Propensity[a][b][c][d] + theHelix->PairWise[a][b][c][d];
        
        //cout << "Tm = " << Tm[a][b][c] << " = " << Propensity[a][b][c] << " + " << PairWise[a][b][c] << endl;
        
        // Determine if this is the best and remember the best & second best registers.
        if (theHelix->Tm[a][b][c][d] >= maxTm)
        {
            secondBestTm = maxTm;
            secondBestReg[0] = bestReg[0];
            secondBestReg[1] = bestReg[1];
            secondBestReg[2] = bestReg[2];
            secondBestReg[3] = bestReg[3];
            
            maxTm = theHelix->Tm[a][b][c][d];
            bestReg[0] = a;
            bestReg[1] = b;
            bestReg[2] = c;
            bestReg[3] = d;
        }
        else
        {
            if (theHelix->Tm[a][b][c][d] >= secondBestTm)
            {
                secondBestTm = theHelix->Tm[a][b][c][d];
                secondBestReg[0] = a;
                secondBestReg[1] = b;
                secondBestReg[2] = c;
                secondBestReg[3] = d;
            }
        }
        
        if (theHelix->numPep == 1)
        {
            if (theHelix->Tm[a][b][c][d] >= theHelix->CCTm)
            {
                //cout << "Setting CCTm for homotrimer" << endl;
                theHelix->CCTm = theHelix->Tm[a][b][c][d];
                theHelix->CCRegister[0] = a;
                theHelix->CCRegister[1] = b;
                theHelix->CCRegister[2] = c;
            }
        }
        
        if ((theHelix->numPep == 2) && ((a != b) || (a != c) || (b != c)))
        {
            
            if (theHelix->Tm[a][b][c][d] >= theHelix->CCTm)
            {
                //cout << "This is an A2B register: " << n << " " <<bestReg[0]<<bestReg[1]<<bestReg[2] << endl;
                theHelix->CCTm = theHelix->Tm[a][b][c][d];
                theHelix->CCRegister[0] = a;
                theHelix->CCRegister[1] = b;
                theHelix->CCRegister[2] = c;
                theHelix->CCRegister[3] = d;
            }
        }
            
        if ((theHelix->numPep == 3) && ((a != b) && (a != c) && (b != c)))
        {
            
            if (theHelix->Tm[a][b][c][d] >= theHelix->CCTm)
            {
                //cout << "This is an ABC register: " << n << " " <<bestReg[0]<<bestReg[1]<<bestReg[2] << endl;
                theHelix->CCTm = theHelix->Tm[a][b][c][d];
                theHelix->CCRegister[0] = a;
                theHelix->CCRegister[1] = b;
                theHelix->CCRegister[2] = c;
                theHelix->CCRegister[3] = d;
            }
        }
        
        
    } // end for abcd
    
    
    // move all values to the helix parameters...
    theHelix->HighTm = maxTm;
    theHelix->bestRegister[0] = bestReg[0];
    theHelix->bestRegister[1] = bestReg[1];
    theHelix->bestRegister[2] = bestReg[2];
    theHelix->bestRegister[3] = bestReg[3];
    theHelix->bestPropensity = theHelix->Propensity[theHelix->bestRegister[0]][theHelix->bestRegister[1]][theHelix->bestRegister[2]][theHelix->bestRegister[3]];
    theHelix->bestPairwise = theHelix->PairWise[theHelix->bestRegister[0]][theHelix->bestRegister[1]][theHelix->bestRegister[2]][theHelix->bestRegister[3]];
    theHelix->secTm = secondBestTm;
    theHelix->secRegister[0] = secondBestReg[0];
    theHelix->secRegister[1] = secondBestReg[1];
    theHelix->secRegister[2] = secondBestReg[2];
    theHelix->secRegister[3] = secondBestReg[3];
    theHelix->specificity = maxTm - secondBestTm;
    // This sets the deviation as between the experimental Tm and the correct composition / register found here.
    // This was originally created to allow known systems that have POGn peptides with high Tm to be scored for
    // their ABC or A2B helix. However it likely introduces some problems in optimization.
    // theHelix->deviation = theHelix->CCTm - theHelix->expTm;
    // // // // // // // // // // // //
    // Using maxTm is more rigorous but has its own set of problems.
    // theHelix->deviation = maxTm - theHelix->expTm;
    // Third option uses maxTm and the differenct between maxTm and CCTm for the difference score
    
    if ((theHelix->CCRegister[0] == theHelix->bestRegister[0]) && (theHelix->CCRegister[1] == theHelix->bestRegister[1]) && (theHelix->CCRegister[2] == theHelix->bestRegister[2]))
    {
        // The best register and expected register (CCregisters) are the same. This does not consider non-canonical registers.
        if (theHelix->expTm == -10)
        {
            // an experimental Tm of -10 indicates that no transition was observed experimentally. Only penalize Tm's higher than 10C.
            if (maxTm <= 10) theHelix->deviation = 0; else theHelix->deviation = maxTm - 10;
        }
        else
            theHelix->deviation = maxTm - theHelix->expTm;
    }
    else
    {
        // The best register and the expected register are not the same.
        // Score as above but with additional penalty for wrong register.
        if (theHelix->expTm == -10)
        {
            if (maxTm <= 10) theHelix->deviation = 0;
            if (maxTm >  10) theHelix->deviation = maxTm - 10;
        }
        else
        {
            theHelix->deviation = theHelix->CCTm - theHelix->expTm;
            if (theHelix->deviation < 0)
            {
                // cout << "Deviation was: " << theHelix->deviation << ". Penalty of " << -0.5*(abs(theHelix->CCTm - maxTm)) << " applied. Final Deviation: ";
                theHelix->deviation -= 0.5*(abs(theHelix->CCTm - maxTm));
                //cout <<  theHelix->deviation << endl;
            }
            else
            {
                // cout << "Deviation was: " << theHelix->deviation << ". Penalty of " << 0.5*(abs(theHelix->CCTm - maxTm)) << " applied. Final Deviation: ";
                theHelix->deviation += 0.5*(abs(theHelix->CCTm - maxTm));
                //cout <<  theHelix->deviation << endl;
            }
        }
    }
    

    
    
}

void ScoreLibrary (short start, short stop, parameterType parameters, TripleHelix * Lib)
{
    short n;
    
    for (n=start;n<stop;n++)
    {
        ScoreHelix(parameters, &Lib[n]);
    }
}

void pause (double wait)
{
    time_t t;
    t = clock();
    while ((double(clock() - t)) < wait)
    {
        // do nothing
    }
}

short readLibrary (TripleHelix * Lib, string Lib_Name)
{
    short x, y, n;
    
    short TotalHelices = 0;
    string killString;
    string seqDate;
    char peptideInput[100]; // holds the peptide sequences as they are input
    ifstream seq_input(Lib_Name);
    
    if (!seq_input.is_open())
    {
        cout << "We couldn't open the seq_input.txt file." << endl;
        return 0;
    }
    getline (seq_input, seqDate);
    cout << "Sequence Library: " << seqDate << endl;
    seq_input >> TotalHelices; // indicates the total number of helices to be expected from the text file input
    
    x = 0;
    n = 0; // n tracks the triple helix we are evaluating.
    for (n=0; n<(TotalHelices); n++)
    {
        x = 0;
        seq_input >> Lib[n].numPep;
        while (Lib[n].numPep == 0)
        {
            getline(seq_input, killString);
            seq_input >> Lib[n].numPep;
            //cout << n << " " << x << " " << killString << endl;
            x++;
            if (x>50)
            {
                cout << "Problem reading " << Lib_Name << ". Likely problem: is either 1) indicated number of helices exceeds actual number of helices or 2) with commenting." << endl;
                return 0;
            }
        }
        x = 0;
        if ((Lib[n].numPep < 1) || (Lib[n].numPep>3))
        {
            cout << "Number of unique peptides in a helix must be 1-3. Value read was Library[" << n << "].numPep = " << Lib[n].numPep << ". Stopping now." <<endl;
            Lib[n-1].dissect();
            Lib[n].dissect();
            return 0;
        }

        seq_input >> Lib[n].numAA;
        
        if ((Lib[n].numAA < 21) || (Lib[n].numAA>48))
        {
            cout << "Number of amino acids in the peptide must be 21-48. Value read was Library[" << n << "].numAA = " << Lib[n].numAA  << ". Stopping now." <<endl;
            Lib[n-1].dissect();
            Lib[n].dissect();
            return 0;
        }
        
        seq_input >> Lib[n].Nterm;
        seq_input >> Lib[n].Cterm;
        seq_input >> Lib[n].expTm;
                
        for (x=0; x<Lib[n].numPep; x++)
        {
            for (y=0;y<Lib[n].numAA;y++)
            {
                seq_input >> peptideInput[y];
                Lib[n].sequences[x][y] = toupper(peptideInput[y]);
            }
        }
    }
    seq_input.close();
    
    // Determine Xaa, Yaa and Gly positions (reptition)
    for (n=0; n<(TotalHelices); n++) Lib[n].determine_reptition();
        
    
    
    
    
    
    
    return TotalHelices;
}

// // // // // // // // // // //
// MAIN STARTS HERE!
// // // // // // // // // // //
int main (int argc, const char * argv[])
{
    cout << "-------------------------------------------" << endl;
    cout << "SCEPTTr" << endl;
    cout << "v1.2 BETA 2021-09-23" << endl;
    cout << "Only canonical compositions/registers examined!" << endl;
    cout << "Based in part on previously published" << endl;
    cout << "SCEPTTr 1.0 and 1.1" <<endl;
    cout << "-------------------------------------------" << endl;

    // cout << "64 == " << char(64) << ", 65 == " << char(65) << endl;
    // 64 == '@', 65 == 'A'
    
    clock_t time;
    
    
    
    TripleHelix Library[500]; // Holds results for all triple helices in the library. Maybe make this dynamic based on read input.
    TripleHelix userHelix;
    TripleHelix userLib[100];
    
    short TotalHelices=0; // This is the total number of helices that will be evaluated from the read-in text file.
    short totalUserHelices = 0;
    
    // various counters
    short n=0, x=0, y=0;
    //short xCount, yCount, zCount;
    
    double sumSquaredDev = 0;
    double NewSumSquaredDev = 0;
    
    
    
    
    // This is the amount that parameters are changed during each optimization round.
    // A value of approximately 0.02 seems to be a good trade off between speed and effective optimization.
    // A value above 0.10 doesn't work all that well.
    double delta = 0.1;
    
    // This is the maximum number of improvement rounds that will be tried before stopping.
    short maxRounds = 25;
    
    // During the optimization this is the largest a value can deviate from experimental parameters.
    double maxDev = 2.0; // negative number should make it impossible to make changes...
    
    // This struct holds all the parameters that are read from file.
    parameterType parameters;
    parameterType countInteractions;
    //parameterType expVal;
    
    // These booleans are used to decide if a parameter should
    // be included in optimization or not.
    
    
    //bool goodPeptide;
    short worstHelix = -1;
    double worstDeviation = -1;
    bool done = false;
    short round = 0;
    bool improved = false;
    short useCase = 3;
    
    // // // // // // // // // // //
    // READ INITIAL PARAMETERS HERE
    // // // // // // // // // // //
    parameters = ReadParameters();
    // DisplayParameters(parameters);
    
    
    
    // // // // // // // // // // //
    // READ PEPTIDE SEQUENCES HERE
    // // // // // // // // // // //
    
    for (x=0;x<500;x++) Library[x].initializeAll();
    
    
    TotalHelices = readLibrary(Library, "seq_input.txt");
    
    cout << "TotalHelices in training library = " << TotalHelices << endl;
    if (TotalHelices == 0)
    {
        cout << "TotalHelices in training library = " << TotalHelices << ". Stopping." << endl;
        return 0;
    }
    
    // cout << "Do you want to (0) optimize parameters against the existing peptide library, (1) manually enter the parameters for a new helix or (2) evaluate user_lib.txt?" << endl;
    cout << "Do you want to (1) manually enter the parameters for a new helix or (2) evaluate user_lib.txt?" << endl;
    while ((useCase != 0) && (useCase != 1) && (useCase != 2))
    {
        cin >> useCase;
    }
    
    //useCase = 1;
    
    // useCase = 1;
    // Ask user for sequence information
    if (useCase == 1)
    {
        cout << "How many distinct peptides are in your helix? (1) Homotrimer, (2) A2B Heterotrimer, or (3) ABC Heterotrimer?" << endl;
        cin >> userHelix.numPep;
        while ((userHelix.numPep < 1) || (userHelix.numPep > 3)) cin >> userHelix.numPep;
        
        cout << "How many amino acids are in each peptide? Acceptable range 18-60." << endl;
        cin >> userHelix.numAA;
        while ((userHelix.numAA < 21) || (userHelix.numAA >48)) cin >> userHelix.numAA;
            
        cout << "How is the N-terminus functionalized? N = free amine, Ac = Acetylated" << endl;
        cin >> userHelix.Nterm;
        
        cout << "How is the C-terminus functionalized? C = Carboxylic acid, Am = Amidated" << endl;
        cin >> userHelix.Cterm;
        
        cout << "Type the sequence of the first peptide using single letter amino acid code." << endl;
        cin >> userHelix.sequences[0];
        
        if ((userHelix.numPep == 2) || (userHelix.numPep == 3))
        {
            cout << "Type the sequence of the second peptide using single letter amino acid code." << endl;
            cin >> userHelix.sequences[1];
        }
        if (userHelix.numPep == 3)
        {
            cout << "Type the sequence of the third peptide using single letter amino acid code." << endl;
            cin >> userHelix.sequences[2];
        }
    
        // Convert all amino acids to upper case
        for (x=0; x<userHelix.numPep; x++)
        {
            for (y=0;y<userHelix.numAA;y++)
            {
                userHelix.sequences[x][y] = toupper(userHelix.sequences[x][y]);
            }
        }
        userHelix.determine_reptition();
    }
    
    
    // // // // // // // // //
    // At this point we have declared our variables, read our parameters and sequences. Next evaluate the training library.
    // // // // // // // // //
    
    // Count of all Xaa, Yaa and pairwise interactions possible in canonical registers for this library
    // Zero interactions
    for (y=0;y<27;y++)
    {
        countInteractions.propensityX[y] = 0;
        countInteractions.propensityY[y] = 0;
        for (x=0;x<27;x++)
        {
            countInteractions.lateral[y][x] = 0;
            countInteractions.axial[y][x] = 0;
        }
    }
    
    // Count interactions
    short a,b,c;
    for (n=0; n<TotalHelices; n++) // doesn't look at the user's helix.
    {
        for (a=0;a<Library[n].numPep;a++)
        {
            for (x=0;x<Library[n].numAA;x++)
            {
                if (Library[n].isXaa(x)) countInteractions.propensityX[short(Library[n].sequences[a][x])-64]++;
                if (Library[n].isYaa(x)) countInteractions.propensityY[short(Library[n].sequences[a][x])-64]++;
            }
            for (b=0;b<Library[n].numPep;b++)for (c=0;c<Library[n].numPep;c++)
            {
                for (x=0;x<Library[n].numAA;x++)
                {
                    if (Library[n].isYaa(x))
                    {
                        if ((x+2) < Library[n].numAA) countInteractions.axial[short(Library[n].sequences[a][x])-64][short(Library[n].sequences[b][x+2])-64]++;
                        if ((x+2) < Library[n].numAA) countInteractions.axial[short(Library[n].sequences[b][x])-64][short(Library[n].sequences[c][x+2])-64]++;
                        if ((x+5) < Library[n].numAA) countInteractions.axial[short(Library[n].sequences[c][x])-64][short(Library[n].sequences[a][x+5])-64]++;
                    
                        if (x>1) countInteractions.lateral[short(Library[n].sequences[a][x])-64][short(Library[n].sequences[b][x-1])-64]++;
                        if (x>1) countInteractions.lateral[short(Library[n].sequences[b][x])-64][short(Library[n].sequences[c][x-1])-64]++;
                        if ((x+2) < Library[n].numAA) countInteractions.lateral[short(Library[n].sequences[c][x])-64][short(Library[n].sequences[a][x+2])-64]++;
                    }
                }
            }
        }
    }
    
    // Set optimization by counts.
    for (y=0;y<27;y++)
    {
        if (countInteractions.propensityX[y] > 25) parameters.optPropX[y] = true;
        if (countInteractions.propensityY[y] > 25) parameters.optPropY[y] = true;
        for (x=0;x<27;x++)
        {
            if (countInteractions.lateral[y][x] > 25) parameters.optLat[y][x] = true;
            if (countInteractions.axial[y][x] > 25) parameters.optAxial[y][x] = true;
        }
    }
    // but force Xaa Pro ax & lat interactions to remain at 0. and Yaa Hyp ax & lat interactions to remain at 0.
    for (x=0;x<27;x++)
    {
        parameters.optPropX[short('P'-64)] = false;
        parameters.optPropY[short('O'-64)] = false;
        
        parameters.optLat[short('O')-64][x] = false;
        parameters.optAxial[short('O')-64][x] = false;
        parameters.optAxial[short('P')-64][x] = false; // This was creating results here that don't make much sense imo.
        
        parameters.optLat[x][short('P')-64] = false;
        parameters.optAxial[x][short('P')-64] = false;
        
    }
    
    //
    // For Axial & Lateral interactions: 0-25 counts = poorly modeled. 26-50 counts = moderately modeled. 51+ counts = well modeled.
    //
    // Optionally display interaction counts.
    /*
    for (x=1;x<27;x++)
    {
        cout << char(x + 64) << "\t" << countInteractions.propensityX[x] << "\t" << countInteractions.propensityY[x] << endl;
    }
    cout << endl;
    cout << "Axial Counts:" << endl;
    for (y=0;y<27;y++)
    {
        cout << char(y + 64);
        for (x=1;x<27;x++)
        {
            if (y==0) cout << "\t" << char(x+64);
            if (y!=0) cout << "\t" << countInteractions.axial[y][x];
        }
        cout << endl;
    }
    cout << endl;
    cout << "Lateral Counts:" << endl;
    
    for (y=0;y<27;y++)
    {
        cout << char(y + 64);
        for (x=1;x<27;x++)
        {
            if (y==0) cout << "\t" << char(x+64);
            if (y!=0) cout << "\t" << countInteractions.lateral[y][x];
        }
        cout << endl;
    }
    */
    
    // Optimize against the training library "seq_input.txt"
    if (useCase == 0)
    {
        // // // // // // // //
        // First scoring of library with default parameters.
        // // // // // // // //
        sumSquaredDev = 0;
        worstHelix = -1;
        worstDeviation = -1;
        double sumDeviation = 0;
        
        thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
        thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
        th1.join();
        th2.join();
        for (n=0; n<(TotalHelices); n++)
        {
            sumDeviation += Library[n].deviation;
            sumSquaredDev += (Library[n].deviation * Library[n].deviation);
            // if (abs(Library[n].deviation) > 1) sumSquaredDev += (abs(Library[n].deviation)-1);
            
            if (abs(Library[n].deviation) > abs(worstDeviation))
            {
                worstDeviation = Library[n].deviation;
                worstHelix = n;
            }
            if (abs(Library[n].deviation) > 9)
            {
                cout << "Helix Number: " << n << endl;
                Library[n].dissect();
                Library[n].userOutput();
                cout << endl;
            }
            // Library[n].dissect(); // show everything!
        } // end for loop of helix evaluation.
        
        /*
        cout << "Total Helices in Standard Library = " << TotalHelices << endl;
        Library[TotalHelices-1].dissect();
        Library[TotalHelices-1].userOutput();
        cout << endl;
        cout << "You are the worst Burr:" << endl;
        cout << "Helix Number: " << worstHelix << endl;
        Library[worstHelix].dissect();
        cout << endl;
        Library[worstHelix].userOutput();
        */
    
        cout << "Initial sumDeviation = " << sumDeviation << ". Average sumDeviation = " << sumDeviation / TotalHelices << endl;
        cout << "Initial sumSquaredDev = " << sumSquaredDev << ". Average sumSquaredDev = " << sumSquaredDev / TotalHelices << endl;
        cout << "Starting Optimization." << endl;
        cout << endl;
        //cout << endl << "Last Helix in Library: " << TotalHelices << endl;
        
        
        
        NewSumSquaredDev = 0;
        done = false;
        round = 0;
        improved = false;
        bool improvedRound = false;
        done = true; // ie don't make changes! Comment this out to allow optimization.
        
        cout << "Maximum Deviation from Experimental Values = " << maxDev << endl;
        cout << "delta (change per test) = " << delta << endl;
        cout << "Max Rounds = " << maxRounds << endl;
        cout << endl;
        
        
        
        
        time = clock();
        
        
        while (not done)
        {
        // // // // // // // //
        // Now optimize parameters where possible.
        // // // // // // // //
        {
        // LENGTH PARAMETERS //
        /* for now, don't try to optimize these since we don't have a lot of different sized peptides in our library. Mostly just 24 & 30 aa peptides.
            // // //
            // A  //
            // // //
            parameters.A -= 0.01;
            NewSumSquaredDev = 0;
            for (n=0; n<TotalHelices; n++)
            {
                Library[n] = ScoreHelix (parameters, Library[n]);
                NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
            }
            if (NewSumSquaredDev < sumSquaredDev)
            {
                // keep this new parameter and move on.
                cout << "A-length. New SSD: " << NewSumSquaredDev << ". Old SSD: " << sumSquaredDev << endl;
                sumSquaredDev = NewSumSquaredDev;
                improved = true;
            }
            else
            {
                parameters.A += 0.02;
                NewSumSquaredDev = 0;
                for (n=0; n<TotalHelices; n++)
                {
                    Library[n] = ScoreHelix (parameters, Library[n]);
                    NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                }
                if (NewSumSquaredDev < sumSquaredDev)
                {
                    // keep this new parameter and move on.
                    cout << "A-length. New SSD: " << NewSumSquaredDev << ". Old SSD: " << sumSquaredDev << endl;
                    sumSquaredDev = NewSumSquaredDev;
                    improved = true;
                }
                else
                {
                    // neither change resulted in an improvement. Go back to original parameter.
                    parameters.A -= 0.01;
                }
            }
        
            // // //
            // B  //
            // // //
            parameters.B -= 0.002;
            NewSumSquaredDev = 0;
            for (n=0; n<TotalHelices; n++)
            {
                Library[n] = ScoreHelix (parameters, Library[n]);
                NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
            }
            if (NewSumSquaredDev < sumSquaredDev)
            {
                // keep this new parameter and move on.
                cout << "B-length. New SSD: " << NewSumSquaredDev << ". Old SSD: " << sumSquaredDev << endl;
                sumSquaredDev = NewSumSquaredDev;
                improved = true;
            }
            else
            {
                parameters.B += 0.004;
                NewSumSquaredDev = 0;
                for (n=0; n<TotalHelices; n++)
                {
                    Library[n] = ScoreHelix (parameters, Library[n]);
                    NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                }
                if (NewSumSquaredDev < sumSquaredDev)
                {
                    // keep this new parameter and move on.
                    cout << "B-length. New SSD: " << NewSumSquaredDev << ". Old SSD: " << sumSquaredDev << endl;
                    sumSquaredDev = NewSumSquaredDev;
                    improved = true;
                }
                else
                {
                    // neither change resulted in an improvement. Go back to original parameter.
                    parameters.B -= 0.002;
                }
            }
        
            // // //
            // C  //
            // // //
            parameters.C -= 0.0002;
            NewSumSquaredDev = 0;
            for (n=0; n<TotalHelices; n++)
            {
                Library[n] = ScoreHelix (parameters, Library[n]);
                NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
            }
            if (NewSumSquaredDev < sumSquaredDev)
            {
                // keep this new parameter and move on.
                cout << "C-length. New SSD: " << NewSumSquaredDev << ". Old SSD: " << sumSquaredDev << endl;
                sumSquaredDev = NewSumSquaredDev;
                improved = true;
            }
            else
            {
                parameters.C += 0.0004;
                NewSumSquaredDev = 0;
                for (n=0; n<TotalHelices; n++)
                {
                    Library[n] = ScoreHelix (parameters, Library[n]);
                    NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                }
                if (NewSumSquaredDev < sumSquaredDev)
                {
                    // keep this new parameter and move on.
                    cout << "C-length. New SSD: " << NewSumSquaredDev << ". Old SSD: " << sumSquaredDev << endl;
                    sumSquaredDev = NewSumSquaredDev;
                    improved = true;
                }
                else
                {
                    // neither change resulted in an improvement. Go back to original parameter.
                    parameters.C -= 0.0002;
                }
            }
            */
        } // hidden length optimization which we are not doing at the moment.
            
        for (x=0;x<27;x++)
        {
            if (parameters.optPropX[x])
            {
                parameters.propensityX[x] -= delta;
                NewSumSquaredDev = 0;
                if (parameters.propensityX[x] >= (parameters.exPropensityX[x] - maxDev)) // Only test this optimization if we are in range (2) of experimental values.
                {
                    thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
                    thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
                    th1.join();
                    th2.join();
                    for (n=0; n<TotalHelices; n++)
                    {
                        //ScoreHelix(parameters, &Library[n]);
                        //cout << n << " ";
                        NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                        //if (abs(Library[n].deviation) > 1) NewSumSquaredDev += (abs(Library[n].deviation)-1);
                    }
                    //cout << endl;
                    if (NewSumSquaredDev < sumSquaredDev)
                    {
                        // keep this new parameter and move on.
                        cout << "Xaa" << char(x+64) << " adjusted to " << parameters.propensityX[x] << ". New SSDev = " << NewSumSquaredDev << endl;
                        sumSquaredDev = NewSumSquaredDev;
                        improved = true;
                        improvedRound = true;
                    }
                }
                if (not improved)
                {
                    parameters.propensityX[x] += delta;
                    parameters.propensityX[x] += delta;
                    NewSumSquaredDev = 0;
                    if (parameters.propensityX[x] <= (parameters.exPropensityX[x] + maxDev)) // Only test this optimization if we are in range (2) of experimental values.
                    {
                        thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
                        thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
                        th1.join();
                        th2.join();
                        for (n=0; n<TotalHelices; n++)
                        {
                            //ScoreHelix(parameters, &Library[n]);
                            NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                        }
                        if (NewSumSquaredDev < sumSquaredDev)
                        {
                            // keep this new parameter and move on.
                            cout << "Xaa" << char(x+64) << " adjusted to " << parameters.propensityX[x] << ". New SSDev = " << NewSumSquaredDev << endl;
                            sumSquaredDev = NewSumSquaredDev;
                            improved = true;
                            improvedRound = true;
                        }
                    }
                }
                if (not improved)
                {
                    // neither change resulted in an improvement. Go back to original parameter.
                    parameters.propensityX[x] -= delta;
                    // cout << "Xaa" << x << ".  No changes resulted in an improvement." << endl;
                }
            }
            improved = false;
            
            if (parameters.optPropY[x])
            {
                parameters.propensityY[x] -= delta;
                NewSumSquaredDev = 0;
                if (parameters.propensityY[x] >= (parameters.exPropensityY[x] - maxDev)) // Only test this optimization if we are in range (2) of experimental values.
                    {
                    thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
                    thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
                    th1.join();
                    th2.join();
                    for (n=0; n<TotalHelices; n++)
                    {
                        //ScoreHelix(parameters, &Library[n]);
                        NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                        //if (abs(Library[n].deviation) > 1) NewSumSquaredDev += (abs(Library[n].deviation)-1);
                    }
                    if (NewSumSquaredDev < sumSquaredDev)
                    {
                        // keep this new parameter and move on.
                        cout << "Yaa" << char(x+64) << " adjusted to " << parameters.propensityY[x] << ". New SSDev = " << NewSumSquaredDev << endl;
                        sumSquaredDev = NewSumSquaredDev;
                        improved = true;
                        improvedRound = true;
                    }
                }
                if (not improved)
                {
                    parameters.propensityY[x] += delta;
                    parameters.propensityY[x] += delta;
                    NewSumSquaredDev = 0;
                    if (parameters.propensityY[x] <= (parameters.exPropensityY[x] + maxDev)) // Only test this optimization if we are in range (2) of experimental values.
                    {
                        thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
                        thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
                        th1.join();
                        th2.join();
                        for (n=0; n<TotalHelices; n++)
                        {
                            //ScoreHelix(parameters, &Library[n]);
                            NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                            //if (abs(Library[n].deviation) > 1) NewSumSquaredDev += (abs(Library[n].deviation)-1);
                        }
                        if (NewSumSquaredDev < sumSquaredDev)
                        {
                            // keep this new parameter and move on.
                            cout << "Yaa" << char(x+64) << " adjusted to " << parameters.propensityY[x] << ". New SSDev = " << NewSumSquaredDev << endl;
                            sumSquaredDev = NewSumSquaredDev;
                            improved = true;
                            improvedRound = true;
                        }
                    }
                }
                if (not improved)
                {
                    // neither change resulted in an improvement. Go back to original parameter.
                    parameters.propensityY[x] -= delta;
                    // cout << "Yaa" << x << ".  No changes to this parameter resulted in an improvement." << endl;
                }
            }
            improved = false;
            
            for (y=0;y<27;y++)
            {
                if (parameters.optAxial[x][y])
                {
                    parameters.axial[x][y] -= delta;
                    NewSumSquaredDev = 0;
                    if (parameters.axial[x][y] >= (parameters.exAxial[x][y] - maxDev)) // Only test this optimization if we are in range (2) of experimental values.
                    {
                        thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
                        thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
                        th1.join();
                        th2.join();
                        for (n=0; n<TotalHelices; n++)
                        {
                            // ScoreHelix(parameters, &Library[n]);
                            NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                        }
                        if (NewSumSquaredDev < sumSquaredDev)
                        {
                            // keep this new parameter and move on.
                            cout << "axial" << char(x+64) << "," << char(y+64) << " adjusted to " << parameters.axial[x][y] << ".  New SSDev =  " << NewSumSquaredDev << endl;
                            sumSquaredDev = NewSumSquaredDev;
                            improved = true;
                            improvedRound = true;
                        }
                    }
                    if (not improved)
                    {
                        parameters.axial[x][y] += delta;
                        parameters.axial[x][y] += delta;
                        NewSumSquaredDev = 0;
                        if (parameters.axial[x][y] <= (parameters.exAxial[x][y] + maxDev)) // Only test this optimization if we are in range (2) of experimental values.
                            {
                                thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
                                thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
                                th1.join();
                                th2.join();
                                for (n=0; n<TotalHelices; n++)
                                {
                                    // ScoreHelix(parameters, &Library[n]);
                                    NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                                }
                                if (NewSumSquaredDev < sumSquaredDev)
                                {
                                    // keep this new parameter and move on.
                                    cout << "axial" << char(x+64) << "," << char(y+64) << " adjusted to " << parameters.axial[x][y] << ".  New SSDev =  " << NewSumSquaredDev << endl;
                                    sumSquaredDev = NewSumSquaredDev;
                                    improved = true;
                                    improvedRound = true;
                                }
                        }
                    }
                    if (not improved)
                    {
                        // neither change resulted in an improvement. Go back to original parameter.
                        parameters.axial[x][y] -= delta;
                        // cout << "axial" << x << "," << y << ".  No changes to this parameter resulted in an improvement." << endl;
                    }
                }
                improved = false;
                
                if (parameters.optLat[x][y])
                {
                    parameters.lateral[x][y] -= delta;
                    NewSumSquaredDev = 0;
                    if (parameters.lateral[x][y] >= (parameters.exLateral[x][y] - maxDev)) // Only test this optimization if we are in range (2) of experimental values.
                        {
                        // optimize this parameter (test slightly lower & slightly higher, pick the best)
                        thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
                        thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
                        th1.join();
                        th2.join();
                        for (n=0; n<TotalHelices; n++)
                        {
                            //ScoreHelix(parameters, &Library[n]);
                            NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                            //if (abs(Library[n].deviation) > 1) NewSumSquaredDev += (abs(Library[n].deviation)-1);
                        }
                        if (NewSumSquaredDev < sumSquaredDev)
                        {
                            // keep this new parameter and move on.
                            cout << "lateral" << char(x+64) << "," << char(y+64) << " adjusted to " << parameters.lateral[x][y] << ".  New SSDev =  " << NewSumSquaredDev << endl;
                            sumSquaredDev = NewSumSquaredDev;
                            improved = true;
                            improvedRound = true;
                        }
                    }
                    if (not improved)
                    {
                        parameters.lateral[x][y] += delta;
                        parameters.lateral[x][y] += delta;
                        NewSumSquaredDev = 0;
                        if (parameters.lateral[x][y] <= (parameters.exLateral[x][y] + maxDev)) // Only test this optimization if we are in range (2) of experimental values.
                        {
                            thread th1(ScoreLibrary, 0, TotalHelices/2, parameters, Library);
                            thread th2(ScoreLibrary, TotalHelices/2, TotalHelices, parameters, Library);
                            th1.join();
                            th2.join();
                            for (n=0; n<TotalHelices; n++)
                            {
                                //ScoreHelix(parameters, &Library[n]);
                                NewSumSquaredDev += (Library[n].deviation * Library[n].deviation);
                                //if (abs(Library[n].deviation) > 1) NewSumSquaredDev += (abs(Library[n].deviation)-1);
                            }
                            if (NewSumSquaredDev < sumSquaredDev)
                            {
                                // keep this new parameter and move on.
                                cout << "lateral" << char(x+64) << "," << char(y+64) << " adjusted to " << parameters.lateral[x][y] << ".  New SSDev =  " << NewSumSquaredDev << endl;
                                sumSquaredDev = NewSumSquaredDev;
                                improved = true;
                                improvedRound = true;
                            }
                        }
                    }
                    if (not improved)
                    {
                        // neither change resulted in an improvement. Go back to original parameter.
                        parameters.lateral[x][y] -= delta;
                        // cout << "lateral" << x << "," << y << ".  No changes to this parameter resulted in an improvement." << endl;
                    }
                }
                improved = false;
            }
        }
            
            if (not improvedRound) done = true;
            improved = false;
            improvedRound = false;
            round++;
            if (round >= maxRounds) done = true;
            cout << "End round #" << round << ". Avg of SSDev = " << sumSquaredDev / TotalHelices << endl << endl;
        } // end while loop
        time = clock() - time;
        
        cout << "Total time required: " << double(time)/CLOCKS_PER_SEC << endl;
        cout << endl;
        
        // End Game Items:
        cout << endl;
        cout << "Sum of Squared Deviation = " << sumSquaredDev << endl;
        cout << "Average of Squared Deviations = " << sumSquaredDev / TotalHelices << endl;
        
        // Write out all new optimized parameters.
        WriteParameters(parameters);
        
        ofstream output("A3.txt");
        if (!output.is_open())
        {
            cout << "Failed to open file." << endl;
            return 0;
        }
        // Write to the file
        output << "n ExpTm A3 HighTm Dev" << endl;
        for (n=0;n<TotalHelices;n++)
        {
            if (Library[n].numPep == 1) output << n << " " << Library[n].expTm << " " << Library[n].CCTm << " " << Library[n].HighTm << " " << Library[n].deviation << endl;
        }
        // Close the file
        output.close();
        
        output.open("A2B.txt");
        if (!output.is_open())
        {
            cout << "Failed to open file." << endl;
            return 0;
        }
        // Write to the file
        output << "n ExpTm A2B HighTm Dev" << endl;
        for (n=0;n<TotalHelices;n++)
        {
            if (Library[n].numPep == 2) output << n << " " << Library[n].expTm << " " << Library[n].CCTm << " " << Library[n].HighTm << " " << Library[n].deviation << endl;
        }
        // Close the file
        output.close();
        
        output.open("ABC.txt");
        if (!output.is_open())
        {
            cout << "Failed to open file." << endl;
            return 0;
        }
        // Write to the file
        output << "n ExpTm ABC HighTm Dev" << endl;
        for (n=0;n<TotalHelices;n++)
        {
            if (Library[n].numPep == 3) output << n << " " << Library[n].expTm << " " << Library[n].CCTm << " " << Library[n].HighTm << " " << Library[n].deviation << endl;
        }
        // Close the file
        output.close();
    }
    
    // Score user peptide
    if (useCase == 1)
    {
        // tell user about their peptide
        ScoreHelix(parameters, &userHelix);
        userHelix.dissect();
        userHelix.userOutput();
        
        short poorPrediction = 0;
        short poorAxList[27][27];
        short poorLatList[27][27];
        for (x=0;x<27;x++)for(y=0;y<27;y++)
        {
            poorAxList[x][y] = 0;
            poorLatList[x][y] = 0;
        }
        short lowConfidenceCut = 25;
        n = TotalHelices;
        for (a=0;a<userHelix.numPep;a++) for (b=0;b<userHelix.numPep;b++)for (c=0;c<userHelix.numPep;c++)
        {
            for (x=0;x<userHelix.numAA;x++)
            {
                if (userHelix.isYaa(x))
                {
                    if ((x+2) < userHelix.numAA)
                    {
                        if (countInteractions.axial[short(userHelix.sequences[a][x])-64][short(userHelix.sequences[b][x+2])-64] < lowConfidenceCut)
                        {
                            poorPrediction++;
                            poorAxList[short(userHelix.sequences[a][x])-64][short(userHelix.sequences[b][x+2])-64]++;
                        }
                    }
                    if ((x+2) < userHelix.numAA)
                    {
                        if (countInteractions.axial[short(userHelix.sequences[b][x])-64][short(userHelix.sequences[c][x+2])-64] < lowConfidenceCut)
                        {
                            poorPrediction++;
                            poorAxList[short(userHelix.sequences[b][x])-64][short(userHelix.sequences[c][x+2])-64]++;
                        }
                    }
                    if ((x+5) < Library[n].numAA)
                    {
                        if (countInteractions.axial[short(userHelix.sequences[c][x])-64][short(userHelix.sequences[a][x+5])-64] < lowConfidenceCut)
                        {
                            poorPrediction++;
                            poorAxList[short(userHelix.sequences[c][x])-64][short(userHelix.sequences[a][x+5])-64]++;
                        }
                    }
                    if (x>1)
                    {
                        if (countInteractions.lateral[short(userHelix.sequences[a][x])-64][short(userHelix.sequences[b][x-1])-64] < lowConfidenceCut)
                        {
                            poorPrediction++;
                            poorLatList[short(userHelix.sequences[a][x])-64][short(userHelix.sequences[b][x-1])-64]++;
                        }
                    }
                    if (x>1)
                    {
                        if (countInteractions.lateral[short(userHelix.sequences[b][x])-64][short(userHelix.sequences[c][x-1])-64] < lowConfidenceCut)
                        {
                            poorPrediction++;
                            poorLatList[short(userHelix.sequences[b][x])-64][short(userHelix.sequences[c][x-1])-64]++;
                        }
                    }
                    if ((x+2) < Library[n].numAA)
                    {
                        if (countInteractions.lateral[short(userHelix.sequences[c][x])-64][short(userHelix.sequences[a][x+2])-64] < lowConfidenceCut)
                        {
                            poorPrediction++;
                            poorLatList[short(userHelix.sequences[c][x])-64][short(userHelix.sequences[a][x+2])-64]++;
                        }
                    }
                    
                }
            }
        }
        
        cout << "Total Number of low confidence interactions in user helix: " << poorPrediction << "." << endl;
        cout << "Note: this count is based on all possible intereactions in all possible canonical registers / compositions." << endl;
        
        if (poorPrediction > 0)
        {
            cout << "Low confidence Axial Interactions (shown as Yaa,Xaa):" << endl;
            for (x=0;x<27;x++)for(y=0;y<27;y++)
            {
                if (poorAxList[x][y] > 0)
                {
                    cout << char(x+64) << "," << char(y+64) << ": " << poorAxList[x][y] << endl;
                }
            }
            cout << endl;
            cout << "Low confidence Lateral Interactions (shown as Yaa,Xaa):" << endl;
            for (x=0;x<27;x++)for(y=0;y<27;y++)
            {
                if (poorLatList[x][y] > 0)
                {
                    cout << char(x+64) << "," << char(y+64) << ": " << poorLatList[x][y] << endl;
                }
            }
        }
        cout << "Would you like to change any amino acids? (Y/N)" << endl;
        char Answer;
        short changeAA;
        short changePep;
        cin >> Answer;
        while ((Answer == 'Y') || (Answer == 'y'))
        {
            /*
            cout << "Peptide 0:" << endl;
            for (x=0;x<userHelix.numAA;x++)
            {
                cout << x/10;
                //cout << " ";
                //y++;
            }
            cout << endl;
            y=0;
            for (x=0;x<userHelix.numAA;x++)
            {
                if (y==10) y = 0;
                cout << y;
                y++;
            }
            cout << endl;
            for (x=0;x<userHelix.numAA;x++)
            {
                if (userHelix.sequences[0][x] == 'K') cout << "\x1b[1m\x1b[34m";
                if (userHelix.sequences[0][x] == 'R') cout << "\x1b[34m";
                if ((userHelix.sequences[0][x] == 'E') || (userHelix.sequences[0][x] == 'D')) cout << "\x1b[1m\x1b[31m";
                if ((userHelix.sequences[0][x] == 'F') || (userHelix.sequences[0][x] == 'Y') || (userHelix.sequences[0][x] == 'W')) cout << "\x1b[1m";
                if (userHelix.sequences[0][x] == 'Q') cout << "\x1b[1m\x1b[32m";
                cout << userHelix.sequences[0][x];
                cout << "\x1b[0m";
            }
            cout << endl;
            */
            for (a=0;a<userHelix.numPep;a++)
            {
                cout << endl;
                cout << "Peptide " << a << ": " << endl;
                for (x=0;x<userHelix.numAA;x++)
                {
                    cout << x/10;
                    //cout << " ";
                    //y++;
                }
                cout << endl;
                y=0;
                for (x=0;x<userHelix.numAA;x++)
                {
                    if (y==10) y = 0;
                    cout << y;
                    y++;
                }
                cout << endl;
                for (x=0;x<userHelix.numAA;x++)
                {
                    if (userHelix.sequences[a][x] == 'K') cout << "\x1b[1m\x1b[34m";
                    if (userHelix.sequences[a][x] == 'R') cout << "\x1b[34m";
                    if ((userHelix.sequences[a][x] == 'E') || (userHelix.sequences[a][x] == 'D')) cout << "\x1b[1m\x1b[31m";
                    if ((userHelix.sequences[a][x] == 'F') || (userHelix.sequences[a][x] == 'Y') || (userHelix.sequences[a][x] == 'W')) cout << "\x1b[1m";
                    if (userHelix.sequences[a][x] == 'Q') cout << "\x1b[1m\x1b[32m";
                    cout << userHelix.sequences[a][x];
                    cout << "\x1b[0m";
                }
                cout << endl;
            }
            
            cout << "Which peptide do you want to change? (0, 1 or 2)" << endl;
            cin >> changePep;
            cout << "Which # amino acid do you want to change?" << endl;
            cin >> changeAA;
            cout << "What should the new amino acid be? (single letter amino acid code)" << endl;
            cin >> userHelix.sequences[changePep][changeAA];
            ScoreHelix(parameters, &userHelix);
                           userHelix.userOutput();
            cout << endl;
            cout << "Another change?" << endl;
            cin >> Answer;
        }
    }
        
    // Score "user_lib.txt"
    if (useCase == 2)
    {
        totalUserHelices = readLibrary(userLib, "user_lib.txt");
        ScoreLibrary(0, totalUserHelices, parameters, userLib);
        cout << "totalUserHelices = " << totalUserHelices << endl;
        if (TotalHelices == 0)
        {
            cout << "totalUserHelices = " << totalUserHelices << ". Stopping." << endl;
            return 0;
        }
        cout << endl;
        for (n=0; n<totalUserHelices; n++)
        {
            // userLib[n].dissect();
            cout << "User Helix #" << n+1 << endl;
            userLib[n].userOutput();
        }
    }
    
} // end main()
