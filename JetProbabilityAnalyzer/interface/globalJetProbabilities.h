// Joshua Hardenbrook 18-11-2015
// 1. encapsulates the relative tag efficiency given some selection
// 2. build the categories (or histograms) to perform the division for fake rate calculation
// 3. holds the binned pdf for signal/background allowing for the calculation of correlated likelihood variables 
// 3. ouputs an xml/JSON containing 'fixed' fake-rates to be used in the probability calculation
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "jetSelector.h"

class globalJetProbabilities {
 public:
  // constructor for a single probability object
  explicit globalJetProbabilities(const std::string &label, 
				  const std::string &stack,
				  const bool & isMC_ , 
				  const bool & isSig_ , 
				  const float& evWeight_, 
				  TTree*& tree,
				  const jetSelector& jetSel,
				  const int & debug_);
  
  // constructor from a stored json 

  explicit globalJetProbabilities(const std::string& label, 
				  const std::string& stack,
				  const bool & isMC_ , 
				  const bool & isSig_ , 
				  const float& evWeight_, 
				  Json::Value probabilities,
				  const int & debug_);

  // add jets to the probabilities
  // void addJet(const jetCandidate&);  

  // when all the jets have been added, perform the division for the efficiency
  //  void generateGlobalProbabilities();

  // create the pdfs for the likelihood
  //void generateLikelihoodPDF();

  // called by jetTagAnalyzer to retrieve probabilities to be included in the outputted JSON
  //void getGlobalProbabilities();
  float getNJetTagProbabilityVector();

  // calculates the local pdf likelihood
  //float getJetLikelihood(const jetCandidate&);

  float getJetFakeProbability(float binVariable, float catVar);
  std::pair<float, float> getJetFakeProbabilityError(float binVariable, float catVar);
  void printHistStatus();

  Json::Value getProbabilitiesJSON();

  // acessor
  std::string getBinningVarName() { return binningVar; }
  std::string getCategoryVarName() { return categoryVar; }

  // accessors
  TGraphAsymmErrors getRatioGraph();
  TH1F		    getTaggedHist();
  TH1F		    getAllHist();
  TH1F              getCentralEffHist();
  TH1F              getCentralEffHistErrUp();
  TH1F              getCentralEffHistErrDn();

  // local variables
  TH1F		    taggedJetHist;
  TH1F		    allJetHist;
  TGraphAsymmErrors ratioGraph;
  // histograms of the assymetric errors
  TH1F              ratioHistEff;
  TH1F              ratioHistEffErrUp;
  TH1F              ratioHistEffErrDn;

  std::vector<double>           histBinVals;
  std::vector<double>           catBinVals;

  // bins and category variables
  std::string   binningVar;
  std::string   categoryVar;
  // selection strings
  std::string	jetCutString;
  std::string	eventCutString;
  std::string	baselineJetCutString;

  // configuration params
  const bool isMC, isSig;
  const float evWeight;
  const int debug;
  int nBins;
};
