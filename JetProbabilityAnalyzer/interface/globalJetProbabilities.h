// Joshua Hardenbrook 18-11-2015
// 1. encapsulates the relative tag efficiency given some selection
// 2. build the categories (or histograms) to perform the division fo
// 3. holds the binned pdf for signal/background allowing for the calculation of correlated likelihood variables 
// 3. ouputs an xml/JSON containing 'fixed' fake-rates to be used in the probability calculation


class globalJetProbabilities {

 public:
  // constructor
  globalJetProbailities(const jetSelection&);

  // re-build the jet probabilities from the outputted format 
  void loadConfiguration();

  // add jets to the probabilities
  void addJet(jetCandidate);  

  // when all the jets have been added, perform the division for the efficiency
  void generateGlobalProbabilities();
  // create the pdfs for the likelihood
  void generateLikelihoodPDF();

  // called by jetTagAnalyzer to retrieve probabilities to be included in the outputted JSON
  void getGlobalProbabilities();
  float getNJetTagProbability(vector<jetCandidates>, const int& ntags); 

  // calculates the local pdf likelihood
  float getJetLikelihood(const jetCandidate&);
}